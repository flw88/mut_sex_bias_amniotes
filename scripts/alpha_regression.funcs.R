#!/usr/bin/env Rscript

##### LIBRARIES #####

require(R.utils)
require(data.table)
require(stringr)
require(ape)
require(ggplot2)
require(RColorBrewer)
require(castor)

##### UTILITY FUNCTIONS #####
FilterAndReport <- function(in.dat, rm.mask, filter.msg="Filtering:\n"){
  keep.count <- in.dat[!rm.mask, table(factor(region, levels=c("A","XZ","YW")))]
  cat(filter.msg)
  cat(str_interp("$[5d]{keep.count['A']} A   regions remaining\n"))
  cat(str_interp("$[5d]{keep.count['XZ']} X/Z regions remaining\n"))
  cat(str_interp("$[5d]{keep.count['YW']} Y/W regions remaining\n"))
  
  in.dat <- in.dat[!rm.mask]
  # cat(str_interp("${in.dat[,.N]} regions remaining\n"))
  return(in.dat)
}

# Map Strong/Weak single nt mutation types to regular types
SWtoType <- function(x){
  if(x == "S>S"){
    y <- "C>G"
  } else if(x == "W>W"){
    y <- "T>A"
  } else if(x == "S>W"){
    y <- c("C>A", "C>T")
  } else if(x == "W>S"){
    y <- c("T>C", "T>G")
  }
  return(y)
}

BGCtoType <- function(x){
  if(x == "gc-BGC-affected"){
    y <- c("C>A", "C>T", "T>C", "T>G")
  } else if(x == "gc-BGC-unaffected"){
    y <- c("C>G", "T>A")
  }
  return(y)
}
# Combine a subset of mutation types (numerator, denominator, rate)
# num.list and den.list should be lists whose names are the new columns
# New columns will be <col> for rate, <col>_num for numerator,
# <col>_den for denominator
CombineMutTypes <- function(in.dat, col.names, num.list, den.list){
  if( length(unique(length(col.names), length(num.list), length(den.list))) != 1){
    stop("col.names, num.list, and den.list must all have the same length")
  }
  
  num.names <- str_c(col.names, "_num")
  den.names <- str_c(col.names, "_den")
  for(i in seq_along(col.names)){
    in.dat[, eval(num.names[i]) := rowSums(.SD[, num.list[[i]], with=FALSE])]
    in.dat[, eval(den.names[i]) := rowSums(.SD[, den.list[[i]], with=FALSE])]
    in.dat[, eval(col.names[i]) := get(num.names[i]) / get(den.names[i])]
  }
  
  return(in.dat)
}

WhichEdgeInTree <- function(phy, targ.edge){
  i <- which(apply(phy$edge, 1, function(x) all(x == targ.edge)))
  return(i)
}

GetDescendants <- function(phy, targ.node){
  return( get_subtree_at_node(phy, targ.node - length(phy$tip.label))$subtree$tip.label )
}

GetAncBranch <- function(phy, tip.labs){
  mrca.node <- get_mrca_of_set(phy, tip.labs)
  root.node <- find_root(phy)
  if(mrca.node == root.node){
    stop("No ancestral branch. MRCA is root.")
  }
  targ.edge <- nodepath(phy, mrca.node, root.node)[1:2]
  if(length(WhichEdgeInTree(phy, targ.edge)) == 0){
    targ.edge <- rev(targ.edge)
    if(length(WhichEdgeInTree(phy, targ.edge)) == 0){
      stop("Something went wrong. Failed to find ancestral branch.")
    }
  }
  return(targ.edge)
}

DistTips2Parent <- function(phy){
  dist.mat <- dist.nodes(phy) # Node distance matrix
  path.list <- nodepath(phy) # List of root to tip node paths
  
  # Get distance from tip to parent
  out.vals <- lapply(path.list, function(x){y <- tail(x, n=2); dist.mat[y[1], y[2]]})
  
  names(out.vals) <- phy$tip.label
  return(out.vals)
}

# Parse newick string into tree object
# Returns list of sub rates at tips
ParseTree <- function(tree.str, edge.index=NULL){
  tree.obj <- read.tree(text=tree.str)
  if(is.null(edge.index)){
    out.list <- DistTips2Parent(tree.obj)
  } else {
    out.list <- as.list(tree.obj$edge.length[edge.index])
    names(out.list) <- str_c("E", edge.index)
  }
  return(out.list) # Returns a list
}

# Parse comma-separated list of variances and apply to edges of phy 
# Returns list of variances at tips 
ParseTreeVariance <- function(var.str, phy){
  phy$edge.length <- as.numeric(str_split(var.str, pattern=",")[[1]])
  return(DistTips2Parent(phy))
}

# Report Confidence Interval
ConfInt <- function(x, interval=0.95, return.list=FALSE){
  y <- quantile(x, c(1-interval, 1+interval)/2)
  names(y) <- c("lwr", "upr")
  if(return.list){ y <- as.list(y) }
  return(y)
}

# Check if model is weighted
IsWeighted <- function(lm.obj){
  return("(weights)" %in% names(lm.obj$model))
}

# Check if model is Poisson
IsPoisson <- function(lm.obj){
  out.bool <- ("glm" %in% class(lm.obj))
  if(out.bool){
    out.bool <- out.bool && (lm.obj$family$family == "poisson")
  }
  return(out.bool)
}

# Return boolean array of "odd" windows, ordered by chromosome size
# Negate to get "even" values.
OddWindows <- function(in.dat){
  return( in.dat[, .I %% 2 == 1, by=region]$V1 )
}

##### REGRESSION FUNCTIONS #####
CalcModelPvalue <- function(lm.obj){
  if(IsPoisson(lm.obj)){
    p.val <- pchisq(lm.obj$deviance, df=lm.obj$df.residual, lower.tail=FALSE)
  } else {
    fstat <- summary(lm.obj)$fstatistic
    p.val <- unname(pf(fstat[1], fstat[2], fstat[3], lower.tail=FALSE))
  }
  return(p.val)
}

# If poisson.glm==FALSE, response.var denotes the column containing sub rates. 
#   If TRUE, then response.var is a prefix and the function looks for 
#   <response.var>_num and <response.var>_den columns
# Set intercept.only=TRUE to run regression with a region effect (X/A) on
#   only the intercept. 
# The arg feat.vars should be a named char vector. These are the features to 
#   to control for (explanatory vars). The names of feat.vars are the names used
#   for the vars in the model formula
# Set variance.var to the name of the variance variable to perform weighted 
# regression. (Weights using 1/variance.var) 
# Default is to perform unweighted regresion allowing for full 
# interaction.
RegressRates <- function(in.dat, response.var, feat.vars, variance.var=NULL, 
                         intercept.only=FALSE, unrestrict.features=FALSE, 
                         gc.squared=FALSE, poisson.glm=FALSE, sex.chrom="XZ"){
  # Build regression data table
  cur.col <- unname(c(response.var, feat.vars, "region"))
  if(poisson.glm){
    count.cols <- str_c(response.var, c("_num", "_den"))
    cur.col <- c(cur.col, count.cols)
  }
  reg.dat <- in.dat[, cur.col, with=FALSE]
  reg.dat <- reg.dat[region %in% c("A", sex.chrom)]
  reg.dat[, region := factor(region, levels=c("A", sex.chrom))]
  
  # Set sub rate name
  setnames(reg.dat, response.var, "sub_rate")
  
  # Set feature variables names
  setnames(reg.dat, feat.vars, names(feat.vars))
  
  # Filter feature variables
  for(nm in names(feat.vars)){
    reg.dat <- FilterAndReport(reg.dat, reg.dat[, is.na(get(nm))], str_interp("Filtering NaN ${nm}:\n"))
  }
  
  # Set count names
  if(poisson.glm){
    setnames(reg.dat, count.cols, c("counts", "size"))
    reg.dat[, counts := round(counts)]
    reg.dat[,   size := round(size)]
  }
  
  # Set variance var
  if(!is.null(variance.var)){
    reg.dat[, w := 1 / in.dat[[variance.var]]]
  }
  
  # Filter to only include overlapping feature range of X and A
  if(!unrestrict.features){
    for(i in names(feat.vars)){
      cur.min <- max(reg.dat[, min(get(i)), by=region]$V1)
      cur.max <- min(reg.dat[, max(get(i)), by=region]$V1)
      rm.mask <- reg.dat[, (get(i) < cur.min) | (get(i) > cur.max)]
      reg.dat <- FilterAndReport(reg.dat, rm.mask, str_interp("Filtering down to overlapping ${i} range of X and A:\n"))
    }
  }
  
  # Create formula
  if(poisson.glm){
    f.str <- "counts ~ offset(log(size)) + "
  } else {
    f.str <- "sub_rate ~ "
  }
  
  if(intercept.only){
    f.str <- str_c(f.str, str_c(c(names(feat.vars), "region"), collapse=" + "))
  } else { # Full model with interaction
    f.str <- str_c(f.str, str_c(names(feat.vars), "region", sep="*", collapse=" + "))
  }
  
  if(gc.squared){
    if(!intercept.only){
      stop("Can only include gc^2 term without full interaction")
    }
    f.str <- str_c(f.str, "+ I(gc^2)")
  }
  
  fmula <- as.formula(f.str)
  
  # Run regression
  if(!is.null(variance.var)){ # weighted
    lm.obj <- lm(fmula, data=reg.dat, weights=w)
  } else { # not weighted
    if(poisson.glm){
      # reg.dat[, counts := round(sub_rate * size)]
      lm.obj <- glm(fmula, data=reg.dat, family=poisson(link=log))
    } else {
      lm.obj <-  lm(fmula, data=reg.dat)
    }
  }
  return(lm.obj)
}

# Compute predicted sub rate from list of explanatory variables
#   using given linear model or regression model coefficients (if provided)
PredictSubRate <- function(lm.obj, explan.vars, lm.coef=NULL){
  tmp.lm <- lm.obj
  if(!is.null(lm.coef)){ # Use the provided coefficients rather than lm.obj
    if(!all(names(tmp.lm$coefficients) == names(lm.coef))){
      stop("Alternative coefficients 'lm.coef' do not match coefficients in model 'lm.obj'")
    }
    tmp.lm$coefficients <- lm.coef
    # Note that we can do the above mainly because predict() just calls the coefficients
    # from the LM object and applies the formula. I think as long as the names are correct
    # this should always work, though it's a bit hacky.
  }
  out.rate <- predict(tmp.lm, newdata=explan.vars)
  if(IsPoisson(lm.obj)){
    out.rate <- exp(out.rate) / explan.vars[["size"]]
  }
  
  return(unname(out.rate))
}

# Repeatedly resample points and recompute coefficients. Returns table of coefficients.
ResampleRegression <- function(lm.obj, n.boot=500, by.region=FALSE){
  SampRegr <- function(lm.obj, by.region){
    if(by.region){ # If by region, then sample X and A separately
      samp.rows <- c(sample(which(lm.obj$model$region == "XZ"), replace=TRUE),
                     sample(which(lm.obj$model$region == "A"), replace=TRUE))
    } else {
      samp.rows <- sample.int(nrow(lm.obj$model), replace=TRUE)
    }
    
    use.dat <- lm.obj$model[samp.rows,]
    
    if(IsWeighted(lm.obj)){ # Weighted regression
      resamp.lm <- lm(lm.obj$terms, data=use.dat, weights=`(weights)`)
    } else {
      # If Poisson
      if(IsPoisson(lm.obj)){
        use.dat[["size"]] <- exp(use.dat[["offset(log(size))"]])
        resamp.lm <- glm(lm.obj$terms, data=use.dat, family=poisson(link=log))
      } else {
        resamp.lm <- lm(lm.obj$terms, data=use.dat)
      }
      
    }
    return(resamp.lm$coefficients)
  }
  
  samp.coef <- t(replicate(n.boot, expr=SampRegr(lm.obj, by.region)))
  return(samp.coef)
}

# Print summary of regression results
CatRegressionSummary <- function(lm.list, out.file){
  sink(out.file)
  for(nm in names(lm.list)){
    hash.str.1 <- "#####################################################"
    hash.str.2 <- str_c(rep("#", str_length(hash.str.1) - str_length(nm) - 14), collapse="")
    cat(str_interp("${hash.str.1}\n"))
    cat(str_interp("##########  ${nm}  ${hash.str.2}\n"))
    print(summary(lm.list[[nm]]))
    
    # Proportion variance (or deviance if Poisson) explained
    a <- anova(lm.list[[nm]])
    if(IsPoisson(lm.list[[nm]])){
      print(a)
    } else{
      tmp.dat <- data.frame("Percent.Var.Explained"=a[["Sum Sq"]] / sum(a[["Sum Sq"]]) * 100)
      rownames(tmp.dat) <- rownames(a)
      print(tmp.dat)
    }
    cat("\n")
  }
  sink()
}


##### ALPHA FUNCTIONS #####
Miyata <- function(x, is.zw=FALSE){
  if(is.zw){
    y <- (4*x + 2) / (3*x + 3)
  } else {
    y <- (2*x + 4) / (3*x + 3)
  }
  return(y)
}

CalcAlpha <- function(xza.rat, g.rat=1, is.zw=FALSE){
  miy.g <- Miyata(g.rat, is.zw)
  if(is.zw){
    alph <- ((3 * xza.rat * miy.g) - 2) /
            (4 - (3 * xza.rat * miy.g))
  } else {
    alph <- ((3 * xza.rat * miy.g) - 4) /
            (2 - (3 * xza.rat * miy.g))
  }
  alph[alph < 0] <- Inf
  
  return(unname(alph))
}

# Compute alpha from regression using mean feature values (default)
#   or provided feature values (feat.values).
#   Returns alpha (alpha) w/ 95% CIs (alpha_lwr, alpha_upr) and
#   X(Z)-A ratios (xz_a)  w/ 95% CIs (xz_a_lwr, xz_a_upr)
AlphaFromRegression <- function(lm.obj, feat.vars, samp.coef=NULL, feat.values=NULL, g.rat=1, is.zw=FALSE){
  if(is.null(feat.values)){ # Get mean feature values across X(Z)
    x.rows <- lm.obj$model$region == "XZ"
    feat.values <- list()
    for(i in feat.vars){
      if(IsWeighted(lm.obj)){
        feat.values[[i]] <- rep(weighted.mean(lm.obj$model[x.rows, i], w=lm.obj$model[x.rows, "(weights)"]), 2)
      } else {
        feat.values[[i]] <- rep(mean(lm.obj$model[x.rows, i]), 2)
      }
    }
    
    if(IsPoisson(lm.obj)){
      feat.values[["size"]] <- rep(exp(mean(lm.obj$model[x.rows, "offset(log(size))"])), 2)
    } 
    
    feat.values[["region"]] <- c("A","XZ")
    feat.values <- as.data.frame(feat.values)
  } else { # Use provided feat.values
    if( !all(c(feat.vars, "region") %in% names(feat.values)) ){
      stop("feat.values provided but names do not include feat.vars and/or 'region")
    }
    
    if( !(all(sort(feat.values$region) == c("A", "XZ"))) ){
      stop("feat.values$region must be c('A','XZ')")
    }
  }
  
  # Compute rates (point estimate)
  est.rates <- PredictSubRate(lm.obj, feat.values)
  names(est.rates) <- feat.values$region
  
  # Compute resampled rates
  if(!is.null(samp.coef)){
    resamp.rates <- apply(samp.coef, 1, function(x) PredictSubRate(lm.obj, feat.values, lm.coef=x))
    rownames(resamp.rates) <- feat.values$region
  }
  
  # Compute alpha
  xza.fit <- c("xz_a"=unname(est.rates["XZ"] / est.rates["A"]))
  xza.ci <- c()
  if(!is.null(samp.coef)){
    xza.ci <- ConfInt(resamp.rates["XZ",] / resamp.rates["A",])
    names(xza.ci) <- str_c("xz_a_", names(xza.ci))
  }
  
  a.fit <- c("alpha"=CalcAlpha(xza.fit, is.zw=is.zw))
  a.ci <- c()
  if(!is.null(samp.coef)){
    a.ci <- sort(CalcAlpha(xza.ci, is.zw=is.zw))
    names(a.ci) <- c("alpha_lwr", "alpha_upr")
  }
  
  out.res <- c(a.fit, a.ci, xza.fit, xza.ci)
  return(out.res)
}


# Compute alpha just using mean across windows (default).
#   Can be run weighted.
#   Returns alpha (alpha) w/ 95% CIs (alpha_lwr, alpha_upr) and
#   X(Z)-A ratios (xz_a)  w/ 95% CIs (xz_a_lwr, xz_a_upr)
AlphaFromMean <- function(dat, subrate.col, is.zw=FALSE, weighted=FALSE, n.boot=500, sex.chrom="XZ"){
  # Build data table
  size.col <- str_c(subrate.col, "_den")
  keep.col <- unname(c(subrate.col, size.col, "region"))
  
  reg.dat <- dat[, keep.col, with=FALSE]
  reg.dat <- reg.dat[region %in% c("A", sex.chrom)]

  # Set sub rate name, size columns
  setnames(reg.dat, c(subrate.col, size.col), c("sub_rate", "size"))
  
  if(!weighted){
    reg.dat[, size := 1]
  }

  # Function for calculating X(Z)-A ratios
  CalcXZA <- function(in.dat, i=in.dat[,.I]){
    tmp <- in.dat[i, weighted.mean(sub_rate, w=size), by=region]
    xza.est <- tmp[region == "XZ", V1] / tmp[region == "A", V1]
    return(xza.est)
  }

  # Calculate point estimates
  xza.est <- c("xz_a"=CalcXZA(reg.dat))
  a.est <- c("alpha"=CalcAlpha(xza.est, is.zw=is.zw))

  # Calculate CIs
  xza.ci <- ConfInt( replicate(CalcXZA(reg.dat, i=reg.dat[, sample(.I, replace=TRUE), by=region]$V1), n=n.boot) )
  names(xza.ci) <- c("xz_a_lwr", "xz_a_upr")

  a.ci <- sort(CalcAlpha(xza.ci, is.zw=is.zw))
  names(a.ci) <- c("alpha_lwr", "alpha_upr")

  out.res <- c(a.est, a.ci, xza.est, xza.ci)
  return(out.res)
}

##### PLOT FUNCTIONS #####

# Plot the regression using feat.var as the variable, all other explanatory
# variables (aside from X/A) held constant at mean value 
PlotRegression <- function(lm.list, feat.var, samp.coef.list, facet.nrow=1, norm.subrate=TRUE){
  MakePlotData <- function(lm.obj, feat.var, samp.coef){
    # Points to plot
    pts.dat <- as.data.table(lm.obj$model)
    if(IsWeighted(lm.obj)){
      pts.dat[, `(weights)` := NULL]
    }
    if(IsPoisson(lm.obj)){
      pts.dat[, size := exp( `offset(log(size))` )]
      pts.dat[, sub_rate := counts / size]
      pts.dat[, `offset(log(size))` := NULL]
    }
    if("I(gc^2)" %in% names(pts.dat)){
      pts.dat[, `I(gc^2)` := NULL]
    }
    
    # Regression line to plot
    x.range <- pts.dat[, as.list(range(get(feat.var))), by=region]
    setnames(x.range, c("V1", "V2"), c("MIN", "MAX"))
    
    # 100 points across range of feature
    line.dat <- x.range[, seq(MIN, MAX, length.out=100), by=region]
    setnames(line.dat, "V1", feat.var)
    const.var <- names(pts.dat)[!(names(pts.dat) %in% c(feat.var, "region", "sub_rate", "counts"))]
    
    # Set other variables constant to mean
    for(i in const.var){
      tmp <- pts.dat[, mean(get(i))]
      line.dat[, eval(i) := tmp]
    }
    
    # Add confidence interval on regression
    # if(do.boot.ci){
    keep.cols <- c(feat.var, const.var, "region")
    if(IsPoisson(lm.obj)){
      line.dat[, sub_rate := exp(predict(lm.obj, newdata=line.dat)) / size]
      keep.cols <- c(keep.cols, "size")
    } else {
      line.dat[, sub_rate := predict(lm.obj, newdata=line.dat)]
    }
    tmp  <- apply(samp.coef, 1, function(x) PredictSubRate(lm.obj, line.dat[, keep.cols, with=FALSE], lm.coef=x))
    tmp2 <- apply(tmp, 1, ConfInt)
    line.dat[, lwr := tmp2["lwr",]]
    line.dat[, upr := tmp2["upr",]]
    
    ### Alternate code for using the regression confidence intervals rather than bootstrap reps
    # } else {
    #   pred.res <- predict(lm.obj, newdata=line.dat, interval="confidence")
    #   for(i in colnames(pred.res)){
    #     line.dat[, eval(i) := pred.res[,i]]
    #   }
    #   setnames(line.dat, "fit", "sub_rate")
    # }
    ####################################### --
    
    # Normalize to mean sub rate
    if(norm.subrate){
      norm.factor <- 1 / mean(pts.dat$sub_rate)
      pts.dat[,  sub_rate := sub_rate * norm.factor]
      line.dat[, sub_rate := sub_rate * norm.factor]
      line.dat[, lwr      := lwr *      norm.factor]
      line.dat[, upr      := upr *      norm.factor]  
    }
    
    
    # Generate text annotation
    r2 <- summary(lm.obj)$r.squared # R-squared
    p.val <- CalcModelPvalue(lm.obj) # P-value
    if(IsPoisson(lm.obj)){
      label.annot <- ""
    } else {
      label.annot <- str_interp("r2 = $[0.2f]{r2}\n")
    }
    label.annot <- str_c(label.annot, str_interp("pval = $[0.1e]{p.val}"))
    if((length(const.var) > 0) && !(all(const.var == "size"))){
      label.annot <- str_c(label.annot, "\nRegression at:")
      for(i in const.var){
        if(i == "size"){ next }
        cur.val <- line.dat[[i]][1]
        label.annot <- str_c(label.annot, str_interp("\n  ${i} = $[0.3f]{cur.val}"))
      }
    }
    
    label.dat <- data.table(x=-Inf, y=Inf, label=label.annot)
    
    return(list(pts.dat, line.dat, label.dat))
  }
  
  # Plot the regression
  use.colors <- brewer.pal(n=3, "Dark2")[1:2]
  pts.dat <- data.table(); line.dat <- data.table(); label.dat <- data.table()
  for(nm in names(lm.list)){
    plt.dat <- MakePlotData(lm.list[[nm]], feat.var, samp.coef.list[[nm]])
    
    # Add species name column
    for(i in seq_along(plt.dat)){ plt.dat[[i]][, species := nm] }
    
    # Append
    pts.dat   <- rbindlist(list(pts.dat   , plt.dat[[1]]))
    line.dat  <- rbindlist(list(line.dat  , plt.dat[[2]]))
    label.dat <- rbindlist(list(label.dat , plt.dat[[3]]))
  }
  gg.obj <- ggplot(pts.dat, aes_string(x=feat.var, y="sub_rate", color="region")) +
    geom_point(alpha=0.2, size=1) + 
    geom_ribbon(aes_string(x=feat.var, ymin="lwr", ymax="upr", fill="region", color=NULL), data=line.dat, alpha=0.5) +
    geom_line(aes_string(x=feat.var, y="sub_rate", color="region"), data=line.dat) + 
    geom_label(aes(x=x, y=y, label=label), data=label.dat, hjust=0, vjust=1, color="black") +
    ylab("Normalized substitution rate") +
    scale_color_manual(values=use.colors) + scale_fill_manual(values=use.colors) + 
    facet_wrap(vars(species), nrow=facet.nrow)
  return(gg.obj)
}

# Plot alpha estimate across the given feature variable
PlotAlphaByFeat <- function(lm.list, feat.var, samp.coef.list=NULL, facet.nrow=1){
  MakePlotData <- function(lm.obj, feat.var, samp.coef){
    pts.dat <- as.data.table(lm.obj$model)
    if(IsWeighted(lm.obj)){
      pts.dat[, `(weights)` := NULL]
    }
    
    ### Alpha line to plot
    x.range <- pts.dat[, as.list(range(get(feat.var))), by=region]
    setnames(x.range, c("V1", "V2"), c("MIN", "MAX"))
    
    # 100 points across range of feature
    tmp <- x.range[region=="XZ", seq(MIN, MAX, length.out=100)]
    line.dat <- data.table()
    line.dat[, eval(feat.var) := tmp]
    const.var <- names(pts.dat)[!(names(pts.dat) %in% c(feat.var, "region", "sub_rate"))]
    
    # Set other variables constant to mean
    for(i in const.var){
      tmp <- pts.dat[, mean(get(i))]
      line.dat[, eval(i) := tmp]
    }
    rm(pts.dat)
    
    # Add confidence interval on alpha
    FeatValue <- function(x){
      out.dat <- data.frame(region=c("A", "XZ"))
      for(i in names(x)){ out.dat[[i]] <- rep(x[i], 2) }
      return(out.dat)
    }
    
    tmp  <- apply(line.dat, 1, function(x) AlphaFromRegression(lm.obj, names(x), samp.coef, feat.values=FeatValue(x)))
    for(nm in rownames(tmp)){
      line.dat[, eval(nm) := tmp[nm,]]
    }
    
    # Generate text annotation
    # r2 <- summary(lm.obj)$r.squared # R-squared
    # p.val <- CalcModelPvalue(lm.obj) # P-value
    # label.annot <- str_interp("r^2 = $[0.2f]{r2}")
    # label.annot <- str_c(label.annot, str_interp("\npval = $[0.1e]{p.val}"))
    # if(length(const.var) > 0){
    #   label.annot <- str_c(label.annot, "\nRegression at:")
    #   for(i in const.var){
    #     cur.val <- line.dat[[i]][1]
    #     label.annot <- str_c(label.annot, str_interp("\n  ${i} = $[0.5f]{cur.val}"))
    #   }
    # }
    # 
    # label.dat <- data.table(x=-Inf, y=Inf, label=label.annot)
    
    return(line.dat)
  }
  
  # Plot the regression
  use.color <- "#999999"
  line.dat <- data.table()
  for(nm in names(lm.list)){
    cat(nm, "\n", sep="")
    plt.dat <- MakePlotData(lm.list[[nm]], feat.var, samp.coef.list[[nm]])
    # Add species name column
    plt.dat[, species := nm]
    
    # Append
    line.dat  <- rbindlist(list(line.dat  , plt.dat))
  }
  gg.obj <- ggplot(line.dat, aes_string(x=feat.var, y="alpha", ymin="alpha_lwr", ymax="alpha_upr")) +
    geom_ribbon(fill=use.color, alpha=0.5) +
    geom_line(color=use.color) + xlab("Alpha") +
    facet_wrap(vars(species), nrow=facet.nrow)
  return(gg.obj)
}

