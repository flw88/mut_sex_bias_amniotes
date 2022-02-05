#!/usr/bin/env Rscript
rm(list=ls())

##### LIBRARIES #####
suppressMessages( library(data.table) )
suppressMessages( library(stringr) )
suppressMessages( library(ggplot2) )
suppressMessages( library(ggnewscale) )
suppressMessages( library(RColorBrewer) )

# ggplot theme
theme_set(theme_bw() +
            theme(axis.text=element_text(size=12), panel.border=element_rect(size = 1.5)))

##### LOAD DATA #####
scripts.dir <- "/Users/felix/mt_mp_lab/projects/male_mutation_bias_XA/mut_sex_bias_amniotes/scripts"
load(str_interp("${scripts.dir}/lm_res/pipeline_fig_data2.RData"))
# lm.hum, samp.coef.hum, lm.horse, samp.coef.horse

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
  return(unname(alph))
}

##### PLOT FUNCTION #####
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
  out.rate <- exp(out.rate) / explan.vars[["size"]]
  
  return(unname(out.rate))
}

# Report Confidence Interval
ConfInt <- function(x, interval=0.95, return.list=FALSE){
  y <- quantile(x, c(1-interval, 1+interval)/2)
  names(y) <- c("lwr", "upr")
  if(return.list){ y <- as.list(y) }
  return(y)
}

PlotRegression <- function(lm.list, feat.var, samp.coef.list, facet.nrow=1, 
                           norm.subrate=TRUE, thin.autosome=1){
  MakePlotData <- function(lm.obj, feat.var, samp.coef){
    # Points to plot
    pts.dat <- as.data.table(lm.obj$model)
    if((thin.autosome < 1) && (thin.autosome > 0)){
      a.rows <- pts.dat[, .I[region == "A"]]
      x.rows <- pts.dat[, .I[region == "XZ"]]
      pts.dat <- pts.dat[c(sample(a.rows, size=round(length(a.rows) * thin.autosome), replace=FALSE), x.rows)]
    } else if(thin.autosome != 1) {
      stop("thin.autosome must be in (0,1]")
    }
    
    pts.dat[, size := exp( `offset(log(size))` )]
    pts.dat[, sub_rate := counts / size]
    pts.dat[, `offset(log(size))` := NULL]
    
    if("I(gc^2)" %in% names(pts.dat)){
      pts.dat[, `I(gc^2)` := NULL]
    }
    
    # Regression line to plot
    x.range <- pts.dat[, as.list(range(get(feat.var))), by=region]
    setnames(x.range, c("V1", "V2"), c("MIN", "MAX"))
    
    # 1000 points across range of feature
    line.dat <- x.range[, seq(MIN, MAX, length.out=1000), by=region]
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
    
    line.dat[, sub_rate := exp(predict(lm.obj, newdata=line.dat)) / size]
    keep.cols <- c(keep.cols, "size")
    
    tmp  <- apply(samp.coef, 1, function(x) PredictSubRate(lm.obj, line.dat[, keep.cols, with=FALSE], lm.coef=x))
    tmp2 <- apply(tmp, 1, ConfInt)
    line.dat[, lwr := tmp2["lwr",]]
    line.dat[, upr := tmp2["upr",]]
    
    # Normalize to mean sub rate
    if(norm.subrate){
      norm.factor <- 1 / mean(pts.dat$sub_rate)
      pts.dat[,  sub_rate := sub_rate * norm.factor]
      line.dat[, sub_rate := sub_rate * norm.factor]
      line.dat[, lwr      := lwr *      norm.factor]
      line.dat[, upr      := upr *      norm.factor]  
    }
    
    
    # # Generate text annotation
    # r2 <- summary(lm.obj)$r.squared # R-squared
    # p.val <- CalcModelPvalue(lm.obj) # P-value
    # if(IsPoisson(lm.obj)){
    #   label.annot <- ""
    # } else {
    #   label.annot <- str_interp("r2 = $[0.2f]{r2}\n")
    # }
    # label.annot <- str_c(label.annot, str_interp("pval = $[0.1e]{p.val}"))
    # if((length(const.var) > 0) && !(all(const.var == "size"))){
    #   label.annot <- str_c(label.annot, "\nRegression at:")
    #   for(i in const.var){
    #     if(i == "size"){ next }
    #     cur.val <- line.dat[[i]][1]
    #     label.annot <- str_c(label.annot, str_interp("\n  ${i} = $[0.3f]{cur.val}"))
    #   }
    # }
    # 
    # label.dat <- data.table(x=-Inf, y=Inf, label=label.annot)
    
    return(list(pts.dat, line.dat))
  }
  
  # Plot the regression
  hi.colors <- c("#524C98", "#B04C00") #brewer.pal(n=3, "Dark2")[c(3,2)] # Purple, orange
  lo.colors <- c("#CFCDEB", "#FF9F55")
  pts.dat <- data.table(); line.dat <- data.table(); #label.dat <- data.table()
  nm <- names(lm.list)[1]
  # for(nm in names(lm.list)){
    plt.dat <- MakePlotData(lm.list[[nm]], feat.var, samp.coef.list[[nm]])
    
    # Add species name column
    # for(i in seq_along(plt.dat)){ plt.dat[[i]][, species := nm] }
    
    # Append
    pts.dat   <- rbindlist(list(pts.dat   , plt.dat[[1]]))
    line.dat  <- rbindlist(list(line.dat  , plt.dat[[2]]))
    # label.dat <- rbindlist(list(label.dat , plt.dat[[3]]))
  # }
  # gg.obj <- ggplot(pts.dat, aes_string(x=feat.var, y="sub_rate", color="region")) +
  #   geom_point(aes_string(alpha=feat.var), size=1, shape=16) + 
  #   geom_ribbon(aes_string(x=feat.var, ymin="lwr", ymax="upr", fill="region", color=NULL), data=line.dat, alpha=0.5) +
  #   geom_line(aes_string(x=feat.var, y="sub_rate", color="region"), data=line.dat) + 
  #   # geom_label(aes(x=x, y=y, label=label), data=label.dat, hjust=0, vjust=1, color="black") +
  #   ylab("Normalized substitution rate") +
  #   scale_color_manual(values=hi.colors) + scale_fill_manual(values=hi.colors) 
  #   # facet_wrap(vars(species), nrow=facet.nrow)
  
  # new.col <- str_interp("${feat.var}_binned")
  # pts.dat[, eval(new.col) := cut_interval(get(feat.var), n=9)]
  feat.range <- range(pts.dat[[feat.var]])
  # feat.range[1] <- feat.range[1] * 0.95
  
  names(hi.colors) <- c("A", "XZ")
  names(lo.colors) <- names(hi.colors)
  gg.obj <- ggplot() 
  for(i in names(hi.colors)){
    gg.obj <- gg.obj + geom_point(aes_string(x=feat.var, y="sub_rate", fill=feat.var), data=pts.dat[region==i], size=2.7, shape=21, color="white", stroke=0.375) +
      scale_fill_gradient(low=lo.colors[i], high=hi.colors[i], limits=feat.range)
    if(i == "A"){
      gg.obj <- gg.obj + new_scale_fill()#scale_color_brewer(palette="Purples")
    } #else {
      # gg.obj <- gg.obj + #scale_color_brewer(palette="Oranges")
    # }
  }
  for(i in names(hi.colors)){
    gg.obj <- gg.obj + geom_ribbon(aes_string(x=feat.var, ymin="lwr", ymax="upr", color=NULL), data=line.dat[region==i], fill="grey", alpha=0.5) +
      geom_line(aes_string(x=feat.var, y="sub_rate"), data=line.dat[region==i], color=hi.colors[i], size=1) 
  }
  if(norm.subrate){
    ylab.txt <- "Normalized substitution rate"
  } else {
    ylab.txt <- "Substitution rate"
  }
  gg.obj <- gg.obj + ylab(ylab.txt)
  
    return(gg.obj)
}

##### CREATE X/A PLOT #####
c.red <- "#FF3939"

gc.mean <- lm.horse[[1]]$data[region=="XZ", mean(gc)]
tmp <- predict(lm.horse[[1]], newdata=data.table(region=c("XZ", "A"), gc=c(gc.mean, gc.mean), size=c(100, 100)))
x.mean <- exp(tmp[1]) / 100
a.mean <- exp(tmp[2]) / 100
xa.mean <- x.mean / a.mean

gg.obj <- PlotRegression(lm.horse, "gc", samp.coef.horse, norm.subrate=FALSE, thin.autosome=0.5) + 
  geom_vline(xintercept=gc.mean, color=c.red, linetype=2) +
  annotate(geom="point", x=c(gc.mean, gc.mean), y=c(x.mean, a.mean), color=c.red, size=2) +
  xlab("GC Content") + 
  theme(legend.position="none", axis.text=element_blank(), axis.ticks=element_blank())
print(gg.obj)

ggsave(gg.obj, file=str_interp("${scripts.dir}/pdfs/pipeline_figure.1.pdf"),
       height=4, width=5, useDingbats=FALSE)

##### CREATE ALPHA PLOT #####
use.xlim <- c(0.75,1)

plt.dat <- data.table(x=seq(0.73, 1.2, length.out=150))
plt.dat[, y := CalcAlpha(x)]
gg.obj <- ggplot(plt.dat, aes(x=x, y=y)) + geom_line(color="#4D4D4D", size=1.2) +
  geom_vline(xintercept=xa.mean, color=c.red, linetype=2) + coord_cartesian(xlim=use.xlim, ylim=sort(CalcAlpha(use.xlim))) +
  annotate(geom="point", x=xa.mean, y=CalcAlpha(xa.mean), color=c.red, size=2) +
  xlab("X/A Substitution rate ratio") + ylab("Alpha")
print(gg.obj)
ggsave(gg.obj, file=str_interp("${scripts.dir}/pdfs/pipeline_figure.2.pdf"),
       height=4, width=4.2, useDingbats=FALSE)
