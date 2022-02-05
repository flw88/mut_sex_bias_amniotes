#!/usr/bin/env Rscript
suppressMessages( library(getopt) )
suppressMessages( library(ape) )
suppressMessages( library(ggplot2) )
suppressMessages( library(data.table) )
suppressMessages( library(caper) )
suppressMessages( library(stringr) )
suppressMessages( library(Cairo) )
suppressMessages( library(factoextra) )
suppressMessages( library(ggrepel) )
suppressMessages( library(ungeviz) )
suppressMessages( library(castor) )
suppressMessages( library(geiger) )
suppressMessages( library(RColorBrewer) )

# Prepare data frame for PGLS
prepareDataframe <- function(d, xvar, yvar, xlog, ylog){
  v <- c("Species", xvar, yvar)
  subd <- d[,..v]
  colnames(subd) <- c("Species", "xvar", "yvar")
  subd <- as.data.frame(subd)
  if (xlog==TRUE){
    subd[,"xvar"] = log10(subd[,"xvar"])
  }
  if (ylog==TRUE){
    subd[,"yvar"] = log10(subd[,"yvar"])
  }
  return(subd)
}

# PGLS
runPGLS <- function(d, xvar, yvar, xlog, ylog, phylogeny, lambda, kappa, delta){
  
  subd <- prepareDataframe(d, xvar, yvar, xlog, ylog)
  
  subd <- subd[!is.na(subd$xvar),]
  
  if (nrow(subd)<3){
    return("NA")
  }
  trimmed_phylogeny = keep.tip(phylogeny, subd$Species)
  cdat <- comparative.data(data = subd, phy = trimmed_phylogeny, names.col = "Species")#,vcv=TRUE, vcv.dim=3
  mod <- pgls(yvar ~ xvar, cdat, lambda = lambda, kappa = kappa, delta = delta)#, 
              #control = list(fnscale=-1, ndeps=1e-3))
  
  df_mod = list("model" = mod, "data" = subd)
  return(df_mod)
}

# PIC pairwise
# Function to produce PICs based on pairwise comparison of tips only
# (no ancestral node reconstruction; no overlaps).
# For N tips, gives N/2 PIC values (rounded down).
# Uses a minimum distance heuristic to iteratively choose pairs.
picPairwiseTips <- function(x, phylogeny, method="minimum"){
  n <- Ntip(phylogeny)
  if(method == "random"){
    tip1 <- sample(phylogeny$tip.label, floor(n/2), replace=FALSE)
    tip2 <- sample(setdiff(phylogeny$tip.label, tip1), floor(n/2), replace=FALSE)
    tip.pairs <- cbind(tip1, tip2)
  } else if(method == "minimum"){
    dist.mat <- dist.nodes(phylogeny)[1:n, 1:n]
    # Set diagonal to infinity
    diag(dist.mat) <- Inf
    rownames(dist.mat) <- colnames(dist.mat) <- phylogeny$tip.label
    tip.pairs <- matrix(as.character(NA), nrow=floor(n/2), ncol=2)
    for(i in 1:nrow(tip.pairs)){ # Iterate through and take minimum tip pair
      min.i <- which(dist.mat == min(dist.mat), arr.ind=TRUE)[1,]
      tip.pairs[i,] <- rownames(dist.mat)[min.i]
      mask <- !(rownames(dist.mat) %in% tip.pairs[i,])
      dist.mat <- dist.mat[mask, mask]
    }
    
  }
  y <- x[tip.pairs[,1]] - x[tip.pairs[,2]]
  names(y) <- apply(tip.pairs, 1, function(s) str_c(s, collapse="-"))
  return(y)
}

# PIC Correlation
picCor <- function(d, xvar, yvar, xlog, ylog, phylogeny, method="spearman", pic.pairwise=FALSE){
  subd <- prepareDataframe(d, xvar, yvar, xlog, ylog)
  subd <- subd[!is.na(subd$xvar),]
  
  if (nrow(subd)<3){
    return("NA")
  }
  row.names(subd) <- subd$Species
  X <- subd[,"xvar"]
  Y <- subd[,"yvar"]
  names(X) <- names(Y) <- row.names(subd)
  subtree <- keep.tip(phylogeny, row.names(subd))
  
  if(pic.pairwise){
    pic.X <- picPairwiseTips(X, subtree, method="minimum")
    pic.Y <- picPairwiseTips(Y, subtree, method="minimum")
  } else {
    pic.X <- pic(X, subtree)
    pic.Y <- pic(Y, subtree)
  }
  cor.res <- cor.test(pic.X, pic.Y, method=method)
  return(list("PIC1"=pic.X, "PIC2"=pic.Y, "cor"=cor.res))
}

# PIC Regression
picModel <- function(d, xvar, yvar, xlog, ylog, phylogeny){
  
  subd <- prepareDataframe(d, xvar, yvar, xlog, ylog)
  subd <- subd[!is.na(subd$xvar),]
  
  if (nrow(subd)<3){
    return("NA")
  }
  row.names(subd) <- subd$Species
  X <- subd[,"xvar"]
  Y <- subd[,"yvar"]
  names(X) <- names(Y) <- row.names(subd)
  subtree <- keep.tip(phylogeny, row.names(subd))
  pic.X <- pic(X, subtree)
  pic.Y <- pic(Y, subtree)
  mod = lm(pic.Y ~ 0 + pic.X)
  df = data.frame(pic.X, pic.Y)
  df_mod = list("model" = mod, "data" = df)
  return(df_mod)
}

# Extract r2 & pval
extract_p <- function(mod, i){
  # coef <- round(summary(mod)$coefficients[i,4], digits = 4)
  coef <- summary(mod)$coefficients[i,4]
  return(coef)
}

# Extract transformed variance-covariance matrix
extract_vcm <- function(mod){
  mx <- mod$Vt
  sp <- mod$data$phy$tip.label
  row.names(mx) <- sp
  colnames(mx) <- sp
  df <- melt(mx)
  return(df)
}


# Extract data for particular comparison(s)
GetExperimentData <- function(xy_data, experiment){
  out.tab <- data.table()
  for(i in experiment){
    out.tab <- rbindlist(list(out.tab, xy_data[[i]]))
  }
  return(out.tab)
}
