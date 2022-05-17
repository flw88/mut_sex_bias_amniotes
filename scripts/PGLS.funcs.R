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

pglsMod <- function (formula, data, lambda = 1, kappa = 1, delta = 1, param.CI = 0.95, 
                     control = list(fnscale = -1), bounds = NULL, optimPar = NULL) 
{
  Dfun <- function(Cmat) {
    iCmat <- solve(Cmat, tol = .Machine$double.eps)
    svdCmat <- La.svd(iCmat)
    D <- svdCmat$u %*% diag(sqrt(svdCmat$d)) %*% t(svdCmat$v)
    return(t(D))
  }
  if (!inherits(data, "comparative.data")) 
    stop("data is not a 'comparative' data object.")
  dname <- deparse(substitute(data))
  call <- match.call()
  miss <- model.frame(formula, data$data, na.action = na.pass)
  miss.na <- apply(miss, 1, function(X) (any(is.na(X))))
  if (any(miss.na)) {
    miss.names <- data$phy$tip.label[miss.na]
    data <- data[-which(miss.na), ]
  }
  m <- model.frame(formula, data$data)
  y <- m[, 1]
  x <- model.matrix(formula, m)
  k <- ncol(x)
  namey <- names(m)[1]
  xVar <- apply(x, 2, var)[-1]
  badCols <- xVar < .Machine$double.eps
  if (any(badCols)) 
    stop("Model matrix contains columns with zero variance: ", 
         paste(names(xVar)[badCols], collapse = ", "))
  if (is.null(data$vcv)) {
    V <- if (kappa == 1) {
      VCV.array(data$phy)
    }
    else {
      VCV.array(data$phy, dim = 3)
    }
    data$vcv <- V
  }
  else {
    V <- data$vcv
  }
  nm <- names(data$data)
  n <- nrow(data$data)
  if (!is.null(param.CI)) {
    if (!is.numeric(param.CI) || param.CI <= 0 || param.CI > 
        1) 
      stop("param.CI is not a number between 0 and 1.")
  }
  usrBounds <- bounds
  bounds <- list(kappa = c(1e-06, 3), lambda = c(1e-06, 1), 
                 delta = c(1e-06, 3))
  if (!is.null(usrBounds)) {
    if (!is.list(usrBounds)) 
      stop("Bounds must be a list of named bounds for any or all of kappa, lambda and delta")
    usrNames <- names(usrBounds)
    badNames <- setdiff(usrNames, c("kappa", "lambda", "delta"))
    if (length(badNames) > 0) 
      stop("The list of bounds contains names other than kappa, lambda and delta")
    for (nm in usrNames) {
      bounds[nm] <- usrBounds[nm]
    }
  }
  parVals <- list(kappa = kappa, lambda = lambda, delta = delta)
  for (i in seq_along(parVals)) {
    p <- parVals[[i]]
    nm <- names(parVals)[i]
    if (length(p) > 1) 
      stop(nm, " not of length one.")
    if (is.character(p) & p != "ML") 
      stop(nm, " is character and not 'ML'.")
    bnds <- bounds[[nm]]
    if (length(bnds) > 2) 
      stop("Bounds specified for ", nm, " not of length one.")
    if (!is.numeric(bnds)) 
      stop("Non-numeric bounds specified for ", nm, ".")
    if (any(bnds < 0)) 
      stop("Negative values in bounds specified for ", 
           nm, ".")
    lb <- bnds[1]
    ub <- bnds[2]
    if (lb > ub) 
      stop("Lower bound greater than upper bound for ", 
           nm, ".")
    if (is.numeric(p) & (p < lb | p > ub)) 
      stop(sprintf("%s value (%0.2f) is out of specified bounds [%0.2f, %0.2f]", 
                   nm, p, lb, ub))
  }
  if (kappa != 1 && length(dim(V)) != 3) 
    stop("3D VCV.array needed for kappa transformation.")
  mlVals <- sapply(parVals, "==", "ML")
  if (any(mlVals)) {
    parVals[mlVals] <- lapply(bounds, mean)[mlVals]
    parVals <- as.numeric(parVals)
    names(parVals) <- c("kappa", "lambda", "delta")
    if(is.null(optimPar)){
      optimPar <- parVals[mlVals]
    }
    fixedPar <- parVals[!mlVals]
    lower.b <- sapply(bounds, "[", 1)[mlVals]
    upper.b <- sapply(bounds, "[", 2)[mlVals]
    optim.param.vals <- optim(optimPar, fn = pgls.likelihood, 
                              method = "L-BFGS-B", control = control, upper = upper.b, 
                              lower = lower.b, V = V, y = y, x = x, fixedPar = fixedPar, 
                              optim.output = TRUE)
    if (optim.param.vals$convergence != "0") {
      stop("Problem with optim:", optim.param.vals$convergence, 
           optim.param.vals$message)
    }
    fixedPar <- c(optim.param.vals$par, fixedPar)
    fixedPar <- fixedPar[c("kappa", "lambda", "delta")]
  }
  else {
    fixedPar <- as.numeric(parVals)
    names(fixedPar) <- c("kappa", "lambda", "delta")
  }
  ll <- pgls.likelihood(optimPar = NULL, fixedPar = fixedPar, 
                        y, x, V, optim.output = FALSE)
  log.lik <- ll$ll
  Vt <- pgls.blenTransform(V, fixedPar)
  aic <- -2 * log.lik + 2 * k
  aicc <- -2 * log.lik + 2 * k + ((2 * k * (k + 1))/(n - k - 
                                                       1))
  coeffs <- ll$mu
  names(coeffs) <- colnames(x)
  varNames <- names(m)
  pred <- x %*% ll$mu
  res <- y - pred
  D <- Dfun(Vt)
  pres <- D %*% res
  fm <- list(coef = coeffs, aic = aic, log.lik = log.lik)
  RMS <- ll$s2
  RSSQ <- ll$s2 * (n - k)
  xdummy <- matrix(rep(1, length(y)))
  nullMod <- pgls.likelihood(optimPar = NULL, fixedPar = fixedPar, 
                             y, xdummy, V, optim.output = FALSE)
  NMS <- nullMod$s2
  NSSQ <- nullMod$s2 * (n - 1)
  errMat <- t(x) %*% solve(Vt) %*% x
  errMat <- solve(errMat) * RMS[1]
  sterr <- diag(errMat)
  sterr <- sqrt(sterr)
  RET <- list(model = fm, formula = formula, call = call, RMS = RMS, 
              NMS = NMS, NSSQ = NSSQ[1], RSSQ = RSSQ[1], aic = aic, 
              aicc = aicc, n = n, k = k, sterr = sterr, fitted = pred, 
              residuals = res, phyres = pres, x = x, data = data, varNames = varNames, 
              y = y, param = fixedPar, mlVals = mlVals, namey = namey, 
              bounds = bounds, Vt = Vt, dname = dname)
  class(RET) <- "pgls"
  if (any(miss.na)) {
    RET$na.action <- structure(which(miss.na), class = "omit", 
                               .Names = miss.names)
  }
  if (!is.null(param.CI) && any(mlVals)) {
    param.CI.list <- list(kappa = NULL, lambda = NULL, delta = NULL)
    mlNames <- names(mlVals)[which(mlVals)]
    for (param in mlNames) {
      param.CI.list[[param]] <- pgls.confint(RET, param, 
                                             param.CI)
    }
    RET$param.CI <- param.CI.list
  }
  return(RET)
}

pglsLambdaML <- function(fmula, cdat, init.lam.values=(1:9)/10, kappa=1, delta=1){
  max.lhood <- -Inf
  for(init.lam in init.lam.values){
    cur.ml <- try(pglsMod(fmula, cdat, lambda = "ML", kappa = kappa, delta = delta,
                      optimPar = c("lambda"=init.lam)))
    if(class(cur.ml) == "try-error"){ next }
    if(logLik(cur.ml)[1] > max.lhood){
      max.lhood <- logLik(cur.ml)[1]
      out.ml <- cur.ml
    }
  }
  
  return(out.ml)
}

# PGLS wrapper
runPGLS <- function(d, xvar, yvar, xlog, ylog, phylogeny, lambda, kappa, delta){
  
  subd <- prepareDataframe(d, xvar, yvar, xlog, ylog)
  
  subd <- subd[!is.na(subd$xvar),]
  
  if (nrow(subd)<3){
    return("NA")
  }
  trimmed_phylogeny = keep.tip(phylogeny, subd$Species)
  cdat <- comparative.data(data = subd, phy = trimmed_phylogeny, names.col = "Species")#,vcv=TRUE, vcv.dim=3
  
  if((lambda == "ML") && (kappa != "ML") && (delta != "ML")){
    mod <- pglsLambdaML(yvar ~ xvar, cdat, kappa=kappa, delta=delta)
  } else {
    mod <- try(pglsMod(yvar ~ xvar, cdat, lambda = lambda, kappa = kappa, delta = delta, optimPar = optimPar))
  }
  if(class(mod) == "try-error"){
    return(NA)
  }
  
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
