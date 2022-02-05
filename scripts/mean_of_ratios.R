#!/usr/bin/env Rscript
rm(list=ls())

##### LIBRARIES #####

suppressMessages( library(getopt) )
library(data.table)
library(stringr)
library(ggplot2)

theme_set(theme_bw() +
            theme(axis.text=element_text(size=12), panel.border=element_rect(size = 1.5)))

##### GET OPTS, RUN PARAMETERS #####
opt.spec <- matrix(c("help",         "h",  0,  "logical",   0, 
                     "exp.name",     "c",  1,  "character", 1,
                     "ref.species",  "r",  1,  "character", 1
), byrow=TRUE, ncol=5)
opt <- getopt(opt.spec[,1:4])
req.args <- opt.spec[as.logical(as.integer(opt.spec[,5])), 1]

main.dir <- "/moto/palab/projects/male_mutation_bias_XA/mut_sex_bias_amniotes"

###--- Run this commented block to test locally ---#
# main.dir <- "/Users/felix/mt_mp_lab/projects/male_mutation_bias_XA/mut_sex_bias_amniotes"
# opt <- list(exp.name="Mammals", ref.species="Homo_sapiens")
###--------------------------------------------#

scripts.dir <- str_interp("${main.dir}/scripts")
data.dir    <- str_interp("${main.dir}/data")
pgls.dir    <- str_interp("${main.dir}/pgls_files")

source(str_interp("${scripts.dir}/alpha_regression.funcs.R"))

if( !is.null(opt$help) || (length(opt) == 0) ) {
  cat(getopt(opt.spec, usage=TRUE))
  quit(status=1)
}

for(arg in req.args){
  if( is.null(opt[[arg]]) ){
    cat(getopt(opt.spec, usage=TRUE))
    stop(str_interp("${arg} is a required argument."))
  }
}


grp <- opt$exp.name
ref <- opt$ref.species

out.prefix <- str_interp("${grp}.${ref}")
make.plots <- FALSE

##### LOAD DATA ######
if(grp == "Mammals"){
  load(str_interp("${scripts.dir}/lm_res/${grp}.${ref}.RData"))
} else if(grp == "Birds"){
  master.list <- list()
  for(i in 1:6){
    load(str_interp("${scripts.dir}/lm_res/${grp}_g${i}.${ref}.RData")) 
    if(i == 1){
      for(mt in names(lm.list)){
        master.list[[mt]] <- list()
      }
    }
    for(mt in names(lm.list)){
      for(sp in names(lm.list[[mt]])){
        if((i != 6) && (sp == "Gallus_gallus")){ next }
        master.list[[mt]][[sp]] <- lm.list[[mt]][[sp]]
      }
    }
  }
  lm.list <- master.list
  print(length(lm.list$mod))
}

##### FUNCTIONS #####

AlphaFromSampledRatios <- function(d, n.boot=1000, is.zw=FALSE){
  xz.sub <- d[region == "XZ", sub_rate]
  a.sub <- d[region == "A", sub_rate]

  # Sample ratios n.boot times
  rat.samp <- sample(xz.sub, size=n.boot, replace=TRUE) / sample(a.sub, size=n.boot, replace=TRUE)

  xza.est <- c("xz_a" = mean(rat.samp))
  a.est <- c("alpha"=CalcAlpha(xza.est, is.zw=is.zw))

  xza.ci <- ConfInt(rat.samp)
  names(xza.ci) <- c("xz_a_lwr", "xz_a_upr")

  a.ci <- sort(CalcAlpha(xza.ci, is.zw=is.zw))
  names(a.ci) <- c("alpha_lwr", "alpha_upr")

  xza.rofm <- c("xz_a_RofM" = mean(xz.sub) / mean(a.sub))
  a.rofm <- c("alpha_RofM"=CalcAlpha(xza.rofm, is.zw=is.zw))

  out.res <- c(a.est, a.ci, xza.est, xza.ci, xza.rofm, a.rofm)
  return(out.res)
}

##### FUNCTIONS #####


rpt.gc <- 0.03
mt <- "mod"

alpha.ratio <- data.table()
for(gc.width in seq(0.02, 0.1, by=0.01)){
  if(gc.width == rpt.gc){
    colA <- "# X(Z) windows"; colB <- "# X(A) windows"
    cat(str_interp("GC band width = $[0.3f]{gc.width}\n"))
    cat(str_interp("$[20s]{colA}$[20s]{colB}\n"))
    cat(str_c(rep("-", 40), collapse=""), "\n", sep="")  
  }

  for(sp in names(lm.list$mod)){
    cur.dat <- lm.list[[mt]][[sp]]$data
    mean.gc <- cur.dat[region == "XZ", mean(gc)]
    gc.range <- mean.gc + c(-gc.width, gc.width)/2
    cur.dat <- cur.dat[(gc > gc.range[1]) & (gc < gc.range[2])]
    
    # Report
    if(gc.width == rpt.gc){
      colA <- cur.dat[, table(region)["XZ"]]
      colB <- cur.dat[, table(region)["A"]]
      cat(str_interp("$[20d]{colA}$[20d]{colB}\n"))
    }

    new.row <- as.list(AlphaFromSampledRatios(cur.dat, is.zw=is.zw))
    new.row[["species"]]  <- sp
    new.row[["mut_type"]] <- mt
    new.row[["gc_width"]] <- gc.width

    alpha.ratio <- rbindlist(list(alpha.ratio, as.data.table(new.row)))
  }
}

fwrite(alpha.ratio, file=str_interp("${scripts.dir}/alphas/${out.prefix}.mean_of_ratios.tsv"),
       sep="\t", quote=FALSE)

##### Compare R2 #####
cat("\n")
colA <- "gc_width"; colB <- "R^2"; colC <- "MAE"
cat(str_interp("$[10s]{colA}$[10s]{colB}$[10s]{colC}\n"))
cat(str_c(rep("-", 30), collapse=""), "\n", sep="")  
for(gc.width in sort(unique(alpha.ratio$gc_width))){
  keep.col <- c("alpha", "species", "mut_type")
  cur.dat <- merge(alpha.dat[mut_type %in% mt, keep.col, with=FALSE], 
                   alpha.ratio[(mut_type %in% mt) & (gc_width == gc.width), keep.col, with=FALSE], 
    by=c("species", "mut_type"), suffixes=c("_regr", "_MofR"))

  m <- lm(as.formula(str_interp("alpha_MofR ~ alpha_regr")), data=cur.dat)
  r2 <- summary(m)$r.squared
  mae <- cur.dat[, sum(abs(alpha_MofR - alpha_regr))/.N]

  cat(str_interp("$[10.2f]{gc.width}$[10.3f]{r2}$[10.3f]{mae}\n"))
}







if(make.plots){
  col.dat <- fread(str_interp("${data.dir}/data/plot_colors.tsv"))
  setkey(col.dat, Order)

##### Plot Ratio vs Regression #####
  alpha.dat <- fread(str_interp("alphas/${grp}.${ref}.LM.tsv"))

  for(plt.gc in c(0.02, 0.03, 0.06)){
    keep.col <- c("alpha", "alpha_lwr", "alpha_upr", "xz_a", "xz_a_lwr", "xz_a_upr", "species", "mut_type")
    plt.dat <- merge(alpha.dat[mut_type %in% mt, keep.col, with=FALSE], 
                     alpha.ratio[(mut_type %in% mt) & (round(gc_width, 2) == plt.gc), keep.col, with=FALSE], 
      by=c("species", "mut_type"), suffixes=c("_regr", "_MofR"))

    use.color <- col.dat[grp, Color]

    for(v in c("alpha", "xz_a")){
      x.str <- str_c(v, "_regr")
      y.str <- str_c(v, "_MofR")

      m <- lm(as.formula(str_interp("${y.str} ~ ${x.str}")), data=plt.dat)
      p.val <- summary(m)$coefficients[2,4]
      r2 <- summary(m)$r.squared
      lab.dat <- data.table(x=-Inf, y=Inf, lab=str_interp("\n p-value = $[0.1e]{p.val}\n r2 = $[0.2f]{r2}"))

      p <- ggplot(aes_string(x=x.str, y=y.str), data=plt.dat) + 
        geom_point(color=use.color) + 
        geom_smooth(method="lm", color=use.color) + 
        geom_abline(intercept=0, slope=1, color="gray", linetype=2) +
        geom_text(aes(x=x, y=y, hjust=0, vjust=1, label=lab), data=lab.dat)

      if(v == "xz_a"){
        p <- p + xlab("E(dX)/E(dA)") + ylab("E(dX/dA)")
      } else {
        p <- p + xlab(bquote(hat(alpha) ~ "(Mean of Ratios)")) + 
                 ylab(bquote(hat(alpha) ~ "(Ratio of Means)"))
      }

      ggsave(p, filename=str_interp("pdfs/${grp}.${v}.mean_of_ratios.gc_$[0.2f]{plt.gc}.pdf"), width=4, height=4)
    }
  }

##### Plot Ratio vs RofM #####
  for(plt.gc in c(0.02, 0.03, 0.06)){
    plt.dat <- alpha.ratio[(mut_type %in% mt) & (round(gc_width, 2) == plt.gc)]
    use.color <- col.dat[grp, Color]
    
    for(v in c("alpha", "xz_a")){
      x.str <- str_c(v, "_RofM")
      y.str <- str_c(v, "")

      m <- lm(as.formula(str_interp("${y.str} ~ ${x.str}")), data=plt.dat)
      p.val <- summary(m)$coefficients[2,4]
      r2 <- summary(m)$r.squared
      lab.dat <- data.table(x=-Inf, y=Inf, lab=str_interp("\n p-value = $[0.1e]{p.val}\n r2 = $[0.2f]{r2}"))

      p <- ggplot(aes_string(x=x.str, y=y.str), data=plt.dat) + 
        geom_point(color=use.color) + 
        geom_smooth(method="lm", color=use.color) + 
        geom_abline(intercept=0, slope=1, color="gray", linetype=2) +
        geom_text(aes(x=x, y=y, hjust=0, vjust=1, label=lab), data=lab.dat)

      if(v == "xz_a"){
        p <- p + xlab("E(dX)/E(dA)") + ylab("E(dX/dA)")
      } else {
        p <- p + xlab(bquote(hat(alpha) ~ "(Mean of Ratios)")) + 
                 ylab(bquote(hat(alpha) ~ "(Ratio of Means)"))
      }

      ggsave(p, filename=str_interp("pdfs/${grp}.${v}.mean_of_ratios.no_regr.gc_$[0.2f]{plt.gc}.pdf"), width=4, height=4)
    }
  }
}