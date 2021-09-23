#!/usr/bin/env Rscript
rm(list=ls())

##### LIBRARIES #####

suppressMessages( library(getopt) )
suppressMessages( library(R.utils) )
suppressMessages( library(data.table) )
suppressMessages( library(stringr) )
suppressMessages( library(ape) )
suppressMessages( library(ggplot2) )
suppressMessages( library(RColorBrewer) )

##### GET OPTS, RUN PARAMETERS #####
# ggplot theme
theme_set(theme_bw() +
            theme(axis.text=element_text(size=12), panel.border=element_rect(size = 1.5)))

### Get options      varName              flag  takesArg  type       isRequired         
opt.spec <- matrix(c("help",                "h",  0,     "logical",   0, 
                     "ref.species",         "r",  1,     "character", 1, 
                     "exp.name",            "c",  1,     "character", 1,
                     "out.exp",             "o",  1,     "character", 1
), byrow=TRUE, ncol=5)
opt <- getopt(opt.spec[,1:4])
req.args <- opt.spec[as.logical(as.integer(opt.spec[,5])), 1]

scripts.dir <- "/moto/palab/projects/male_mutation_bias_XA/mut_sex_bias_amniotes/scripts"

###--- Run this commented block to test locally ---#
# scripts.dir <- "/Users/felix/mt_mp_lab/projects/male_mutation_bias_XA/mut_sex_bias_amniotes/scripts"
# opt <- list(ref.species="Gallus_gallus",
#          exp.name="Aves"
#      )
###--------------------------------------------#

if( !is.null(opt$help) || (length(opt) == 1) ) {
  cat(getopt(opt.spec, usage=TRUE))
  quit(status=1)
}

for(arg in req.args){
  if( is.null(opt[[arg]]) ){
    cat(getopt(opt.spec, usage=TRUE))
    stop(str_interp("${arg} is a required argument."))
  }
}

ref.species <- opt$ref.species
exp.name <- opt$exp.name
out.exp <- opt$out.exp # Output experiments. Sub groupings of a particular order
out.exp <- str_split(out.exp, ",")[[1]]

##### LOAD DATA #####
tmp.lm <- NULL
if(length(out.exp) > 1){
  for(oe in out.exp){
    cur.fn <- str_interp("${scripts.dir}/lm_res/${oe}.${ref.species}.RData")
    load(cur.fn)
    
    if(is.null(tmp.lm)){
      tmp.lm <- lm.list
      tmp.coef <- samp.coef.list
    } else {
      for(j in names(lm.list)){
        for(k in names(lm.list[[j]])){
          tmp.lm[[j]][[k]]   <- lm.list[[j]][[k]]
          tmp.coef[[j]][[k]] <- samp.coef.list[[j]][[k]]
        }
      }
    }
  }
  lm.list <- tmp.lm
  samp.coef.list <- tmp.coef
} else {
  load(str_interp("${scripts.dir}/lm_res/${exp.name}.${ref.species}.RData"))
}

# Get orderr information
order.dat <- data.table(Species = names(lm.list$mod), Order="all")
if(exp.name == "Mammals"){
  meta.dat <- fread(str_interp("${scripts.dir}/../data/zoonomia_assembly_metadata.csv")) # Information on orders
  order.dat <- rbindlist(list(order.dat,
                              meta.dat[Species %in% names(lm.list$mod), .(Species, Order)]))
} else if(length(out.exp) > 1){
  for(oe in out.exp){
    cur.sp <- str_split(readLines(str_interp("${scripts.dir}/../data/${oe}.txt")), ",")[[1]]
    order.dat <- rbindlist(list(order.dat, data.table(Species=cur.sp, Order=oe)))
  }
}

# Latin to common names conversion table
latin.dat <- fread(str_interp("${scripts.dir}/../data/latin2common.txt"))
setkey(latin.dat, Species)

# Plot colors
color.dat <- fread(str_interp("${scripts.dir}/../data/plot_colors.tsv"))
setkey(color.dat, Order)

##### FUNCTION FOR TESTING ALPHA ESTIMATE #####
TestAlphaConstancy <- function(lm.list, samp.coef.list, feat.vars, targ.rate.ratio, include.sp=NULL){
  # is.targ <- names(lm.list) == targ.species
  # if(!any(is.targ)){
  #   stop("Target species ${targ.species} not found")
  # }
  p.dat <- list()
  
  if(is.null(include.sp)){
    include.sp <- names(lm.list)
  }
  
  for(sp in include.sp){
    lm.obj <- lm.list[[sp]]
    samp.coef <- samp.coef.list[[sp]]
    
    x.rows <- lm.obj$model$region == "XZ"
    feat.values <- list()
    for(i in feat.vars){
      feat.values[[i]] <- rep(mean(lm.obj$model[x.rows, i]), 2)
    }
    if(IsPoisson(lm.obj)){
      feat.values[["size"]] <- rep(exp(mean(lm.obj$model[x.rows, "offset(log(size))"])), 2)
    } 
    feat.values[["region"]] <- c("A","XZ")
    feat.values <- as.data.frame(feat.values)
    
    resamp.rates <- apply(samp.coef, 1, function(x) PredictSubRate(lm.obj, feat.values, lm.coef=x))
    rownames(resamp.rates) <- feat.values$region
    
    resamp.ratios <- resamp.rates["XZ",] / resamp.rates["A",]
    p <- sum(targ.rate.ratio > resamp.ratios) / length(resamp.ratios)
    if(p > 0.5){ p <- 1 - p }
    
    p.dat[[sp]] <- data.table(species=sp, pval=p, xz_a_comparison=targ.rate.ratio)
  }
  p.dat <- rbindlist(p.dat)
  return(p.dat)
}

##### GET P VALUES FOR MEAN OF EACH ORDER #####
pval.dat <- list()
for(ord in unique(order.dat$Order)){
  cur.sp <- order.dat[Order == ord, Species]
  if(length(cur.sp) < 3){ next }
  
  targ.rat <- mean(alpha.dat[(mut_type == "mod") & (species %in% cur.sp), xz_a])
  pval.dat[[ord]] <- TestAlphaConstancy(lm.list$mod, samp.coef.list$mod, 
                                       feat.vars, targ.rat, 
                                       include.sp=cur.sp)
  pval.dat[[ord]][, order := ord]
}

##### CONVERT TO DATA.TABLE #####
pval.dat <- rbindlist(pval.dat)
pval.dat[, alpha_comparison := CalcAlpha(xz_a_comparison, is.zw=is.zw)]

##### WRITE OUT #####
fwrite(pval.dat, file=str_interp("${scripts.dir}/alphas/${exp.name}.${ref.species}.pvals_by_order.tsv"), sep="\t", quote=FALSE,
       col.names=TRUE, row.names=FALSE)

##### PLOT P-VALUES #####
# Function for making species label string
StringifySpecies <- function(sp.names, width=30){
  out.str <- str_wrap(str_c(sort(sp.names), collapse=", "), width=width)
  return(out.str)
}

if(exp.name %in% color.dat$Order){
  hist.col <- color.dat[exp.name, Color]
} else {
  hist.col <- "#606060"
}

pval.dat[, common := latin.dat[species, Common_names]] # Get common names

n.orders <- length(unique(pval.dat$order))
log.break <- floor(log10(min(pval.dat[pval != 0, pval]))) # Where to switch log axis back to linear
x.trans <- scales::pseudo_log_trans(sigma=10^log.break, base=10)

lab.dat <- pval.dat[, .(a  = sprintf("H0: alpha = %0.2f", alpha_comparison[1]), 
                        small_p_species = sprintf("Species with p=0:\n%s", StringifySpecies(common[pval==0])),
                        n_small_p = sum(pval == 0)), by=order]
lab.dat[(n_small_p == 0), small_p_species := "Species with p=0:\n[none]"]
lab.dat[, label := str_c(a, "\n", small_p_species)]

n.row <- floor(sqrt(n.orders))
n.col <- ceiling(n.orders / n.row)

# per order
if(length(unique(pval.dat$order)) > 1){
  p <- ggplot(aes(x=pval), data=pval.dat[order != "all"]) + 
    geom_histogram(fill=hist.col) +
    geom_label(aes(x = Inf, y = Inf, label=label), data=lab.dat[order != "all"], hjust=1, vjust=1, color="black") +
    scale_x_continuous(trans=x.trans, breaks=c(0, 10^c(log.break:0))) +
    facet_wrap(vars(order), nrow=n.row) +
    theme(axis.text.x = element_text(angle=45, hjust=1))
  
  ggsave(p, file=str_interp("${scripts.dir}/pdfs/${exp.name}.${ref.species}.pvals_by_order.pdf"), width=n.row*4, height=n.col*3)
}

# all
p <- ggplot(aes(x=pval), data=pval.dat[order == "all"]) + 
  geom_histogram(fill=hist.col) +
  geom_label(aes(x = Inf, y = Inf, label=label), data=lab.dat[order == "all"], hjust=1, vjust=1, color="black", size=2.5) +
  scale_x_continuous(trans=x.trans, breaks=c(0, 10^c(log.break:0))) +
  theme(axis.text.x = element_text(angle=45, hjust=1))

ggsave(p, file=str_interp("${scripts.dir}/pdfs/${exp.name}.${ref.species}.pvals_all.pdf"), width=4, height=3)

