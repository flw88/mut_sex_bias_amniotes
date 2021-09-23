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

### Get options      varName              flag  takesArg  type       isRequired         
opt.spec <- matrix(c("help",                "h",  0,     "logical",   0, 
                     "ref.species",         "r",  1,     "character", 1, 
                     "experiment.name",     "c",  1,     "character", 1, 
                     "poisson.glm",         "p",  0,     "logical",   0, # Use Poisson regression (DEFAULT: linear regression)
                     "gc.squared",          "g",  0,     "logical",   0, # Include a gc^2 term in regression (DEFAULT: no gc^2 term)
                     "weighted",            "w",  0,     "logical",   0, # Perform weighted regression (DEFAULT: unweighted)
                     "unrestrict.features", "u",  0,     "logical",   0, # Don't restrict features to X range (DEFAULT: restrict)
                     "out.prefix",          "o",  1,     "character", 0, # Output prefix (DEFAULT: experiment.name)
                     "features",            "f",  1,     "character", 1, # Comma-separated list of features
                     "in.file",             "i",  1,     "character", 0, # Specify the merged_features file explicitly. (DEFAULT: "${scripts.dir}/merged_features/${experiment.name}.features.txt.gz")
                     "seed",                "s",  1,     "character", 1, # Seed for randomization
                     "is.zw",               "z",  0,     "logical",   0  # ZW system (not XY)
                     ), byrow=TRUE, ncol=5)
opt <- getopt(opt.spec[,1:4])
req.args <- opt.spec[as.logical(as.integer(opt.spec[,5])), 1]

scripts.dir <- "/moto/palab/projects/male_mutation_bias_XA/mut_sex_bias_amniotes/scripts"

###--- Run this commented block to test locally ---#
#scripts.dir <- "/Users/felix/mt_mp_lab/projects/male_mutation_bias_XA/scripts"
#opt <- list(ref.species="Bos_taurus",
#          experiment.name="Artiodactyla",
#          # in.file="/Users/felix/mt_mp_lab/male_mutation_bias_XA/scripts/effect_of_regions/merged_features/Euarchonta.features.txt.gz",
#          poisson.glm=TRUE,
#          gc.squared=TRUE,
#          unrestrict.features=TRUE,
#          weighted=FALSE,
#          seed=14712,
#          features="gc",
#          out.prefix="test.Artiodactyla")
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
experiment.name <- opt$experiment.name
feat.vars <- str_split(opt$features, ",")[[1]]
set.seed(as.integer(opt$seed))

poisson.glm <- !is.null(opt$poisson.glm) && opt$poisson.glm
gc.squared <- !is.null(opt$gc.squared) && opt$gc.squared
weighted <- !is.null(opt$weighted) && opt$weighted
unrestrict.features <- !is.null(opt$unrestrict.features) && opt$unrestrict.features
is.zw <- !is.null(opt$is.zw) && opt$is.zw

if(weighted && poisson.glm){
  stop("Poisson GLM has not been implemented for weighted regression")
}

if( is.null(opt$out.prefix) ){
  out.prefix <- experiment.name
} else{
  out.prefix <- opt$out.prefix
}

if( is.null(opt$in.file) ){
  in.file <- str_interp("${scripts.dir}/merged_features/${experiment.name}.${ref.species}.features.txt.gz")
} else {
  in.file <- opt$in.file
}



# Get chrX(Z) and chrY(W) names
sp.to.chrom <- fread(str_interp("../data/Species_to_chromosomes.txt"), header=FALSE,
                     col.names=c("REF", "XZ_NAME", "YW_NAME"), na.strings="NaN")
setkey(sp.to.chrom, REF)

chrXZ <- sp.to.chrom[(ref.species), XZ_NAME]
chrYW <- sp.to.chrom[(ref.species), YW_NAME]

# Min sequence size
min.seq.size <- 1e4


mut.types <- c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G") # Mutation types
sw.types  <- c("S>W", "W>S", "S>S", "W>W") # Strong/weak types
bgc.types <- c("gc-BGC-affected", "gc-BGC-unaffected") # BGC




##### SET UP ENVIORNMENT #####
# Source functions
source(str_interp("${scripts.dir}/alpha_regression.funcs.R"))

# ggplot theme
theme_set(theme_bw() +
            theme(axis.text=element_text(size=12), panel.border=element_rect(size = 1.5)))


##### READ IN DATA #####
dat <- fread(in.file, header=TRUE, fill=TRUE, na.strings=c("", "NaN", "nan"))
dat[, INDEX := 1:.N]

dat[,               region := "A"]
dat[chrom == chrXZ, region := "XZ"]
dat[chrom == chrYW, region := "YW"]

sp.names <- read.tree(text=dat[!is.na(tree), tree][1])$tip.label
cat(str_interp("Starting with ${dat[,.N]} regions\n"))

##### FILTER MISSING DATA #####
dat <- FilterAndReport(dat, dat[, is.na(tree)], "Filtering NaN trees:\n")
# for(i in feat.vars){
#   dat <- FilterAndReport(dat, dat[, is.na(get(i))], str_interp("Filtering NaN ${i}:\n"))
# }
dat <- FilterAndReport(dat, dat[, size < min.seq.size], 
                          str_interp("Filtering sequences < ${min.seq.size} bp:\n"))
dat <- FilterAndReport(dat, dat[, region == "YW"], "Filtering Y regions:\n")

##### PARSE TREE STRINGS #####
# Tip sub rates (also convert to counts)
cat("Obtaining sub rate values at tips\n")
mod.rate.names <- str_c(sp.names, "_mod")
dat[, (mod.rate.names) := ParseTree(tree), by=INDEX]
mod.num.names <- str_c(mod.rate.names, "_num")
mod.den.names <- str_c(mod.rate.names, "_den")
dat[, (mod.den.names) := as.list(rep(size, length(mod.den.names))), by=INDEX]
dat[, (mod.num.names) := as.list(size * .SD[, mod.rate.names, with=FALSE]), by=INDEX]


# Variance in tip sub rates
var.names <- str_c(sp.names, "_var", sep="")
if(weighted){
  cat("Obtaining sub rate variance\n")
  base.tree <- read.tree(text=dat[!is.na(tree), tree][1])
  dat[, (var.names) := ParseTreeVariance(tree_var, base.tree), by=INDEX]
} else {
  cat("Running unweighted.\n")
}

# Substitution counts
for(i in c(mut.types, "C", "T")){
  count.names <- str_c(sp.names, i, sep="_")
  if(str_length(i) == 1){
    msg <- str_interp("Obtaining ${i} base counts\n")
  } else {
    msg <- str_interp("Obtaining ${i} sub counts\n")
  }
  cat(msg)
  base.tree <- read.tree(text=dat[!is.na(tree), tree][1])
  cur.col <- str_c("tree_", i)
  dat[, (count.names) := ParseTree(get(cur.col)), by=INDEX]
}

# 1) Rename columns such that e.g. Homo_sapiens_C>T denotes rate,
#    Homo_sapiens_C>T_num denotes numerator, Homo_sapiens_C>T_den denotes denominator
# 2) Collate mutation type totals
# 3) Collate Strong/Weak types
# 4) Collate gBGC/non-gBGC types
for(sp in sp.names){
  # (1)
  # Rename numerator
  rate.names <- str_c(sp, mut.types, sep="_")
  num.names  <- str_c(rate.names, "num", sep="_")
  den.names  <- str_c(rate.names, "den", sep="_")
  setnames(dat, rate.names, num.names)
  
  # Calculate denominator, rate
  base.names <- str_c(sp, str_sub(mut.types, 1, 1), sep="_")
  for(i in seq_along(rate.names)){
    dat[, eval(den.names[i])  := get(base.names[i])]
    dat[, eval(rate.names[i]) := get(num.names[i]) / get(den.names[i])]
  }
  
  # (2)
  cat("Collating total mutation counts\n")
  ets.rate.names <- str_c(sp, "exptotsub", sep="_")
  ets.num.list   <- list(num.names)
  ets.den.list   <- list(str_c(sp, c("C", "T"), sep="_"))
  dat <- CombineMutTypes(dat, ets.rate.names, ets.num.list, ets.den.list)
  
  
  # (3)
  cat("Collating strong/weak mutation counts\n")
  sw.rate.names <- str_c(sp, sw.types, sep="_")
  sw.num.list   <- sapply(sw.types, function(x) str_c(sp, SWtoType(x), "num", sep="_"))
  sw.den.list   <- as.list(str_c(sp, str_replace(str_sub(sw.types, 1, 1), c("S", "W"), c("C", "T")), sep="_"))
  dat <- CombineMutTypes(dat, sw.rate.names, sw.num.list, sw.den.list)
  
  # (4)
  cat("Collating BGC/non-BGC mutation counts\n")
  bgc.rate.names <- str_c(sp, bgc.types, sep="_")
  bgc.num.list   <- sapply(bgc.types, function(x) str_c(sp, BGCtoType(x), "num", sep="_"))
  bgc.den.list   <- as.list(rep(str_c(sp, "exptotsub_den", sep="_"), 2))
  dat <- CombineMutTypes(dat, bgc.rate.names, bgc.num.list, bgc.den.list)
  
}


### Filter
# Filter sub rates >= 1
rm.mask <- apply(dat[, mod.rate.names, with=FALSE], 1, function(x) any(x >= 1))
dat <- FilterAndReport(dat, rm.mask, "Filtering sub rates >= 1:\n")

# Filter negative variance
if(weighted){
  rm.mask <- apply(dat[, var.names, with=FALSE], 1, function(x) any(x <= 0))
  dat <- FilterAndReport(dat, rm.mask, "Filtering sub rate variance <= 0:\n")
}

##### REGRESSION #####
# save.image(file=str_interp("${scripts.dir}/lm_res/${out.prefix}.RData"))
# stop("debug")
### Perform regression
sum.types <- c("mod", "exptotsub")

lm.list <- list(); samp.coef.list <- list()
for(mt in c(sum.types, mut.types)){ lm.list[[mt]] <- list(); samp.coef.list[[mt]] <- list() }

for(sp in sp.names){
  cat(str_interp("Performing regression on ${sp}\n"))
  
  if(weighted){ 
    variance.var <- str_interp("${sp}_var") 
  } else {
    variance.var <- NULL
  }
  
  # Use species-specific gc
  if("gc" %in% feat.vars){
    gc.col <- str_interp("${sp}-gc")
    if(!(gc.col %in% names(dat))){ gc.col <- "gc" }
    use.vars <- c(gc.col, feat.vars[feat.vars != "gc"])
    names(use.vars) <- c("gc", use.vars[-1])
  } else {
    use.vars <- feat.vars
    names(use.vars) <- feat.vars
  }
  
  # By type
  for(mt in c(sum.types, mut.types, sw.types, bgc.types)){
    cat(str_interp("Regressing ${sp}, mutation type ${mt}\n"))
    response.var <- str_c(sp, "_", mt)
    lm.list[[mt]][[sp]] <- RegressRates(dat, response.var, use.vars, variance.var=variance.var, 
                                        intercept.only=TRUE, unrestrict.features=unrestrict.features,
                                        poisson.glm=poisson.glm, gc.squared=gc.squared)
    
    samp.coef.list[[mt]][[sp]] <- ResampleRegression(lm.list[[mt]][[sp]], by.region=TRUE)
  }
}

##### COLLECT ALPHA RESULTS #####
alpha.dat <- data.table()
for(mt in c(sum.types, mut.types, sw.types, bgc.types)){
  for(sp in sp.names){
    new.row <- as.list(AlphaFromRegression(lm.list[[mt]][[sp]], feat.vars, samp.coef.list[[mt]][[sp]], is.zw=is.zw))
    new.row[["species"]]  <- sp
    new.row[["mut_type"]] <- mt
    alpha.dat <- rbindlist(list(alpha.dat, as.data.table(new.row)))
  }
}

fwrite(alpha.dat, file=str_interp("${scripts.dir}/alphas/${out.prefix}.LM.tsv"),
       sep="\t", quote=FALSE)


### Print summaries of regression results
CatRegressionSummary(lm.list[["mod"]], out.file=str_interp("${scripts.dir}/lm_res/${out.prefix}.lm_summaries.txt"))

### Save workspace for debugging, etc.
save.image(file=str_interp("${scripts.dir}/lm_res/${out.prefix}.RData"))


##### PLOT #####
plt.nrow <- floor(sqrt(length(sp.names)))
plt.ncol <- ceiling(length(sp.names) / plt.nrow)
for(i in feat.vars){
  # Sub rate vs. features
  gg.obj <- PlotRegression(lm.list[["mod"]], i, samp.coef.list[["mod"]], facet.nrow=plt.nrow)

  # !! Adjust height and width as needed depending on number of species included
  ggsave(gg.obj, filename=str_interp("${scripts.dir}/pdfs/${out.prefix}.rate_vs_${i}.pdf"),
         height=plt.nrow*3, width=plt.ncol*4)
}


