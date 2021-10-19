#!/usr/bin/env Rscript
rm(list=ls())

##### Load packages #####
suppressMessages( library(getopt) )
suppressMessages( library(ape) )
suppressMessages( library(data.table) )
suppressMessages( library(stringr) )

##### Get options #####
opt.spec <- matrix(c("help",         "h",  0,  "logical",   0, 
                     "ref.species",  "r",  1,  "character", 1, 
                     "exp.name",     "c",  1,  "character", 1,
                     "out.exp",      "o",  1,  "character", 1,
                     "overwrite",    "w",  0,  "logical",   0
), byrow=TRUE, ncol=5)
opt <- getopt(opt.spec[,1:4])
req.args <- opt.spec[as.logical(as.integer(opt.spec[,5])), 1]

main.dir <- "/moto/palab/projects/male_mutation_bias_XA/mut_sex_bias_amniotes"

###--- Run this commented block to test locally ---#
# main.dir <- "/Users/felix/mt_mp_lab/projects/male_mutation_bias_XA/mut_sex_bias_amniotes"
# opt <- list(ref.species="Gallus_gallus",
#             exp.name="Aves",
#             out.exp="Aves_g1,Aves_g2,Aves_g3,Aves_g4,Aves_g5,Aves_g6")
###--------------------------------------------#

scripts.dir <- str_interp("${main.dir}/scripts")
data.dir    <- str_interp("${main.dir}/data")
pgls.dir    <- str_interp("${main.dir}/pgls_files")
alpha.dir   <- str_interp("${scripts.dir}/alphas")

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

if(is.null(opt$overwrite)){
  opt[["overwrite"]] <- FALSE
}

ref.species <- opt$ref.species
exp.name <- opt$exp.name
out.exp <- opt$out.exp # Output experiments. Sub groupings of a particular order
out.exp <- str_split(out.exp, ",")[[1]]
overwrite <- opt$overwrite

##### Set up file names #####
if(exp.name != "Aves"){
  stop("Currently only needed to run 'Aves' right now")
}

out.fn <- str_interp("${data.dir}/Aves_metadata_and_alphas.csv")
if(file.exists(out.fn)){
  if(overwrite){
    cat(str_interp("This script will overwrite ${out.fn}\n"))
  } else {
    stop(str_interp("Stopping script so as not to overwrite ${out.fn}"))
  }
}

alpha.fn <- str_interp("${scripts.dir}/alphas/${exp.name}.${ref.species}.LM.tsv")
gen.fn <- str_interp("${data.dir}/Aves_traits.tsv")
meta.fn <- str_interp("${data.dir}/Aves_assembly_metadata.csv")

##### LOAD DATA #####
a.dat <- fread(alpha.fn)
g.dat <- fread(gen.fn)
m.dat <- fread(meta.fn)
m.dat <- m.dat[!duplicated(Species)]

##### JOIN DATA #####
a.dat <- a.dat[mut_type == "mod"]

a.dat <- merge(a.dat, g.dat, by.x="species", by.y="Species")
a.dat <- merge(a.dat, m.dat, by.x="species", by.y="Species")

##### ADD COLUMN #####
a.dat[, GenerationLength_d := `Generation time (years)` * 365]

##### Write table to output #####
# Output file
fwrite(a.dat, file=out.fn, sep=",", quote=TRUE)
