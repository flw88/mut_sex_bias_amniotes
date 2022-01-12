#!/usr/bin/env Rscript
rm(list=ls())

##### Load packages #####
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

opt.spec <- matrix(c("help",         "h",  0,  "logical",   0, 
                     "ee.mut",       "e",  1,  "character", 1  # How many early embryonic mutations to test for
), byrow=TRUE, ncol=5)
opt <- getopt(opt.spec[,1:4])
req.args <- opt.spec[as.logical(as.integer(opt.spec[,5])), 1]

main.dir <- "/Users/felix/mt_mp_lab/projects/male_mutation_bias_XA/mut_sex_bias_amniotes"
# main.dir <- "/moto/palab/projects/male_mutation_bias_XA/mut_sex_bias_amniotes"

scripts.dir <- str_interp("${main.dir}/scripts")
data.dir    <- str_interp("${main.dir}/data")
pgls.dir    <- str_interp("${main.dir}/pgls_files")

source(str_interp("${scripts.dir}/PGLS.funcs.R"))

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


ee.mut <- as.integer(opt$ee.mut)

exp.name <- "DNMs"


##### Set up file names #####
alpha.file <- str_interp("${data.dir}/Table_S1_ee${ee.mut}.csv")

tree.file   <- str_interp("${main.dir}/trees/DNMs.TimeTree.nwk")

config.file <- str_interp("${scripts.dir}/pgls_config_file.${exp.name}.txt")

colors.file    <- str_interp("${main.dir}/data/plot_colors.tsv")

out.prefix <- str_interp("${scripts.dir}/pgls_res/ee${ee.mut}/${exp.name}")

# where to cat output to
cat.file    <- str_interp("${out.prefix}.output.txt")
cat("", file=cat.file, append=FALSE)

##### Load data #####
# Read input files
alpha <- fread(alpha.file, sep=",", header=TRUE, na.strings=c("NaN", "", "Na", "NA", "."))
alpha <- alpha[!is.na(Alpha_dnm) & !(Group %in% c("Birds", "Snakes"))] # Filter birds and snakes and non-dnm

full.tree <- read.tree(tree.file)


# Filter the tree to match species in alpha table
tree.sp <- intersect(full.tree$tip.label, alpha$Species)
tree <- keep.tip(full.tree, tree.sp)

tree[["node.label"]] <- str_c("A", Ntip(tree) + (1:Nnode(tree)))

config <- fread(cmd = str_interp("grep -v '^#' ${config.file}"))

##### Add data columns #####
for(x in config$x){
  if(!(x %in% names(alpha))){
    alpha[[x]] <- alpha$Predicted_alpha_dnms
  }
}

##### Run PGLS & PIC #####
xy_data <- list()
pic_data <- list()
models <- list()

cor.method <- "spearman" # kendall or spearman

# Convert AssemblyStatus to integers
alpha[, AssemblyStatus := as.factor(AssemblyStatus)]
config[, skipped := FALSE]

# Configure each trait variable (to be compared with alpha)
hdr <- c("Comparison", "pgls_ml_pvalue")
cat(str_interp("\n$[30s]{hdr[1]}   ${hdr[2]}\n"), file=cat.file, append=TRUE)

for (row in config[, .I[!skipped]]){
  # X and Y
  x <- config[row, x]
  y <- config[row, y]
  xlog <- config[row, xlog]
  ylog <- config[row, ylog]
  experiment <- str_c(x, ".", y)
  
  rm.species <- str_split(config[row, rm_species], ",")[[1]]
  
  cur.alpha <- alpha[!(Species %in% rm.species)][]
  
  # Run pgls and pic
  pgls_ml <- runPGLS(cur.alpha, x, y, xlog, ylog, tree, lambda="ML", kappa=1, delta=1)
  if (pgls_ml[1]=="NA"){
    cat(str_interp("Skipping ${experiment}, no observations left"), file=cat.file, append=TRUE)
    config[row, skipped := TRUE]
    next
  }
  pgls <- runPGLS(cur.alpha, x, y, xlog, ylog, tree, lambda=1, kappa=1, delta=1)
  pic  <- picModel(cur.alpha, x, y, xlog, ylog, tree)
  if(!(class(cur.alpha[[x]]) %in% c("factor", "character"))){
    pic.cor <- picCor(cur.alpha, x, y, xlog, ylog, tree, method=cor.method, pic.pairwise=FALSE)
  }
  
  pgls$data$experiment <- experiment
  
  # Keep x/y data
  xy_data[[experiment]] <- as.data.table(pgls$data) #rbind(xy_data, pgls$data)
  
  # Keep pics from correlation test
  pic_data[[experiment]] <- data.table(PIC1 = pic.cor$PIC1,
                                       PIC2 = pic.cor$PIC2,
                                       PIC_names = names(pic.cor$PIC1),
                                       experiment = experiment)
  
  # And keep model info
  pgls_ml_pval <- extract_p(pgls_ml$model, 2)
  pgls_pval    <- extract_p(pgls$model   , 2)
  
  pic_pval <- extract_p(pic$model, 1)
  pic_cor_pval <- pic.cor$cor$p.value
  
  pgls_ml_rsq <- summary(pgls_ml$model)$adj.r.squared[1]
  pgls_rsq    <- summary(pgls$model   )$adj.r.squared[1]
  
  ml_lambda <- unname(summary(pgls_ml$model)$param["lambda"])
  cor_rho <- unname(pic.cor$cor$estimate)
  
  cat(str_interp("$[30s]{experiment}:  $[0.4f]{pgls_ml_pval}\n"), file=cat.file, append=TRUE)
  
  pgls_ml_intercept <- coef(pgls_ml$model)[[1]]
  pgls_ml_slope <- coef(pgls_ml$model)[[2]]
  pgls_intercept <- coef(pgls$model)[[1]]
  pgls_slope <- coef(pgls$model)[[2]]
  
  models[[experiment]] <- data.table(pgls_ml_intercept=pgls_ml_intercept, 
                                     pgls_ml_slope=pgls_ml_slope, 
                                     pgls_intercept=pgls_intercept, pgls_slope=pgls_slope,
                                     pgls_ml_pval=pgls_ml_pval, pgls_pval=pgls_pval, 
                                     pic_pval=pic_pval, pic_cor_pval=pic_cor_pval,
                                     pgls_ml_rsq=pgls_ml_rsq, pgls_rsq=pgls_rsq,
                                     ml_lambda=ml_lambda, cor_rho=cor_rho,
                                     experiment=experiment)
}


##### Process model summary stats #####

models <- rbindlist(models)

# Compute adjusted p-value by Benjamini-Hochberg
models[, pgls_ml_padj := p.adjust(pgls_ml_pval, "BH")]
models[, pgls_padj    := p.adjust(pgls_pval, "BH")]
models[, pic_padj     := p.adjust(pic_pval, "BH")]


##### SAVE MODEL RESULTS #####
save.cols <- str_subset(names(models), pattern="_str$", negate=TRUE)
fwrite(models[, save.cols, with=FALSE], 
       file=str_interp("${out.prefix}.model_params.tsv"), 
       col.names=TRUE, sep="\t")

##### SAVE XY DATA #####
if("AssemblyStatus.Alpha" %in% names(xy_data)){
  xy_data$`AssemblyStatus.Alpha`[, xvar := as.integer(xvar)]
}

out.dat <- rbindlist(xy_data)
fwrite(out.dat, 
       file=str_interp("${out.prefix}.xy_data.tsv"), 
       col.names=TRUE, sep="\t")
# warnings()
