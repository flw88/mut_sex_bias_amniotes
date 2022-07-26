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
# opt <- list(ee.mut="5")

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


ee.mut <- opt$ee.mut
exp.name <- "DNMs"


##### Set up file names #####
alpha.file <- str_interp("${data.dir}/Table_S2.csv")

tree.file   <- str_interp("${main.dir}/trees/DNMs.TimeTree.nwk")

config.file <- str_interp("${scripts.dir}/pgls_config_file.${exp.name}.txt")

colors.file    <- str_interp("${main.dir}/data/plot_colors.tsv")

out.prefix <- str_interp("${scripts.dir}/pgls_res/Me${ee.mut}/${exp.name}")

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
    alpha[[x]] <- alpha[[str_interp("Predicted_alpha_dnm_Me${ee.mut}")]]
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
  
  # if (pgls_ml[1]=="NA"){
  #   cat(str_interp("Skipping ${experiment}, no observations left"), file=cat.file, append=TRUE)
  #   config[row, skipped := TRUE]
  #   next
  # }
  pgls <- runPGLS(cur.alpha, x, y, xlog, ylog, tree, lambda=1, kappa=1, delta=1)
  pic  <- picModel(cur.alpha, x, y, xlog, ylog, tree)
  if(!(class(cur.alpha[[x]]) %in% c("factor", "character"))){
    pic.cor <- picCor(cur.alpha, x, y, xlog, ylog, tree, method=cor.method, pic.pairwise=FALSE)
  }
  
  # Ordinary least squares
  ols <- lm(yvar ~ xvar, data=pgls$data)
  
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
  ols_pval     <- extract_p(ols          , 2)
  
  pic_pval <- extract_p(pic$model, 1)
  pic_cor_pval <- pic.cor$cor$p.value
  
  pgls_ml_rsq <- summary(pgls_ml$model)$r.squared[1]
  pgls_rsq    <- summary(pgls$model   )$r.squared[1]
  ols_rsq     <- summary(ols          )$r.squared[1]
  
  ml_lambda <- unname(summary(pgls_ml$model)$param["lambda"])
  
  # LRT for lambda=0 versus ML lambda estimate
  ml_lambda_pval <- pgls.profile(pgls_ml$model)$ci$bounds.p[1] 
  
  cor_rho <- unname(pic.cor$cor$estimate)
  
  cat(str_interp("$[30s]{experiment}:  $[0.4f]{pgls_ml_pval}\n"), file=cat.file, append=TRUE)
  
  pgls_ml_intercept <- coef(pgls_ml$model)[[1]]
  pgls_ml_slope <- coef(pgls_ml$model)[[2]]
  
  pgls_intercept <- coef(pgls$model)[[1]]
  pgls_slope <- coef(pgls$model)[[2]]
  
  
  ols_intercept <- coef(ols)[[1]]
  ols_slope <- coef(ols)[[2]]
  
  models[[experiment]] <- data.table(pgls_ml_intercept=pgls_ml_intercept, pgls_ml_slope=pgls_ml_slope, 
                                     pgls_intercept=pgls_intercept,       pgls_slope=pgls_slope,
                                     ols_intercept=ols_intercept,         ols_slope=ols_slope,
                                     pgls_ml_pval=pgls_ml_pval, pgls_pval=pgls_pval, ols_pval=ols_pval,
                                     pic_pval=pic_pval, pic_cor_pval=pic_cor_pval,
                                     pgls_ml_rsq=pgls_ml_rsq, pgls_rsq=pgls_rsq, ols_rsq,
                                     ml_lambda=ml_lambda, ml_lambda_pval=ml_lambda_pval, cor_rho=cor_rho,
                                     experiment=experiment)
}


##### Process model summary stats #####
models <- rbindlist(models)

# Compute adjusted p-value by Benjamini-Hochberg
models[, pgls_ml_padj := p.adjust(pgls_ml_pval, "BH")]
models[, pgls_padj    := p.adjust(pgls_pval, "BH")]
models[, pic_padj     := p.adjust(pic_pval, "BH")]


print(models[, .(experiment, pgls_better=(ols_pval > pgls_pval), 
                 ols_pval=round(ols_pval, digits=3), 
                 pgls_pval=round(pgls_pval, digits=3), 
                 ols_rsq=round(ols_rsq, digits=2), 
                 pgls_rsq=round(pgls_rsq, digits=2),
                 ml_lambda_pval=round(ml_lambda_pval, digits=3) )])

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

##### TEST #####
# require(ggtree)
# # i <- "Predicted_alpha_dnms_min2trios.Alpha_dnm"
# i <- "Predicted_alpha_dnms_min8trios.Alpha_dnm"
# pic.dat <- pic_data[[i]]
# raw.dat <- xy_data[[i]]
# 
# pic.lm <- lm(PIC2 ~ PIC1 + 0, data=pic.dat)
# raw.lm <- lm(yvar ~ xvar, data=raw.dat)
# 
# p <- ggplot(aes(x=PIC1, y=PIC2), data=pic.dat) + geom_point() + geom_text_repel(aes(label=PIC_names)) + 
#   geom_smooth(method="lm", formula= y ~ x + 0)
# print(p)
# 
# p <- ggplot(aes(x=xvar, y=yvar), data=raw.dat) + geom_point() +
#   geom_smooth(method="lm")
# print(p)
# 
# print(summary(pic.lm))
# print(summary(raw.lm))
# 
# # Tree
# rm.species <- str_split(config[str_c(x, y, sep=".") == i, rm_species], ",")[[1]]
# sub.tree <- keep.tip(tree, setdiff(tree$tip.label, rm.species))
# tmp.dat <- alpha[Species %in% sub.tree$tip.label, .(Species, Predicted_alpha_dnms, Alpha_dnm)]
# setkey(tmp.dat, Species)
# sub.tree$tip.label <- tmp.dat[sub.tree$tip.label, str_c(Species, "\n(", round(Predicted_alpha_dnms, 2), ", ", round(Alpha_dnm, 2), ")")]
# 
# p <- ggplot(sub.tree) + geom_tree() + theme_tree() + geom_tiplab() + geom_nodelab() + xlim(0, 125)
# print(p)
# 
# 
# print(tmp.dat[c("Homo_sapiens", "Pan_troglodytes", "Gorilla_gorilla")])
