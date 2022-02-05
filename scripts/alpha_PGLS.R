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


##### Get options #####
opt.spec <- matrix(c("help",         "h",  0,  "logical",   0, 
                     "exp.name",     "c",  1,  "character", 1,
                     "ee.mut",       "e",  1,  "character", 1  # How many early embryonic mutations to test for
), byrow=TRUE, ncol=5)
opt <- getopt(opt.spec[,1:4])
req.args <- opt.spec[as.logical(as.integer(opt.spec[,5])), 1]

# main.dir <- "/moto/palab/projects/male_mutation_bias_XA/mut_sex_bias_amniotes"

###--- Run this commented block to test locally ---#
main.dir <- "/Users/felix/mt_mp_lab/projects/male_mutation_bias_XA/mut_sex_bias_amniotes"
# opt <- list(exp.name="Birds", ee.mut="5")
# opt <- list(exp.name="Mammals", ee.mut="5")

###--------------------------------------------#

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


exp.name <- opt$exp.name
ee.mut <- as.integer(opt$ee.mut)
# out.exp <- opt$out.exp # Output experiments. Sub groupings of a particular order
# out.exp <- str_split(out.exp, ",")[[1]]
# covariate <- opt$covariate


##### Set up file names #####
alpha.file <- str_interp("${data.dir}/Table_S2.csv")
if(exp.name == "Mammals"){
  # alpha.file  <- str_interp("${main.dir}/data/XA_2exposure-model_alphas.csv")
  tree.file   <- str_interp("${main.dir}/trees/mammals241.TimeTree.nwk")
  phast.file  <- str_interp("${main.dir}/trees/241-mammalian-2020v2.phast-242.nh")
} else if(exp.name == "Birds"){
  # alpha.file  <- str_interp("${main.dir}/data/Birds_metadata_and_alphas.csv")
  tree.file   <- str_interp("${main.dir}/trees/Birds.TimeTree.nwk")
  phast.file  <- str_interp("${main.dir}/trees/363-avian-2020-phast.nh")
} else {
  stop("Currently can only handle 'Mammals' and 'Birds' right now")
}

# nm.map.file <- str_interp("${main.dir}/trees/TimeTree.spname_map.tab")
config.file <- str_interp("${scripts.dir}/pgls_config_file.${exp.name}.txt")

colors.file    <- str_interp("${main.dir}/data/plot_colors.tsv")

out.prefix <- str_interp("${scripts.dir}/pgls_res/Me${ee.mut}/${exp.name}")

# where to cat output to
cat.file    <- str_interp("${out.prefix}.output.txt")
cat("", file=cat.file, append=FALSE)

##### Load data #####
# Read input files
alpha <- fread(alpha.file, sep=",", header=TRUE, na.strings=c("NaN", "", "Na", "NA", "."))
# setnames(alpha, "G", "Generation_time_y")
alpha <- alpha[Group == exp.name]
full.tree <- read.tree(tree.file)
phast.tree <- read.tree(phast.file)



# Replace Alignment species with close timetree equivalent
for(i in 1:nrow(alpha)){
  sp.t <- alpha[i, TimeTree_Species]
  sp.s <- alpha[i, Species]
  
  if(sp.s != sp.t){
    full.tree$tip.label[which(full.tree$tip.label == sp.t)] <- sp.s 
  }
}

# Filter the tree to match species in alpha table
tree.sp <- intersect(full.tree$tip.label, alpha$Species)
tree <- keep.tip(full.tree, tree.sp)

tree[["node.label"]] <- str_c("A", Ntip(tree) + (1:Nnode(tree)))

config <- fread(cmd = str_interp("grep -v '^#' ${config.file}"))
pca.traits     <- config[trait_type == "pca",     x]
pred.traits    <- config[trait_type == "pred",    x]
gnmstat.traits <- config[trait_type == "gnmstat", x]
extra.traits   <- config[trait_type == "other",   x]

# Convert MutPerYearUCSC to units of per Gb per y (default per Mb)
alpha[, MutPerYearUCSC := MutPerYearUCSC * 1e3]


##### Functions #####
## PCA functions
RotateVarimaxPCA <- function(traits, ncomp, rotate_method="varimax"){
  pca <- prcomp(traits, center=TRUE, scale=TRUE)
  
  sink(cat.file, append=TRUE); print(summary(pca)); sink()
  
  rawLoadings     <- pca$rotation[,1:ncomp] %*% diag(pca$sdev, ncomp, ncomp)
  if(rotate_method == "varimax"){
    rotatedLoadings <- varimax(rawLoadings, normalize=TRUE)$loadings
  } else {
    require(GPArotation)
    rotatedLoadings <- quartimax(rawLoadings, normalize=TRUE)$loadings
  }
  invLoadings     <- t(pracma::pinv(rotatedLoadings))
  
  pcs_rotated    <- as.data.frame(scale(traits) %*% invLoadings)
  colnames(pcs_rotated) <- str_c("RC", 1:ncomp)
  
  pcs <- as.data.frame(pca$x[,1:ncomp])
  pcs_rotated$Species <- rownames(pcs)
  pc_og_rotated <- cbind(pcs, pcs_rotated)
  
  loadings <- as.data.frame(rawLoadings)
  colnames(loadings) <- str_c("PC", 1:ncomp)
  loadings_rotated <- as.data.frame(rotatedLoadings[,1:ncomp])
  colnames(loadings_rotated) <- paste("RC", seq.int(1, ncomp), sep="")
  return(list("pc_rc" = pc_og_rotated, "loadings" = loadings, "loadings_rotated" = loadings_rotated, "pca"= pca))
}

# PCA Plotting with ggplot
BiplotPCA <- function(pca, loadings, latin2comm=TRUE){
  scalar <- 1
  loadings$traits <- row.names(loadings)
  colnames(loadings) <- c("C1","C2","traits")
  colnames(pca) <- c("Species", "C1", "C2")
  if(latin2comm){
    for(i in seq_along(pca$Species)){
      pca$Species[i] <- alpha[Species == pca$Species[i], Common_name]
    }
  }
  
  gg.obj <- ggplot(data=pca, aes(x=C1, y=C2)) + 
    geom_point() + geom_text_repel(aes(label=Species)) +
    geom_segment(data=loadings, 
                 aes(x=0, y=0, xend=C1*scalar, yend=C2*scalar), 
                 arrow = arrow(length = unit(0.03, "npc")),
                 color="cornflowerblue") + 
    geom_text_repel(data=loadings, aes(x=C1*scalar,y=C2*scalar,label=traits),
                    color="cornflowerblue") +
    geom_hline(yintercept = 0, linetype="dashed", alpha=0.5) + 
    geom_vline(xintercept = 0, linetype="dashed", alpha=0.5) +
    theme(axis.line=element_blank(), panel.border = element_blank(),
          panel.background = element_blank())
  
  return(gg.obj)
}



##### Add extra data columns #####


# Compute divergence from nearest species in order, nearest chromosome level assembly
MinDiv2Tips <- function(phy, target.tip, comparison.tips){
  div.vals <- c()
  if(target.tip %in% comparison.tips){
    out.list <- list("tip"=target.tip, "div"=0)
  } else {
    for(i in comparison.tips){
      div.vals[i] <- sum(keep.tip(phy, c(target.tip, i))$edge.length)
    }
    
    j <- which.min(div.vals)
    out.list <- list("tip"=comparison.tips[j], "div"=unname(div.vals[j]))
  }
  
  return(out.list)
}

chrom.level.sp <- alpha[AssemblyStatus == "Chromosome", Species]
alpha[, c("Nearest_chrom_assembly", "Divergence_to_NCA") := MinDiv2Tips(phast.tree, Species, chrom.level.sp), by=Species]


##### Perform PCA, plot #####
if(length(pca.traits) > 0){
  # Choose PCA functions
  pca_dat <- as.data.frame(na.omit(alpha[, c("Species", pca.traits), with=FALSE])) # drop NAs
  cat(str_interp("\nKeeping ${nrow(pca_dat)} of ${nrow(alpha)} species with complete trait info for PCA.\n"),
      file=cat.file, append=TRUE)
  row.names(pca_dat) <- pca_dat$Species
  pca_dat <- log10(pca_dat[, pca.traits])
  
  rotate_method <- "varimax"
  # NB: Played around with quartimax rotation but it doesn't seem to differ all that
  # much from the raw PCs
  
  pca_and_varimax <- RotateVarimaxPCA(pca_dat, 2, rotate_method)
  p <- fviz_pca_biplot(pca_and_varimax$pca, repel=TRUE)
  print(p)
  
  # Plot
  pca_prop_var <- summary(pca_and_varimax$pca)$importance[2,]
  p <- BiplotPCA(pca_and_varimax$pc_rc[,c("Species","PC1","PC2")], pca_and_varimax$loadings) +
    xlab(str_interp("Component 1 ($[0.1f]{pca_prop_var[1]*100}%)")) +
    ylab(str_interp("Component 2 ($[0.1f]{pca_prop_var[2]*100}%)"))
  print(p)
  ggsave(p, filename=str_interp("${out.prefix}.PCA.pdf"), 
         width=6, height=6)
  
  p <- BiplotPCA(pca_and_varimax$pc_rc[,c("Species","RC1","RC2")], pca_and_varimax$loadings_rotated) +
    xlab(str_interp("${str_to_title(rotate_method)} Component 1")) +
    ylab(str_interp("${str_to_title(rotate_method)} Component 2"))
  print(p)
  ggsave(p, filename=str_interp("${out.prefix}.PCA.${rotate_method}_rotated.pdf"),
         width=6, height=6)
  
  # Merge PCA with alpha table
  alpha <- merge(alpha, pca_and_varimax$pc_rc, by.x='Species', by.y="Species", all.x=TRUE)
  
  # Add PCs and RCs to configuration table
  tmp <- data.table(x=c("PC1", "PC2", "RC1", "RC2"), xlog=FALSE, 
                    y="Alpha", ylog=FALSE,
                    rm_species=NA, trait_type="pca")
  config <- rbindlist(list(config, tmp))
  pca.traits <- c(pca.traits, c("PC1", "PC2", "RC1", "RC2"))
}

##### Run PGLS & PIC #####
xy_data <- list()
pic_data <- list()
models <- list()

cor.method <- "spearman" # kendall or spearman

for(cur.trait in pred.traits){
  if(!(cur.trait %in% names(alpha))){
    alpha[[cur.trait]] <- alpha[[str_interp("Predicted_alpha_evo_Me${ee.mut}")]]
  }
}


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
models[, ols_padj     := p.adjust(ols_pval, "BH")]
models[, pic_padj     := p.adjust(pic_pval, "BH")]

# Create plot strings
lambda_char <- "\u03BB"
sq_char <- "\u00B2"
cor_char <- c("spearman"="\u03C1", "kendall"="\u03C4")
include_padj <- FALSE # Include adjusted pvalues in plot?

MakeModelStr <- function(lambda.val, p.val, rsq.val, model.type = "ML", prefix="\n"){
  if(model.type == "ML"){
    out.str <- str_interp("${prefix} ML(${lambda_char})=$[0.2f]{lambda.val}")
  } else if(model.type == "full"){
    out.str <- str_interp("${prefix} ${lambda_char}=$[d]{round(lambda.val)}")
  }
  out.str <- str_interp("${out.str},  p = $[0.3f]{p.val}, r${sq_char} = $[0.2f]{rsq.val}")
  return(out.str)
}


models[, pgls_ml_str := MakeModelStr(ml_lambda, pgls_ml_pval, pgls_ml_rsq,   "ML", prefix="\n"),     by=experiment]
models[, pgls_str    := MakeModelStr(        1,    pgls_pval,    pgls_rsq, "full", prefix="\n\n"), by=experiment]


# PIC stats
models[, cor_rho_str        := sprintf("\n %s %s = %0.2f", str_to_title(cor.method), cor_char[cor.method], cor_rho)]
models[, pic_cor_pval_str := sprintf("\n\n p = %0.4f", pic_cor_pval)]



##### Perform node height tests #####
x <- alpha[, Alpha]
names(x) <- alpha$Species
x.pic <- pic(x, tree)
nheights <- dist.nodes(tree)[find_root(tree), as.integer(str_remove(names(x.pic), "A"))]
nh.res <- lm(abs(x.pic) ~ nheights)
nh.res.coef <- coef(summary(nh.res))
nh.pval <- nh.res.coef[2,4]
nh.r2 <- summary(nh.res)$r.squared

plt.dat <- data.frame(x=nheights, y=abs(x.pic), 
                      nodenames=names(x.pic))
lab.dat <- data.frame(nh_r2_str = str_interp("\n r${sq_char} = $[0.3f]{nh.r2}"),
                      nh_pval_str = str_interp("\n\n p = $[0.3f]{nh.pval}"))

p <- ggplot(aes(x=x, y=y), data=plt.dat) +
  geom_point(alpha=0.7, col="#333333") +
  geom_smooth(method="lm") +
  geom_text(aes(label=nh_pval_str, x = -Inf, y = Inf), data=lab.dat, color="black", 
            hjust = 0, vjust = 1, size = 3) +
  geom_text(aes(label=nh_r2_str, x = -Inf, y = Inf), data=lab.dat, color="black", 
            hjust = 0, vjust = 1, size = 3) +
  xlab("Branching times") + ylab("Absolute value of contrasts")
print(p)

ggsave(p, filename=str_interp("${out.prefix}.nheight.pdf"), 
       width=4.5, height=4)

# geiger::tips(tree, 148)
# geiger::tips(tree, 150)


##### TEST #####
# for(i in names(pic_data)[1]){
#   pic.dat <- pic_data[[i]]
#   raw.dat <- xy_data[[i]]
#   
#   pic.lm <- lm(PIC2 ~ PIC1 + 0, data=pic.dat)
#   raw.lm <- lm(yvar ~ xvar, data=raw.dat)
#   
#   p <- ggplot(aes(x=PIC1, y=PIC2), data=pic.dat) + geom_point() + 
#     geom_smooth(method="lm", formula= y ~ x + 0)
#   print(p)
#   
#   p <- ggplot(aes(x=xvar, y=yvar), data=raw.dat) + geom_point() + 
#     geom_smooth(method="lm")
#   print(p)
#   
#   print(summary(pic.lm))
#   print(summary(raw.lm))
# }
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








