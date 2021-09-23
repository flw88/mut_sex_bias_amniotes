#!/usr/bin/env Rscript
rm(list=ls())

##### Load packages #####
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


##### Get options #####
opt.spec <- matrix(c("help",         "h",  0,  "logical",   0, 
                     "ref.species",  "r",  1,  "character", 1, 
                     "exp.name",     "c",  1,  "character", 1,
                     "out.exp",      "o",  1,  "character", 1
), byrow=TRUE, ncol=5)
opt <- getopt(opt.spec[,1:4])
req.args <- opt.spec[as.logical(as.integer(opt.spec[,5])), 1]

main.dir <- "/moto/palab/projects/male_mutation_bias_XA/mut_sex_bias_amniotes"

###--- Run this commented block to test locally ---#
main.dir <- "/Users/felix/mt_mp_lab/projects/male_mutation_bias_XA/mut_sex_bias_amniotes"
opt <- list(ref.species="Homo_sapiens",
            exp.name="Mammals",
            out.exp="Mammals")
###--------------------------------------------#

scripts.dir <- str_interp("${main.dir}/scripts")
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

ref.species <- opt$ref.species
exp.name <- opt$exp.name
out.exp <- opt$out.exp # Output experiments. Sub groupings of a particular order
out.exp <- str_split(out.exp, ",")[[1]]

##### Set up file names #####
if(exp.name != "Mammals"){
  stop("Currently can only handle 'Mammals' right now")
}

alpha.file  <- str_interp("${main.dir}/data/XA_2exposure-model_alphas.csv")
tree.file   <- str_interp("${main.dir}/trees/mammals241.TimeTree.nwk")
ucsc.file   <- str_interp("${main.dir}/trees/241-mammalian-2020v2.phast-242.nh")
nm.map.file <- str_interp("${main.dir}/trees/mammals241.TimeTree.spname_map.tab")
config.file <- str_interp("${scripts.dir}/pgls_config_file.txt")
latin.file  <- str_interp("${main.dir}/data/latin2common.txt")

anage.file     <- str_interp("${main.dir}/pgls_files/anage_data.txt")
anage.hdr.file <- str_interp("${main.dir}/pgls_files/anage_simplified_header.txt")
impute.file    <- str_interp("${main.dir}/pgls_files/impute_species.txt")

# where to cat output to
cat.file    <- str_interp("${scripts.dir}/pgls_res/${exp.name}.output.txt")
cat("", file=cat.file, append=FALSE)

##### Load data #####
# Read input files
alpha <- fread(alpha.file)
full.tree <- read.tree(tree.file)
ucsc.tree <- read.tree(ucsc.file)

latin.dat <- fread(latin.file)
setkey(latin.dat, Species)

anage <- fread(anage.file, sep="\t", header=TRUE, 
               col.names=readLines(anage.hdr.file))
anage[, FullSpecies := str_c(Genus, Species, sep="_")]
setkey(anage, FullSpecies)

impute.dat <- fread(impute.file)
setkey(impute.dat, Species)

# Replace Zoonomia species with close timetree equivalent
nm.map <- fread(nm.map.file)
for(i in 1:nrow(nm.map)){
  sp.t <- nm.map[i, TimeTree]
  sp.z <- nm.map[i, Zoonomia]
  if(is.na(sp.t)){
    next
  }
  
  full.tree$tip.label[which(full.tree$tip.label == sp.t)] <- sp.z
}

# Filter the tree to match species in alpha table
tree.sp <- intersect(full.tree$tip.label, alpha$species)
tree.sp <- setdiff(tree.sp, nm.map[is.na(TimeTree), Zoonomia])

cat(str_interp("Pruning tree. Keeping ${length(tree.sp)} of ${length(alpha$species)} species in alpha table\n\n"),
    file=cat.file, append=TRUE)
tree <- keep.tip(full.tree, tree.sp)
tree[["node.label"]] <- str_c("A", Ntip(tree) + (1:Nnode(tree)))
alpha <- alpha[species %in% tree.sp,]

config <- fread(config.file)

##### IMPUTE ANAGE DATA IF NEEDED #####
alpha[, AnAge_presence := !is.na(AnAge_presence)]
alpha[, AnAgeSpecies := species]
n.impute <- alpha[, sum(!AnAge_presence)]
cat(str_interp("\nImputing AnAge data for ${n.impute} species:\n"), file=cat.file, append=TRUE)

imp.cols <- setdiff(intersect(names(alpha), names(anage)), c("Order", "Genus")) # Data to impute
for(i in alpha[, .I[!AnAge_presence]]){ # Iterate over rows to impute
  sp <- alpha[i, species]
  a.sp <- impute.dat[sp, anage_Species] # Species to take data from in anage dataset
  
  cat(str_interp("$[30s]{sp} --> ${a.sp}\n"), file=cat.file, append=TRUE)
  
  alpha[i, AnAgeSpecies := a.sp]
  anage.row <- anage[a.sp,]
  for(j in imp.cols){
    new.val <- anage.row[[j]]
    alpha[i, eval(j) := new.val]
  }
}


##### Functions #####
# Make output suffix based on orders kept or removed
# MakeSuffix <- function(ext, kp.order=NULL, rm.order=NULL, grp="thinned"){
#   kp.str <- NULL; rm.str <- NULL
#   if(!is.null(kp.order)){
#     kp.str <- str_c("kp_", str_c(kp.order, collapse="_"))
#   }
#   
#   if(!is.null(rm.order)){
#     rm.str <- str_c("rm_", str_c(rm.order, collapse="_"))
#   }
#   
#   out.str <- str_c(kp.str, rm.str, sep=".")
#   if(grp != "thinned"){
#     out.str <- str_c(out.str, grp, sep=".")
#   }
#   out.str <- str_c(out.str, ext, sep=".")
#   
#   return(out.str)
# }

# Prepare data frame for PGLS
prepareDataframe <- function(d, xvar, yvar, xlog, ylog){
  v <- c("species", "Order", xvar, yvar)
  subd <- d[,..v]
  colnames(subd) <- c("Species", "Order", "xvar", "yvar")
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
  cdat <- comparative.data(data = subd, phy = phylogeny, names.col = "Species")#,vcv=TRUE, vcv.dim=3
  mod <- pgls(yvar ~ xvar, cdat, lambda = lambda, kappa = kappa, delta = delta)
  
  df_mod = list("model" = mod, "data" = subd)
  return(df_mod)
}

# PIC pairwise
# Function to produce PICs based on pairwise comparison of tips only
# (no ancestral node reconstruciton; no overlaps).
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
  coef <- round(summary(mod)$coefficients[i,4], digits = 4)
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

## PCA functions
RotateVarimaxPCA <- function(traits, ncomp, rotate_method="varimax"){
  pca <- prcomp(traits, center=TRUE, scale=TRUE)
  
  sink(cat.file); print(summary(pca)); sink()
  
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
  pcs_rotated$species <- rownames(pcs)
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
  colnames(pca) <- c("species", "C1", "C2")
  if(latin2comm){
    pca$species <- latin.dat[pca$species, Common_names]
  }
  
  gg.obj <- ggplot(data=pca, aes(x=C1, y=C2)) + 
    geom_point() + geom_text_repel(aes(label=species)) +
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

# Extract data for particular comparison(s)
GetExperimentData <- function(xy_data, experiment){
  out.tab <- data.table()
  for(i in experiment){
    out.tab <- rbindlist(list(out.tab, xy_data[[i]]))
  }
  return(out.tab)
}

##### Add extra data columns #####
alpha[, Metabolic_rate_gram := Metabolic_rate / Adult_weight ]

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

chrom.level.sp <- alpha[AssemblyStatus == "Chromosome", species]
alpha[, c("Nearest_chrom_assembly", "Divergence_to_NCA") := MinDiv2Tips(ucsc.tree, species, chrom.level.sp), by=species]


##### Perform PCA, plot #####
# Choose PCA functions
pca_traits <- c("Adult_weight", "GenerationLength_d", "Gestation", "Birth_weight")  #
pca_dat <- as.data.frame(na.omit(alpha[, c("species", pca_traits), with=FALSE])) # drop NAs
cat(str_interp("Keeping ${nrow(pca_dat)} of ${nrow(alpha)} species with complete trait info.\n"),
    file=cat.file, append=TRUE)
row.names(pca_dat) <- pca_dat$species
pca_dat <- log(pca_dat[, pca_traits])

rotate_method <- "varimax"
# NB: Played around with quartimax rotation but it doesn't seem to differ all that
# much from the raw PCs

pca_and_varimax <- RotateVarimaxPCA(pca_dat, 2, rotate_method)
p <- fviz_pca_biplot(pca_and_varimax$pca, repel=TRUE)
print(p)

# Plot
pca_prop_var <- summary(pca_and_varimax$pca)$importance[2,]
p <- BiplotPCA(pca_and_varimax$pc_rc[,c("species","PC1","PC2")], pca_and_varimax$loadings) +
  xlab(str_interp("Component 1 ($[0.1f]{pca_prop_var[1]*100}%)")) +
  ylab(str_interp("Component 2 ($[0.1f]{pca_prop_var[2]*100}%)"))
print(p)
ggsave(p, filename=str_interp("${scripts.dir}/pgls_res/PCA.${exp.name}.pdf"), 
       width=6, height=6)

p <- BiplotPCA(pca_and_varimax$pc_rc[,c("species","RC1","RC2")], pca_and_varimax$loadings_rotated) +
  xlab(str_interp("${str_to_title(rotate_method)} Component 1")) +
  ylab(str_interp("${str_to_title(rotate_method)} Component 2"))
print(p)
ggsave(p, filename=str_interp("${scripts.dir}/pgls_res/PCA.${exp.name}.${rotate_method}_rotated.pdf"),
       width=6, height=6)

# Merge PCA with alpha table
alpha <- merge(alpha, pca_and_varimax$pc_rc, by.x='species', by.y="species", all.x=TRUE)

# Add PCs and RCs to configuration table
tmp <- data.table(x=c("PC1", "PC2", "RC1", "RC2"), y="alpha", xlog=FALSE, ylog=FALSE)
config <- rbindlist(list(config, tmp))

##### Run PGLS & PIC #####
xy_data <- list()
pic_data <- list()
models  <- c()

cor.method <- "spearman" # kendall or spearman

gnmstat.traits <- c("Divergence_to_NCA", "ContigN50", "ScaffoldN50", "AssemblyStatus")
pca_traits <- c(pca_traits, "PC1", "PC2")
keep.traits <- c(pca_traits, gnmstat.traits, "predicted_alpha")

# Convert AssemblyStatus to integers
alpha[, AssemblyStatus := as.factor(AssemblyStatus)]
config[, skipped := !(x %in% keep.traits)]

# Configure each trait variable (to be compared with alpha)
hdr <- c("Comparison", "pgls_ml_pvalue")
cat(str_interp("\n$[30s]{hdr[1]}   ${hdr[2]}\n"), file=cat.file, append=TRUE)
for (row in config[, .I[!skipped]]){
  # X and Y
  x <- config[row, x]
  y <- config[row, y]
  xlog <- config[row, xlog]
  ylog <- config[row, ylog]
  experiment <- x #config[row, experiment]
  
  # Run pgls and pic
  pgls_ml <- runPGLS(alpha, x, y, xlog, ylog, tree, lambda="ML", kappa=1, delta=1)
  if (pgls_ml[1]=="NA"){
    cat(str_interp("Skipping ${experiment}, no observations left"), file=cat.file, append=TRUE)
    config[row, skipped := TRUE]
    next
  }
  pgls    <- runPGLS(alpha, x, y, xlog, ylog, tree, lambda=1, kappa=1, delta=1)
  pic     <- picModel(alpha, x, y, xlog, ylog, tree)
  if(!(class(alpha[[x]]) %in% c("factor", "character"))){
    pic.cor <- picCor(alpha, x, y, xlog, ylog, tree, method=cor.method, pic.pairwise=FALSE)
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
  
  models <- rbind(models, c(pgls_ml_intercept, pgls_ml_slope, 
                            pgls_intercept, pgls_slope,
                            pgls_ml_pval, pgls_pval, pic_pval, pic_cor_pval,
                            pgls_ml_rsq, pgls_rsq,
                            ml_lambda, cor_rho))
}
rm(pgls_ml_pval, pgls_pval, pic_pval, pgls_ml_rsq, ml_lambda)

##### Process model summary stats #####
# Convert to data.table, add experiment code
colnames(models) <- c("pgls_ml_intercept", "pgls_ml_slope",
                      "pgls_intercept", "pgls_slope",
                      "pgls_ml_pval", "pgls_pval", "pic_pval", "pic_cor_pval",
                      "pgls_ml_rsq", "pgls_rsq",
                      "ml_lambda", "cor_rho")
models <- as.data.table(models)
models[, experiment := config[!(skipped), x]]

# Compute adjusted p-value by Benjamini-Hochberg
models[, pgls_ml_padj := p.adjust(pgls_ml_pval, "BH")]
models[, pgls_padj    := p.adjust(pgls_pval, "BH")]
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

# if(include_padj){
#   models[, pgls_ml_str   := sprintf("\n ML(%s)=%0.2f,  p [adj] = %0.4f [%0.4f]", lambda_char, ml_lambda, pgls_ml_pval, pgls_ml_padj), by=experiment]
#   models[, pgls_str    := sprintf("\n\n %s=1,  p [adj] = %0.4f [%0.4f]", lambda_char, pgls_pval, pgls_padj), by=experiment]
# } else {
  models[, pgls_ml_str := MakeModelStr(ml_lambda, pgls_ml_pval, pgls_ml_rsq,   "ML", prefix="\n"),     by=experiment]
  models[, pgls_str    := MakeModelStr(        1,    pgls_pval,    pgls_rsq, "full", prefix="\n\n"), by=experiment]
# }

# PIC stats
models[, cor_rho_str        := sprintf("\n %s %s = %0.2f", str_to_title(cor.method), cor_char[cor.method], cor_rho)]
models[, pic_cor_pval_str := sprintf("\n\n p = %0.4f", pic_cor_pval)]


##### Plot life history PGLS results #####
# ggplot theme
theme_set(theme_bw() +
            theme(axis.text=element_text(size=12), panel.border=element_rect(size = 1.5)))

# Line/text colors
c1 <- "#4363d8" # blue
c2 <- "#e6194B" # red

lightblue <- brewer.pal(3, "Paired")[1]
darkblue  <- brewer.pal(3, "Paired")[2]

purple <- brewer.pal("PRGn", n=3)[1]
green  <- brewer.pal("PRGn", n=3)[3]

w.factor <- 2.4 # width factor (inches)
h.factor <- 2.5 # height factor

# Function to plot PGLS results from given data, models
PlotPGLS <- function(xy_data, models, plt.experiments=names(xy_data), plt.type="PGLS", color.by="Order",
                     x.cis=NULL, y.cis=NULL){
  text.size <- 2.5
  
  nr <- round(sqrt(length(plt.experiments)))
  nc <- ceiling(length(plt.experiments) / nr)
  
  plt_data <- GetExperimentData(xy_data, plt.experiments)
  plt_data[, experiment := factor(experiment, levels=unique(experiment), labels=unique(experiment))]
  
  mod_data <- models[experiment %in% plt.experiments]
  mod_data[, experiment := factor(experiment, levels=experiment, labels=experiment)]
  
  if(length(plt.experiments) > 1){
    xlab.str <- "Covariate"
  } else {
    xlab.str <- plt.experiments
  }
  
  if(!is.null(x.cis)){ # x.cis must have column anmes Species, experiment, x_lwr, x_upr
    plt_data <- merge(plt_data, x.cis, by=c("Species", "experiment"))
  }
  
  if(!is.null(y.cis)){ # y.cis must have column anmes Species, experiment, y_lwr, y_upr
    plt_data <- merge(plt_data, y.cis, by=c("Species", "experiment"))
  }
  
  
  if(plt.type == "PGLS"){
    use.colors <- brewer.pal("Dark2", n=length(unique(plt_data$Order)))
    
    p <-ggplot(data = plt_data, aes(x=xvar, y=yvar))
    
    if(!is.null(x.cis)){
      p <- p + geom_errorbarh(aes(xmin=x_lwr, xmax=x_upr), alpha=0.75, col=lightblue, height=0)
    }
    
    if(!is.null(y.cis)){
      p <- p + geom_errorbar(aes(ymin=y_lwr, ymax=y_upr), alpha=0.75, col=lightblue, width=0)
    }
      
    if(!is.null(color.by) && (color.by != "")){ # Color points by Order (or something else)
      p <- p + geom_point(aes_string(col=color.by), alpha=0.9, size=1.2)
    } else {
      p <- p + geom_point(alpha=0.9, size=1.2, col=darkblue)#"#404040")
    }
    p <- p +
      geom_abline(aes(slope = pgls_ml_slope, intercept = pgls_ml_intercept), mod_data, linetype = "dashed", color=c1) + 
      geom_abline(aes(slope = pgls_slope, intercept = pgls_intercept), mod_data, linetype = "dashed", color=c2) + 
      geom_text(data = mod_data, color=c1, aes(label=pgls_ml_str, color="gray", x = -Inf, y = Inf),
                hjust = 0, vjust = 1, size = text.size) +
      geom_text(data = mod_data, color=c2, aes(label=pgls_str, color="gray", x = -Inf, y = Inf),
                hjust = 0, vjust = 1, size = text.size) +
      scale_color_manual(values=use.colors) +
      xlab(xlab.str) + ylab("Alpha")
  } else if(plt.type == "PIC") {
    p <-ggplot(data = plt_data, aes(PIC1, PIC2))
    if(!is.null(color.by) && (color.by != "")){ # Color points by Order (or something else)
      p <- p + geom_point(aes_string(col=color.by), alpha=0.8, size=1.2)
    } else {
      p <- p + geom_point(alpha=0.8, size=1.2, col=darkblue)
    }
    p <- p + 
      geom_text(data = mod_data, color="black", aes(label=cor_rho_str, color="gray", x = -Inf, y = Inf),
                hjust = 0, vjust = 1, size = text.size) +
      geom_text(data = mod_data, color="black", aes(label=pic_cor_pval_str, color="gray", x = -Inf, y = Inf),
                hjust = 0, vjust = 1, size = text.size) +
      xlab(str_interp("PIC (${xlab.str})")) + ylab("PIC (Alpha)")
  }
  
  if(length(plt.experiments) > 1){
    p <- p + facet_wrap(~experiment, nrow=nr, scales = "free_x", )
  }
  
  return(list("p"=p, "nr"=nr, "nc"=nc))
}

# Plot PGLS results for PCA traits + PCs only with both ML lambda and lambda=1 results
cur.plt <- PlotPGLS(xy_data, models, pca_traits, plt.type="PGLS")
print(cur.plt$p)
ggsave(cur.plt$p, filename=str_interp("${scripts.dir}/pgls_res/PGLS.${exp.name}.pca_traits.pdf"), 
       width=(w.factor*cur.plt$nc)+2, height=h.factor*cur.plt$nr, device=cairo_pdf)

# Plot PGLS results for genome quality traits
cur.plt <- PlotPGLS(xy_data, models, setdiff(gnmstat.traits, "AssemblyStatus"), plt.type="PGLS")
print(cur.plt$p)
ggsave(cur.plt$p, filename=str_interp("${scripts.dir}/pgls_res/PGLS.${exp.name}.gnmstat_traits.pdf"), 
       width=(w.factor*cur.plt$nc)+2, height=h.factor*cur.plt$nr, device=cairo_pdf)

# Plot PGLS results for predicted alpha vs alpha
x.cis <- alpha[, .(Species=species, experiment="predicted_alpha", 
                   x_lwr=predicted_alpha_lwr, x_upr=predicted_alpha_upr)]
y.cis <- alpha[, .(Species=species, experiment="predicted_alpha", 
                   y_lwr=alpha_lwr, y_upr=alpha_upr)]

cur.plt <- PlotPGLS(xy_data, models, "predicted_alpha", plt.type="PGLS", 
                    x.cis=x.cis, y.cis=y.cis, color.by=NULL)
print(cur.plt$p)
ggsave(cur.plt$p, filename=str_interp("${scripts.dir}/pgls_res/PGLS.${exp.name}.predicted_alpha.pdf"), 
       width=4, height=3, device=cairo_pdf)

##### Plot life history PIC results (non-parametric correlation) #####
cur.plt <- PlotPGLS(pic_data, models, setdiff(c(pca_traits, gnmstat.traits), "AssemblyStatus"), plt.type="PIC", color.by=NULL)
print(cur.plt$p)

ggsave(cur.plt$p, filename=str_interp("${scripts.dir}/pgls_res/PIC.${cor.method}.${exp.name}.all_traits.pdf"), 
       width=(w.factor*cur.plt$nc)+2, height=h.factor*cur.plt$nr, device=cairo_pdf)

##### Plot AssemblyStatus vs alpha results #####
if("AssemblyStatus" %in% names(xy_data)){
  plt_data <- GetExperimentData(xy_data, "AssemblyStatus")
  
  mod_data <- models[experiment == "AssemblyStatus"]
  
  hline_data <- mod_data[, .(ml_est = pgls_ml_intercept + c(0, pgls_ml_slope),
                             pic_est = pgls_intercept + c(0, pgls_slope),
                             x=1:2)]
  
  p <-ggplot(data = plt_data, aes(xvar, yvar)) + 
    geom_jitter(aes(col=Order), alpha=0.9, size=1.2, width=0.25) + 
    geom_violin(alpha=0.3) +
    geom_hpline(aes(x=x, y=ml_est), data=hline_data, width=0.2, size=1, col=c1, alpha=0.8) +
    geom_hpline(aes(x=x, y=pic_est), data=hline_data, width=0.2, size=1, col=c2, alpha=0.8) +
    geom_text(data = mod_data, color=c1, aes(label=pgls_ml_str, color="gray", x = -Inf, y = Inf),
              hjust = 0, vjust = 1, size = 2.5) +
    geom_text(data = mod_data, color=c2, aes(label=pgls_str, color="gray", x = -Inf, y = Inf),
              hjust = 0, vjust = 1, size = 2.5) +
    scale_color_manual(values=brewer.pal("Dark2", n=length(unique(plt_data$Order)))) +
    xlab("AssemblyStatus") + ylab("Alpha")
  print(p)
  
  ggsave(p, filename=str_interp("${scripts.dir}/pgls_res/PGLS.${exp.name}.AssemblyStatus.pdf"), 
         width=w.factor+2, height=h.factor, device=cairo_pdf)
}

##### Perform node height tests #####
x <- alpha[, alpha]
names(x) <- alpha$species
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

ggsave(p, filename=str_interp("${scripts.dir}/pgls_res/nheight.${exp.name}.pdf"), 
       width=4.5, height=4)

# geiger::tips(tree, 148)
# geiger::tips(tree, 150)


##### Test Divergence_to_NCA + GenerationLength_d #####
if(FALSE){
  subd <- alpha[, .(Species=species, 
                    experiment,
                    alpha, Divergence_to_NCA, GenerationLength_d)]
  for(i in c("Divergence_to_NCA", "GenerationLength_d")){
    if(config[x==i, xlog]){
      subd[, eval(i) := log10(get(i))]
    }
  }
    
  subd <- as.data.frame(na.omit(subd))
    
  if (nrow(subd)<3){
    stop("Can't analyze: alpha ~ Divergence_to_NCA + GenerationLength_d")
  }
  phylogeny <- keep.tip(tree, subd$Species)
  cdat <- comparative.data(data = subd, phy = phylogeny, names.col = "Species")
  mod <- pgls(alpha ~ Divergence_to_NCA + GenerationLength_d, cdat, lambda="ML", kappa=1, delta=1)
    
  df_mod <- list("model" = mod, "data" = subd)
}

