#!/usr/bin/env Rscript
rm(list=ls())

##### LIBRARIES #####
#library(ggplot2)
library(data.table)
library(stringr)
library(ape)
#library(RColorBrewer)

##### SET UP ENV #####
scripts.dir <- "/moto/palab/projects/male_mutation_bias_XA/mut_sex_bias_amniotes/scripts"
dir <- "/moto/palab/projects/male_mutation_bias_XA/mut_sex_bias_amniotes/"
#scripts.dir <- "/Users/felix/mt_mp_lab/projects/male_mutation_bias_XA/scripts"

# ggplot theme
#theme_set(theme_bw() +
#            theme(axis.text=element_text(size=12), panel.border=element_rect(size = 1.5)))

##### FUNCTIONS #####
GetMLFile <- function(template.str, n.reps=6){
  # template.str should have ${i} somewhere in the string to denote the rep number
  ml.vals <- c()
  for(i in 1:n.reps){
    cur.fn <- str_interp(template.str)
    ml.vals[i] <- fread(cmd=str_interp("grep TRAINING_LNL ${cur.fn}"), header=FALSE)$V2
  }
  
  i <- which.max(ml.vals)
  return(str_interp(template.str))
}

DistTips2Parent <- function(phy){
  dist.mat <- dist.nodes(phy) # Node distance matrix
  path.list <- nodepath(phy) # List of root to tip node paths
  
  # Get distance from tip to parent
  out.vals <- lapply(path.list, function(x){y <- tail(x, n=2); dist.mat[y[1], y[2]]})
  
  names(out.vals) <- phy$tip.label
  return(out.vals)
}

##### GET ORDER DATA #####
outgroup.tab <- fread(str_interp("${dir}/data/outgroups_to_discard_orders.tab"))
setkey(outgroup.tab, experiment)

#param.tab <- fread(str_interp("${dir}/data/SuperOrders_to_assemblies-X.tab"), 
#                   col.names = c("experiment", "ref", "X_name"), header=FALSE)
#setkey(param.tab, experiment)

#super.orders <- param.tab$experiment
#super.orders <- str_subset(super.orders, "Control", negate=TRUE)

##### LOAD TREES #####
ucsc.phy <- read.tree(str_interp("${dir}/trees/241-mammalian-2020v2.phast-242.nh"))

time.phy <- read.tree(str_interp("${dir}/trees/mammals241.TimeTree.nwk"))
time.map <- fread(str_interp("${dir}/data/mammals241.TimeTree.spname_map.tab"))
setkey(time.map, Zoonomia)
time.unmap.sp <- time.map[is.na(TimeTree), Zoonomia]

#noMnC.sp.list <- list()
#thinned.sp <- c(); ultrathinned.sp <- c()
#for(so in super.orders){
#  
#  tmp <- fread(str_interp("${scripts.dir}/species_list/${so}_noMale_noComplex.txt"), sep=",", header=FALSE)
#  noMnC.sp.list[[so]] <- setdiff(unname(unlist(tmp)), c(outgroup.tab[so, species], time.unmap.sp))
#}

##### FUNCTIONS #####
# Drop outgroups and species missing from TimeTree if need be
GetSpeciesNoOutgroup <- function(fn, so, outgroup.tab, time.unmap.sp){
  all.sp <- unname(unlist(fread(fn, sep=",", header=FALSE)))
  keep.sp <- setdiff(all.sp, c(outgroup.tab[so, species], time.unmap.sp))
  return(keep.sp)
}

AppendSpecies <- function(dat, keep.sp, so, subset.nm){
  tmp.dat <- data.table(Species=keep.sp, SuperOrder=so, Subset=subset.nm)
  return(rbindlist(list(dat, tmp.dat)))
}


##### BUILD SUBRATE TABLE #####

# Mammals
out.dat <- data.table()
exp.str <- "AllMammals"
fn <- str_interp("${dir}/data/${exp.str}.txt")
keep.sp <- GetSpeciesNoOutgroup(fn, exp.str, outgroup.tab, time.unmap.sp)
out.dat <- AppendSpecies(out.dat, keep.sp, exp.str, exp.str)

# UCSC rates
out.dat[, SubrateUCSC := unlist(DistTips2Parent(keep.tip(ucsc.phy, Species)))[Species], by=Subset]

# TimeTree species names
out.dat[!(Species %in% time.map$Zoonomia), TimeTreeSpecies := Species]
for(sp in time.map$Zoonomia){
  out.dat[Species == sp, TimeTreeSpecies := time.map[sp, TimeTree]]
}

# TimeTree times
out.dat[, DivTime := unlist(DistTips2Parent(keep.tip(time.phy, TimeTreeSpecies)))[TimeTreeSpecies], by=Subset]

##### USING FULL 241 MAMMAL DATASET #####
keep.sp <- setdiff(ucsc.phy$tip.label, time.unmap.sp)
sub.phy <- keep.tip(ucsc.phy, keep.sp)

cur.dat <- out.dat[Subset == "noMale_noComplex",]
# cur.dat[, SubratePhylofit := NA]
cur.dat[, SubrateUCSC := unlist(DistTips2Parent(sub.phy))[Species]]
cur.dat[, DivTime := unlist(DistTips2Parent(time.phy))[TimeTreeSpecies]]
cur.dat[, Subset := "full"]

out.dat <- rbindlist(list(out.dat, cur.dat))

#### YEARLY RATE #####
# out.dat[, MutPerYearPhylofit := SubratePhylofit / DivTime]
out.dat[, MutPerYearUCSC     := SubrateUCSC     / DivTime]

#### WRITE RESULTS #####
fwrite(out.dat, file=str_interp("${dir}/data/subrate_div_data_allMammals.txt"),
       sep="\t", col.names=TRUE)


##### PLOT #####
#plt.dat <- dcast(out.dat[Subset %in% c("thinned", "full")], Species + SuperOrder ~ Subset, value.var="MutPerYearUCSC")
#plt.dat <- plt.dat[!is.na(thinned)]
#gg.obj <- ggplot(aes(x=full, y=thinned, color=SuperOrder), data=plt.dat) + 
#  geom_point(alpha=0.6) #+ scale_color_manual(values=brewer.pal(length(super.orders), "Set3"))
#print(gg.obj)


# Violin plots of rates
#gg.obj <- ggplot(aes(x=MutPerYearUCSC, y=SuperOrder), data=out.dat) + 
#  geom_violin(fill="gray") + facet_wrap(vars(Subset), nrow=1)
#print(gg.obj)

##### SAVE TREES #####
#for(s in c("Mammals")){
#  tmp.dat <- out.dat[Subset == s]
#  out.phy <- keep.tip(time.phy, tmp.dat[, TimeTreeSpecies])
  
#  setkey(tmp.dat, TimeTreeSpecies)
#  out.phy$tip.label <- tmp.dat[(out.phy$tip.label), Species]
  
#  write.tree(out.phy, file=str_interp("${scripts.dir}/pgls_files/trees/TimeTree.${s}.nwk"))
#}
