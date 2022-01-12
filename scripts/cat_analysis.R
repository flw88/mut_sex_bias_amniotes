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
suppressMessages( library(flexmix) )

scripts.dir <- "/moto/palab/projects/male_mutation_bias_XA/mut_sex_bias_amniotes/scripts"

theme_set(theme_bw() +
            theme(axis.text=element_text(size=12), panel.border=element_rect(size = 1.5)))

keep.species <- c("Felis_catus", "Bos_taurus", "Homo_sapiens", "Dicerorhinus_sumatrensis")

##### LOAD DATA #####
load(str_interp("${scripts.dir}/lm_res/Mammals.Homo_sapiens.RData"))
clust.dat <- fread(str_interp("${scripts.dir}/mix_res/Mammals.Homo_sapiens.gc2.plot_data.tsv"), 
                   sep="\t", na.strings="NA", header=TRUE)[Species %in% keep.species]

##### PREP DATA #####
x.dat <- dat[(region == "XZ"),]
for(sp in keep.species){
  # Which cluster is the max one?
  hi.clust <- clust.dat[Species == sp, mean(sub_rate), by=.(clust)][which.max(V1), clust]
  lo.clust <- clust.dat[Species == sp, mean(sub_rate), by=.(clust)][which.min(V1), clust]
  
  # Which windows are max clust windows?
  hi.w <- clust.dat[(Species == sp) & (clust == hi.clust), start]
  lo.w <- clust.dat[(Species == sp) & (clust == lo.clust), start]
  
  new.col <- str_interp("${sp}_clust")
  x.dat[start %in% hi.w, eval(new.col)  := "hi"]
  x.dat[start %in% lo.w, eval(new.col)  := "lo"]
}

##### FUNCTIONS #####
GetRatio <- function(x.dat, mut.col, clust.col, ci.int=c(0.025, 0.975)){
  num.col <- str_interp("${mut.col}_num")
  den.col <- str_interp("${mut.col}_den")
  
  tmp <- x.dat[, sum(get(num.col)) / sum(get(den.col)), by=clust.col]
  
  out.list <- list(hi_clust_mean = tmp[get(clust.col) == "hi",V1],
                   lo_clust_mean = tmp[get(clust.col) == "lo",V1])
  
  max.samp <- replicate(sample(x.dat[get(clust.col) == "hi", INDEX], replace=TRUE), n=250)
  min.samp <- replicate(sample(x.dat[get(clust.col) == "lo", INDEX], replace=TRUE), n=250)
  
  rat.samp <- apply(max.samp, 2, function(x) x.dat[INDEX %in% x, sum(get(num.col)) / sum(get(den.col))]) / 
    apply(min.samp, 2, function(x) x.dat[INDEX %in% x, sum(get(num.col)) / sum(get(den.col))])
  tmp <- quantile(rat.samp, ci.int)
  
  out.list[["ratio"]] <- out.list$hi_clust_mean / out.list$lo_clust_mean
  out.list[["ratio_lwr"]] <- tmp[1]
  out.list[["ratio_upr"]] <- tmp[2]
  
  out.list[["sub_rate_mean"]] <- x.dat[, sum(get(num.col)) / sum(get(den.col))]
  
  tmp <- x.dat[, sum(size), by=clust.col]
  out.list[["hi_size"]] <- tmp[get(clust.col) == "hi", V1]
  out.list[["lo_size"]] <- tmp[get(clust.col) == "lo", V1]
  
  return(out.list)
}

##### ANALYZE MUTATION TYPES #####
mut.types <- names(lm.list)
type.dat <- list()
for(sp in keep.species){
  clust.col <- str_interp("${sp}_clust")
  for(i in mut.types){
    cat(i, "\n", sep="")
    mut.col <- str_interp("${sp}_${i}")
    
    cur.dat <- as.data.table(GetRatio(x.dat, mut.col, clust.col))
    cur.dat[, mut_type := i]
    cur.dat[, Species := sp]
    
    type.dat[[str_c(sp, ".", i)]] <- cur.dat
  }
}
type.dat <- rbindlist(type.dat)

type.dat[, mut_group := "all"]
type.dat[str_detect(mut_type, "C>|T>"), mut_group := "single nt"]
type.dat[str_detect(mut_type, "S>|W>"), mut_group := "strong/weak"]
type.dat[str_detect(mut_type, "BGC"), mut_group := "BGC"]

setorder(type.dat, mut_group)
type.dat[, mut_type := factor(mut_type, levels=unique(mut_type))]

# Plot
type.dat[, Label := sprintf("%s (%0.2f Mb, %0.2f Mb)", Species[1], hi_size[1] / 1e6, lo_size[1] / 1e6), by=Species]
p <- ggplot(aes(x=mut_type, y=ratio, ymin=ratio_lwr, ymax=ratio_upr), data=type.dat) + 
  geom_linerange(aes(color=mut_group)) + 
  geom_point(aes(color=mut_group, size=sub_rate_mean), fill="white", shape=21, stroke=1) +
  geom_hline(yintercept = 1, color="red", alpha=0.5) + 
  theme(axis.text.x = element_text(angle=-30, hjust=0), panel.grid.major.x=element_blank()) + 
  xlab("") + ylab("Sub Rate Ratio (High clust / Low clust)") + 
  scale_size_continuous(range = c(1, 5)) + 
  facet_wrap(~ Label, nrow=length(keep.species)) 

print(p)

ggsave(p, file=str_interp("${scripts.dir}/mix_res/cat_analysis.pdf"), 
       width=10, height=4*length(keep.species))




