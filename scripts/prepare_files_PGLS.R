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
                     "out.exp",      "o",  1,  "character", 1
), byrow=TRUE, ncol=5)
opt <- getopt(opt.spec[,1:4])
req.args <- opt.spec[as.logical(as.integer(opt.spec[,5])), 1]

main.dir <- "/moto/palab/projects/male_mutation_bias_XA/mut_sex_bias_amniotes"

###--- Run this commented block to test locally ---#
main.dir <- "/Users/felix/mt_mp_lab/projects/male_mutation_bias_XA/mut_sex_bias_amniotes"
opt <- list(ref.species="Homo_sapiens",
            experiment.name="Mammals",
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
experiment.name <- opt$experiment.name
out.exp <- opt$out.exp # Output experiments. Sub groupings of a particular order
out.exp <- str_split(out.exp, ",")[[1]]

##### Set up file names #####
if(experiment.name != "Mammals"){
  stop("Currently can only handle 'Mammals' right now")
}

replace.file     <- str_interp("${pgls.dir}/replace_species.txt")
output.file      <- str_interp("${pgls.dir}/${experiment.name}.csv")
tree.file        <- str_interp("${main.dir}/trees/241-mammalian-2020v2.phast-242.nh")
alpha.file       <- str_interp("${main.dir}/data/XA_2exposure-model_alphas.csv")
sr.file          <- str_interp("${pgls.dir}/subrate_div_data.txt")

# anage.file       <- str_interp("${pgls.dir}/anage_data.txt")
# anage.hdr.file   <- str_interp("${pgls.dir}/anage_simplified_header.txt")
# qual.file        <- str_interp("${main.dir}/data/zoonomia_assembly_metadata.csv")
# gl.file          <- str_interp("${main.dir}/data/XA_2exposure-model_alphas.csv")

##### Load alpha, tree data #####
# Read alpha estimates & nwk tree
alpha.dat <- fread(alpha.file)

ucsc.tree <- keep.tip(read.tree(tree.file), c(alpha.dat$species))

##### Load genome assembly quality data #####
# qual <- fread(qual.file)
# qual.cols <- c("Species", "AssemblyStatus", "Coverage", "ContigN50", "ScaffoldN50")

##### Load life history data #####
# # AnAge database
# anage <- fread(anage.file, sep="\t", header = TRUE, col.names=readLines(anage.hdr.file))
# anage[, FullSpecies := str_c(Genus, Species, sep="_")]

# Generation time and predicted alpha
# gl <- fread(gl.file, sep=",")
# setnames(gl, "Species", "FullSpecies")
# gl[, c("V1", "Order", "Family", "Genus") := NULL] # Delete unneeded columns

##### Load sub rate data file #####
# Subrate
# sr <- fread(sr.file)
# setnames(sr, "Species", "FullSpecies")
# sr <- sr[Subset=="thinned", .(FullSpecies, MutPerYearUCSC)]
# sr[, mutrate_yearly := MutPerYearUCSC / 1e6]

##### Fill in missing info #####
# # If species are missing use closely related species data (if present)
# replace.species <- fread(replace.file, header=TRUE, na.strings="NaN")
# 
# # Merge alpha and replace.species tables
# alpha.traits <-  merge(alpha.dat, replace.species, by="Species", all.x=TRUE)
# setnames(alpha.traits, "Species", "FullSpecies")
# 
# # Merge alpha and quality data based on species names
# alpha.traits <- merge(alpha.traits, qual[, qual.cols, with=FALSE], by.x="FullSpecies", by.y="Species", all.x=TRUE)
# 
# # Merge alpha and AnAge data based on species names
# alpha.traits <- merge(alpha.traits, anage, by.x="anage_Species", by.y="FullSpecies", all.x=TRUE)
# 
# # Merge alpha and subrate data based on species names
# alpha.traits <- merge(alpha.traits, sr, by="FullSpecies", all.x=TRUE)
# 
# # Merge alpha and generation length data based on species names
# alpha.traits <- merge(alpha.traits, gl, by.x="gl_Species", by.y="FullSpecies", all.x=TRUE)
# 
# setkey(alpha.traits, FullSpecies)

# Add some remaining missing info for Homo Sapiens (30y gen time)
# alpha.traits["Homo_sapiens", GenerationLength_d := 30*365]

# Add metabolic rate per gram
# alpha.traits[, Metabolic_rate_gram := Metabolic_rate / Adult_weight]

##### Report how many species required values from closely related species #####
for(i in c("anage", "gl")){
  cur.col <- str_c(i, "_Species")
  n.replace <- alpha.traits[!is.na(get(cur.col)), sum(FullSpecies != get(cur.col))]
  n.miss <- alpha.traits[, sum(is.na(get(cur.col)))]
  cat(str_interp("${i} dataset: ${n.replace} species replaced, ${n.miss} missing\n"))
}
x <- alpha.traits[, sum(is.na(anage_Species) & is.na(gl_Species))]
cat(str_interp("${x} species lack any trait data\n"))


##### Compute divergence from nearest species in order, nearest chromosome level assembly #####
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

chrom.level.sp <- alpha.traits[AssemblyStatus == "Chromosome", FullSpecies]
alpha.traits[, c("Nearest_chrom_assembly", "Divergence_to_NCA") := MinDiv2Tips(ucsc.tree, FullSpecies, chrom.level.sp), by=FullSpecies]


##### Write table to output #####
# Output file
fwrite(alpha.dat, file=output.file, sep=",")



