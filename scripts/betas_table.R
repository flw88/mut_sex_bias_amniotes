#!/usr/bin/env Rscript
rm(list=ls())

##### LIBRARIES #####

suppressMessages( library(getopt) )
library(data.table)
library(stringr)
library(ggplot2)

theme_set(theme_bw() +
            theme(axis.text=element_text(size=12), panel.border=element_rect(size = 1.5)))

##### GET OPTS, RUN PARAMETERS #####
#opt.spec <- matrix(c("help",         "h",  0,  "logical",   0,
#                     "exp.name",     "c",  1,  "character", 1,
#                     "ref.species",  "r",  1,  "character", 1
#), byrow=TRUE, ncol=5)
#opt <- getopt(opt.spec[,1:4])
#req.args <- opt.spec[as.logical(as.integer(opt.spec[,5])), 1]

main.dir <- "/moto/palab/projects/male_mutation_bias_XA/mut_sex_bias_amniotes"

###--- Run this commented block to test ---#
# main.dir <- "/Users/felix/mt_mp_lab/projects/male_mutation_bias_XA/mut_sex_bias_amniotes"
#opt <- list(exp.name="Mammals", ref.species="Homo_sapiens")
###--------------------------------------------#

scripts.dir <- str_interp("${main.dir}/scripts")
data.dir    <- str_interp("${main.dir}/data")
pgls.dir    <- str_interp("${main.dir}/pgls_files")

source(str_interp("${scripts.dir}/alpha_regression.funcs.R"))

# if( !is.null(opt$help) || (length(opt) == 0) ) {
#   cat(getopt(opt.spec, usage=TRUE))
#   quit(status=1)
# }

# for(arg in req.args){
#   if( is.null(opt[[arg]]) ){
#     cat(getopt(opt.spec, usage=TRUE))
#     stop(str_interp("${arg} is a required argument."))
#   }
# }


grp.and.ref <- c("Mammals"="Homo_sapiens",
                 "Birds"="Gallus_gallus",
                 "Snakes"="Thamnophis_elegans")

n.digits <- 3 # how many digits to round to

##### LOAD DATA ######

ReportBetas <- function(grp, ref, mt="mod"){
  if(grp %in% c("Mammals", "Snakes")){
    load(str_interp("${scripts.dir}/lm_res/${grp}.${ref}.RData"))
  } else if(grp == "Birds"){
    lm.tmp <- list()
    sc.tmp <- list()
    for(i in 1:6){
      load(str_interp("${scripts.dir}/lm_res/${grp}_g${i}.${ref}.RData"))
      if(i == 1){
        for(mt in names(lm.list)){
          lm.tmp[[mt]] <- list()
          sc.tmp[[mt]] <- list()
        }
      }
      for(mt in names(lm.list)){
        for(sp in names(lm.list[[mt]])){
          if((i != 6) && (sp == "Gallus_gallus")){ next }
          lm.tmp[[mt]][[sp]] <- lm.list[[mt]][[sp]]
          sc.tmp[[mt]][[sp]] <- samp.coef.list[[mt]][[sp]]
          }
      }
    }
    lm.list <- lm.tmp
    samp.coef.list <- sc.tmp
    rm(lm.tmp, sc.tmp)
  }

  # Compile beta values
  betas <- list()
  for(sp in names(lm.list[[mt]])){
    est <- round(coef(lm.list[[mt]][[sp]]), n.digits)
    cis <- round(apply(samp.coef.list[[mt]][[sp]], 2, ConfInt), n.digits) # CIs
    
    betas[[sp]] <- data.table("Species"      = sp,
                              "Group"        = grp,
                              "beta0"        = str_c(unname(est["(Intercept)"])),
                              "beta0_95%_CI" = str_c("(", str_c(cis[c("lwr","upr"),"(Intercept)"], collapse=","), ")"),
                              "beta1"        = str_c(unname(est["regionXZ"])),
                              "beta1_95%_CI" = str_c("(", str_c(cis[c("lwr","upr"),"regionXZ"], collapse=","), ")"),
                              "beta2"        = str_c(unname(est["gc"])),
                              "beta2_95%_CI" = str_c("(", str_c(cis[c("lwr","upr"),"gc"], collapse=","), ")"),
                              "beta3"        = str_c(unname(est["I(gc^2)"])),
                              "beta3_95%_CI" = str_c("(", str_c(cis[c("lwr","upr"),"I(gc^2)"], collapse=","), ")")
                              )
  }

  return(rbindlist(betas))
}

# Compile table
beta.tab <- list()
for(grp in names(grp.and.ref)){
  ref <- grp.and.ref[grp]
  beta.tab[[grp]] <- ReportBetas(grp, ref)
}
beta.tab <- rbindlist(beta.tab)

# Write results
out.fn <- str_interp("${data.dir}/Table_S6.tsv")
fwrite(beta.tab, file=out.fn, quote=FALSE, sep="\t", col.names=TRUE, row.names=FALSE)

