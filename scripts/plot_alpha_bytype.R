#!/usr/bin/env Rscript
rm(list=ls())

library(data.table)
library(ggplot2)
library(RColorBrewer)
library(stringr)

theme_set(theme_bw() +
            theme(axis.text    = element_text(size=12), 
                  panel.border = element_rect(size = 1.5)))

##### LOAD DATA #####
# scripts.dir <- "/moto/palab/projects/male_mutation_bias_XA/mut_sex_bias_amniotes/scripts"
work.dir <- "/Users/felix/mt_mp_lab/projects/male_mutation_bias_XA/mut_sex_bias_amniotes"
scripts.dir <- str_interp("${work.dir}/scripts")

cur.fn <- str_interp("${scripts.dir}/alphas/Mammals.Homo_sapiens.LM.tsv")
a.est <- fread(cur.fn, header=TRUE)
hum.est <- a.est[(species == 'Homo_sapiens') & !(mut_type %in% c("mod", "exptotsub"))]

cur.fn <- str_interp("${work.dir}/data/decode.alpha_by_type.tsv")
a.type <- fread(cur.fn, header=TRUE)

##### COMPARE HUMAN ESTIMATES BY TYPE #####
plt.dat <- merge(a.type, hum.est, by.x="Type", by.y="mut_type")
# plt.dat <- plt.dat[Category != "Mut type"]
setorder(plt.dat, Category)
plt.dat[, Type := factor(Type, levels=Type)]

plt.lim <- c(3.2, 8.2)

use.colors <- brewer.pal("Paired", n=nrow(plt.dat))
if(length(use.colors) >= 11){
  use.colors[11] <- "#BCBCBC"
}

if(length(use.colors) >= 12){
  use.colors[12] <- "#4D4D4D"
}

p <- ggplot(plt.dat, aes(x=Alpha, y=alpha, color=Type)) + 
  geom_abline(intercept=0, slope=1, color="gray", alpha=0.5) + 
  geom_errorbarh(aes(xmin=Alpha_lwr, xmax=Alpha_upr), height=0) +
  geom_errorbar( aes(ymin=alpha_lwr, ymax=alpha_upr), width=0) +
  geom_point() + 
  scale_color_manual(values=use.colors) +
  xlab("Alpha (DNMs)") + ylab("Alpha (sub rates)") + 
  coord_cartesian(xlim=plt.lim, ylim=plt.lim) + facet_wrap(~ Category, nrow=2)
print(p)

out.fn <- str_interp("${scripts.dir}/pdfs/Homo_sapiens.alpha_bytype.pdf")
ggsave(p, file=out.fn, width=6, height=4)


