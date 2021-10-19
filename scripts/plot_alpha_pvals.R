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


##### SETUP #####
# ggplot theme
theme_set(theme_bw() +
            theme(axis.text=element_text(size=12), panel.border=element_rect(size = 1.5)))

##### LOAD DATA #####
