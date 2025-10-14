#!/usr/bin/env Rscript

options(stringsAsFactor=F)

library(rmarkdown)
library(pandoc)
args <- as.vector(commandArgs(T))

rmarkdown::render(args[1])
