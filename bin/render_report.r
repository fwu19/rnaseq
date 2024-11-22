#!/usr/bin/env Rscript

Sys.setenv(RSTUDIO_PANDOC="/app/software/RStudio-Server/1.4.1717-foss-2021b-Java-11-R-4.1.2/bin/pandoc")

options(stringsAsFactor=F)

library(rmarkdown)
library(pandoc)
args <- as.vector(commandArgs(T))

rmarkdown::render(args[1])
