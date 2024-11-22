#!/usr/bin/env Rscript

options(stringsAsFactor=F)
options(scipen = 99)

library(dplyr)

file.copy("report/rnaseq_regular.Rmd", "_Analysis_report.Rmd")
