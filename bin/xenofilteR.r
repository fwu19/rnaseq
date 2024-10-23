#!/usr/bin/env Rscript

## run xenofilteR

options(stringsAsFactors = F)
library(XenofilteR)

args <- commandArgs(T)
sample.list <- matrix(args[1:2],nrow=1)
output.names <- args[3]
mm.threshold <- as.integer(args[4])
nworkers <- as.integer(args[5])
  
XenofilteR(
  sample.list = sample.list,
  destination.folder = './',
  bp.param = SnowParam(workers = nworkers, type = 'SOCK'),
  output.names = output.names,
  MM_threshold = mm.threshold
)

