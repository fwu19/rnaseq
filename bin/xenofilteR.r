#!/usr/bin/env Rscript

## run xenofilteR

options(stringsAsFactors = F)
library(XenofilteR)

# graft_dir=graft host_dir=host out_prefix=${out_prefix} mm_threshold=$mm_threshold nworkers=${task.cpus}
args <- commandArgs(T)
for (arg in strsplit(args, split = '=')){
  assign(trimws(arg[1]), trimws(arg[2]))
}

graft_bams <- list.files(graft, pattern = 'bam', full.names = T)
host_bams <- list.files(host, pattern = 'bam', full.names = T)

XenofilteR(
  sample.list = cbind(graft_bams, host_bams),
  destination.folder = './',
  bp.param = SnowParam(workers = nworkers, type = 'SOCK'),
  output.names = out_prefix,
  MM_threshold = mm_threshold
)

