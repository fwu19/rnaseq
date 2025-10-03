#!/usr/bin/env Rscript

# Author: @fwu19

options(stringsAsFactors = F)
options(scipen = 99)
library(dplyr)

## read arguments ####
args <- as.vector(commandArgs(T))
for (arg in strsplit(args, split = '=')){
    assign(arg[1], arg[2])
}

stopifnot(exists('csv') & exists('outdir'))
if(!grepl('dummy', csv) & !dir.exists(file.path(outdir, 'csv'))){
    stop("No valid input found in --input or params.outdir/csv/")
}

if (step %in% c('expression_quantification')){
    if (grepl('dummy', csv) & file.exists(file.path(outdir, 'csv/align_fastq.csv'))){
        csv <- file.path(outdir, 'csv/align_fastq.csv')
    }
    ss <- read.csv(in_csv)
    if(!'bam' %in% colnames(ss)){
        stop(paste("Column bam is not found in", in_csv))
    }else if (sum(!file.exists(ss$bam))>0){
        stop(paste("Some bam files in", in_csv, "do not exist!"))
    }
    file.copy(in_csv, 'samplesheet.checked.csv')
}


