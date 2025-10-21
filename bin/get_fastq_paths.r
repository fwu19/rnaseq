#!/usr/bin/env Rscript

# Author: @fwu19

options(stringsAsFactors = F)
library(dplyr)

## functions ####
get_fastqs <- function(
    fq_dirs, 
    r1_pattern = "_S[0-9]+(_L[0-9]+)?_R1_",
    r2_pattern = "_S[0-9]+(_L[0-9]+)?_R2_"
){
    ## get paths to fastq files
    fqs <- normalizePath(grep('undetermined', list.files(fq_dirs, recursive = T, full.names = T, pattern = "fastq.gz"), invert = T, value = T, ignore.case = T))
    fqs1 <- sort(grep(r1_pattern, fqs, value = T))
    fqs2 <- sort(grep(r2_pattern, fqs, value = T))
    
    nfq1 <- length(fqs1)
    nfq2 <- length(fqs2)
    if (nfq1 == 0){
        stop ('No read 1 files found!')
    }
    
    # if fastq_2 exist, should have the same number as fastq_1
    if (nfq1 > 0 & nfq2 > 0 & nfq1 != nfq2){
        stop ('Found different numbers of read 1 and read 2 files!')
    }
    
    ## generate a sample sheet without metadata
    ss <- data.frame(
        id = gsub(paste0(r1_pattern, ".*"), "", basename(fqs1)),
        fastq_1 = fqs1
    )
    
    if (nfq2 > 0){
        ss <- ss %>% 
            mutate(
                fastq_2 = fqs2,
                single_end = 'false'
            )
    } else{
        ss <- ss %>% 
            mutate(
                fastq_2 = "",
                single_end = 'true'
            )
    }
    
    ss <- ss %>% 
        mutate(
            id = gsub(' +|&', '-', id)
        )
    
    ## 
    return(ss)
}

## read arguments ####
args <- as.vector(commandArgs(T))
lst <- strsplit(args, split = '=')
for (x in lst){
    assign(x[1],x[2])
} # read arguments: r1_pattern r2_pattern
rm(lst)

if(!exists('r1_pattern')){r1_pattern <- "_S[0-9]+(_L[0-9]+)?_R1_"}
if(!exists('r2_pattern')){r2_pattern <- "_S[0-9]+(_L[0-9]+)?_R2_"}

fqs <- list.files('fastq', full.names = T)
if(file_test('-f', fqs[1])){
    input_dirs <- scan(fqs[1], what = 'character')
}else if (file_test('-d', fqs[1])){
    input_dirs <- fqs[1]
}else{
    stop("--input_dir takes either path/to/fastq/dir or a file containing paths/to/fastq/dir (one path in each row)")
}
ss <- get_fastqs(input_dirs)

ss %>% 
    write.table('input.csv', sep = ',', quote = F, row.names = F)

