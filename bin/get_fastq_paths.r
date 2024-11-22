#!/usr/bin/env Rscript

# Author: @fwu19

options(stringsAsFactors = F)
library(dplyr)

## functions ####
get_fastqs <- function(fq_dirs){
    ## get paths to fastq files
    fqs <- normalizePath(grep('undetermined', list.files(fq_dirs, recursive = T, full.names = T, pattern = "fastq.gz"), invert = T, value = T, ignore.case = T))
    fqs1 <- sort(grep("_S[0-9]+(_L[0-9]+)?_R1_", fqs, value = T))
    fqs2 <- sort(grep("_S[0-9]+(_L[0-9]+)?_R2_", fqs, value = T))
    
    nfq1 <- length(fqs1)
    nfq2 <- length(fqs2)
    if (nfq1 == 0){
        stop ('No read 1 files found!')
    }
    
    if (nfq1 > 0 & nfq2 > 0 & nfq1 != nfq2){
        stop ('Found different numbers of read 1 and read 2 files!')
    }
    
    ## generate a sample sheet without metadata
    ss <- data.frame(
        fastq_1 = fqs1,
        single_end = ifelse(nfq2 == 0, 'true', 'false')
    ) %>% 
        mutate(
            id = gsub("_S[0-9]+(_L[0-9]+)?_R1_.*", "", basename(fastq_1)),
            fastq_2 = ifelse(single_end, "", fqs2)
        ) %>% 
        mutate(
            sample_group = id
        )
    
    ## 
    return(ss)
}

## read arguments ####
args <- as.vector(commandArgs(T))
input_dirs <- args

ss <- bind_rows(lapply(
    input_dirs, get_fastqs
))

ss %>% 
    write.table('input.csv', sep = ',', quote = F, row.names = F)

