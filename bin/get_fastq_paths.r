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
    
    # if (nfq1 > 0 & nfq2 > 0 & nfq1 != nfq2){
    #     stop ('Found different numbers of read 1 and read 2 files!')
    # }
    
    ## generate a sample sheet without metadata
    ss <- data.frame(
        id = gsub("_S[0-9]+(_L[0-9]+)?_R1_.*", "", basename(fqs1)),
        fastq_1 = fqs1
    )
    
    if (nfq2 > 0){
        ss <- ss %>% 
            mutate(
                match_id = gsub("_S[0-9]+(_L[0-9]+)?_R1_.*", "", fastq_1)
            ) %>% 
            left_join(
                data.frame(
                    fastq_2 = fqs2,
                    match_id = gsub("_S[0-9]+(_L[0-9]+)?_R2_.*", "", fqs2)
                ),
                by = 'match_id'
            ) %>% 
            mutate(
                single_end = ifelse(is.na(fastq_2), 'true', 'false'),
                fastq_2 = ifelse(is.na(fastq_2), "", fastq_2)
            ) %>% 
            dplyr::select(-match_id)
    } else{
        ss <- ss %>% 
            mutate(
                fastq_2 = "",
                single_end = 'true'
            )
    }
    
    ss <- ss %>% 
        mutate(
            sample_group = id
        )
    
    ## 
    return(ss)
}

## read arguments ####
args <- as.vector(commandArgs(T))

if(file_test('-f', args[1])){
    input_dirs <- scan(args[1], what = 'character')
}else if (file_test('-d', args[1])){
    input_dirs <- args[1]
}else{
    stop("--input_dir takes either path/to/fastq/dir or a file containing paths/to/fastq/dir (one path in each row)")
}
ss <- get_fastqs(input_dirs)

ss %>% 
    write.table('input.csv', sep = ',', quote = F, row.names = F)

