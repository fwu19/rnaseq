#!/usr/bin/env Rscript

# Author: @fwu19

options(stringsAsFactors = F)
options(scipen = 99)
library(tidyverse)

## functions ####
add_col <- function(df, col, default.value){
    if(!col %in% colnames(df)){
        df[,col] <- default.value
    }
    return(df[,col])
}

add_metadata <- function(ss, meta_csv){
    ## update with metadata if provided
    # metadata contains a required columns id and optional columns: sample_group, sample_replicate, target, control, call_peak, call_rep_peak, call_con_peak
    
    if (file_test('-f', meta_csv) & grepl('.csv$', meta_csv) & !grepl('dummy', meta_csv)){
        meta <- read.csv(meta_csv)
        ss <- ss %>% 
            inner_join(
                meta, by = 'id', suffix = c(".x", "")
            ) %>% 
            dplyr::select(!ends_with(".x")) 
    }
    
    ## add missing columns
    ss$sample_group <- add_col(ss, 'sample_group', ss$id)
    
    if (sum(duplicated(ss$id)) > 0){
        ss <- ss %>% 
            group_by(id) %>% 
            mutate(
                fastq_1 = paste(basename(fastq_1), collapse = ';'),
                fastq_2 = paste(basename(fastq_2), collapse = ';')
            ) %>% 
            unique.data.frame()
    }
    ss %>% 
        dplyr::relocate(id, fastq_1, fastq_2, sample_group)
        
}

## read arguments ####
args <- as.vector(commandArgs(T))
in_csv <- args[1]
meta_csv <- ifelse(length(args) > 1, args[2], '')


## generate sample sheet ####
ss <- read.csv(in_csv)

## check single end fastq
if (!'id' %in% colnames(ss)){
    stop ( 'Missing column id!' )
}else if ( ! 'fastq_1' %in% colnames(ss)){
    stop ( 'Missing column fastq_1!' )
}else if ( ! 'fastq_2' %in% colnames(ss)){
    ss$fastq_2 <- ''
}
ss$fastq_2 <- add_col(ss, 'fastq_2', "")
ss$single_end <- ifelse(ss$fastq_2 == "", 'true', 'false')


## add metadata and collapse by id if necessary ####
ssv <- add_metadata(ss, meta_csv)
ssv %>% 
    write.table('samplesheet.valid.csv', sep = ',', quote = F, row.names = F)

## write fastq paths ####
ss %>% 
    filter(id %in% ssv$id) %>% 
    write.table('fq.csv', sep = ',', quote = F, row.names = F)

