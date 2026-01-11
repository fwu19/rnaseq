#!/usr/bin/env Rscript

# Author: @fwu19

options(stringsAsFactors = F)
options(scipen = 99)
library(dplyr)

## functions ####
## add default value if a column is missing
fill_column <- function(df, colv, default.value, na.value = NULL, missing.value = NULL){
    if (!colv %in% colnames(df)){
        df[,colv] <- default.value
    }
    
    if (!is.null(na.value)){
        df[,colv] <- ifelse(is.na(df[,colv]), na.value, df[,colv])
    }
    
    if (!is.null(missing.value)){
        df[,colv] <- ifelse(df[,colv] == "", missing.value, df[,colv])
    }
    
    return(df)
}

add_metadata <- function(ss, meta_csv){
    ## update with metadata if provided
    # metadata contains a required columns id and optional columns: sample_group, sample_replicate, target, control, call_peak, call_rep_peak, call_con_peak
    
    if (file_test('-f', meta_csv) & grepl('.csv$', meta_csv) & !grepl('dummy', meta_csv)){
        meta <- read.csv(meta_csv) %>% 
            mutate(
                id = gsub(' +|&', '-', id)
            )
        
        ss <- ss %>% 
            inner_join(
                meta, by = 'id', suffix = c(".x", "")
            ) %>% 
            dplyr::select(!ends_with(".x")) 
    }
    
    ss %>% 
        dplyr::relocate(id, fastq_1, fastq_2, sample_group) %>% 
        mutate(
            sample_group = gsub(' +|&', '-', sample_group)
        )
        
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
}
ss <- ss %>% 
    fill_column('sample_group', ss$id, ss$id, ss$id) %>% 
    fill_column('fastq_2', "", "", "") %>% 
    mutate(
        id = gsub(' +|&', '-', id),
        sample_group = gsub(' +|&', '-', sample_group),
        single_end = ifelse(fastq_2 == "", 'true', 'false')
    )


## add metadata and collapse by id if necessary ####
ss <- add_metadata(ss, meta_csv)
if(nrow(ss) == 0){
  stop(paste0('No input in the sample sheet!\nCheck ', in_csv, '\nIf ', meta_csv, ' is supplied, Column id should share at least one value with the same column in ', meta_csv, '.\n'))
}

## collapse sample sheet
if (sum(duplicated(ss$id)) > 0){
    ssu <- ss %>% 
        group_by(id) %>% 
        mutate(
            fastq_1 = paste(unique(basename(fastq_1)), collapse = ';'),
            fastq_2 = paste(unique(basename(fastq_2)), collapse = ';')
        ) %>% 
        unique.data.frame()
}else{
    ssu <- ss
}
ssu %>% 
    write.table('samplesheet.collapsed.csv', sep = ',', quote = F, row.names = F)


## write fastq paths
ss %>% 
    write.table('samplesheet.valid.csv', sep = ',', quote = F, row.names = F)

