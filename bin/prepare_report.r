#!/usr/bin/env Rscript

options(stringsAsFactor=F)
options(scipen = 99)

library(dplyr)

args <- as.vector(commandArgs(T))
workflow <- args[1]
report_dir <- args[2]

dat <- list()
dat$ss <- read.csv('sample_sheet.csv')

## QC metrics ####
if (dir.exists('multiqc_data')){
    fname <- 'multiqc_data/multiqc_general_stats.txt'
    if(file.exists(fname)){
        dat$stat <- read.delim(fname)
    }

    fname <- 'multiqc_data/multiqc_star.txt'
    if(file.exists(fname)){
        dat$star <- read.delim(fname)
    }
    
    fname <- 'multiqc_data/multiqc_rna_seqc.txt'
    if(file.exists(fname)){
        dat$rnaseqc <- read.delim(fname)
    }
    
}

## collect_hs_metrics ####
if (dir.exists('hs_metrics/')){
    flist <- list.files('hs_metrics/', full.names = T)
    dat$gatk <- data.table::rbindlist(lapply(
        flist, function(fname){
            read.delim(fname, comment.char = '#', nrows = 1) %>% 
                mutate(Sample = gsub('.hs_metrics.txt', '', basename(fname)))
        }
    ), fill = T, use.names = T)
    
}


## save data.rds ####
saveRDS(dat, 'data.rds')

## copy Rmd files ####
if (workflow == 'regular'){
    file.copy(file.path(report_dir, "rnaseq_regular.Rmd"), "00_RNAseq_analysis_report.Rmd")
}else if (workflow == 'exome'){
    file.copy(file.path(report_dir, "rnaseq_exome.Rmd"), "00_RNAexome_analysis_report.Rmd")
}else if (workflow == 'pdx'){
    file.copy(file.path(report_dir, "rnaseq_pdx.Rmd"), "00_PDX_RNAseq_analysis_report.Rmd")
}else{
    cat(workflow, ' is not recognized!')
}
