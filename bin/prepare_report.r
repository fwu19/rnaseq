#!/usr/bin/env Rscript

options(stringsAsFactor=F)
options(scipen = 99)

library(dplyr)

args <- as.vector(commandArgs(T))
workflow <- args[1]

dat <- list()
dat$ss <- read.csv('sample_sheet.csv')

## RNA-SeQC ####
if (dir.exists('multiqc_data')){
    if(file.exists('multiqc_data/multiqc_rna_seqc.txt')){
        dat$rnaseqc <- read.delim('multiqc_data/multiqc_rna_seqc.txt')
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
    file.copy("report/rnaseq_regular.Rmd", "00_RNAseq_analysis_report.Rmd")
}else if (workflow == 'exome'){
    file.copy("report/rnaseq_exome.Rmd", "00_RNAexome_analysis_report.Rmd")
}else if (workflow == 'pdx'){
    file.copy("report/rnaseq_pdx.Rmd", "00_PDX_RNAseq_analysis_report.Rmd")
}else{
    cat(workflow, ' is not recognized!')
}
