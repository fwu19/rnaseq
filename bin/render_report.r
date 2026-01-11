#!/usr/bin/env Rscript

options(stringsAsFactor=F)
options(scipen = 99)

args <- as.vector(commandArgs(T))
for (arg in strsplit(args, split = '=')){
    assign(arg[1], arg[2])
}
if (sum(c('workflow','fdr','fdr2','fc','fc2') %in% ls()) < 5){
    stop(paste('The following arguments are missing!\n', setdiff(c('workflow','fdr','fdr2','fc','fc2'), ls()), '\n'))    
}

## copy Rmd files ####
if (workflow == 'regular'){
    src <- "report/rnaseq_regular.Rmd"
    dst <- "00_RNAseq_analysis_report.Rmd"
    
}else if (workflow == 'exome'){
    src <- "report/rnaseq_exome.Rmd"
    dst <- "00_RNAexome_analysis_report.Rmd"
    
}else if (workflow == 'pdx'){
    src <- "report/rnaseq_pdx.Rmd"
    dst <- "00_PDX_RNAseq_analysis_report.Rmd"
}else{
    cat(workflow, ' is not recognized!')
}

file.copy(src, dst)
rmarkdown::render(
    dst, 
    params = list(
        fdr = as.numeric(fdr), fc = as.numeric(fc), fdr2 = as.numeric(fdr2), fc2 = as.numeric(fc2)
    )
)
