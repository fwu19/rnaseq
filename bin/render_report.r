#!/usr/bin/env Rscript

options(stringsAsFactor=F)
options(scipen = 99)

args <- as.vector(commandArgs(T))
# report_rmd=${report_rmd} workflow=${params.workflow} fdr=${params.fdr} fc=${params.fc} fdr2=${params.fdr2} fc2=${params.fc2}
for (arg in strsplit(args, split = '=')){
    assign(trimws(arg[1]), trimws(arg[2]))
}

if (sum(c('workflow','fdr','fdr2','fc','fc2') %in% ls()) < 5){
    stop(paste('The following arguments are missing!\n', setdiff(c('workflow','fdr','fdr2','fc','fc2'), ls()), '\n'))    
}

## copy Rmd files ####
src <- report_rmd
dst <- gsub('.html$', '.Rmd', report_html)
file.copy(src, dst)
rmarkdown::render(
    dst, 
    params = list(
        fdr = as.numeric(fdr), fc = as.numeric(fc), fdr2 = as.numeric(fdr2), fc2 = as.numeric(fc2)
    )
)
