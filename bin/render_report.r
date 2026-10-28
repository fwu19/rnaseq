#!/usr/bin/env Rscript

options(stringsAsFactor=F)
options(scipen = 99)

args <- as.vector(commandArgs(T))

lst <- list()
for (arg in strsplit(args, split = '=')){
    key <- trimws(arg[1])
    value <- trimws(arg[2])
    
    if (key %in% c('report_rmd', 'report_html')){
        assign(key, value)
    }else{
        lst[[key]] <- value 
    }
    
}


for (i in c('fdr', 'fdr2', 'fc', 'fc2', 'min_samples')){
    if(i %in% names(lst)){lst[[i]] <- as.numeric(lst[[i]])}
}

## copy Rmd files ####
src <- report_rmd
dst <- gsub('.html$', '.Rmd', report_html)
file.copy(src, dst)
rmarkdown::render(
    dst, 
    params = lst
)
