#!/usr/bin/env Rscript

## process GTF files

options(stringsAsFactors = F)
library(rtracklayer)
library(GenomicRanges)
library(dplyr)

## read arguments
args <- as.vector(commandArgs(T)) 
gtf_file <- args[1]

## process gtf ####
stopifnot(file.exists(gtf_file))
gtf <- rtracklayer::import(gtf_file)

## optional
# gtf <- keepStandardChromosomes(gtf, pruning.mode = 'coarse')


## output transcript.txt
tx <- as.data.frame(gtf[gtf$type == 'transcript']) %>% 
    mutate(star = start - 1) %>% 
    dplyr::rename("chrom" = "seqnames") 
tx <- tx[intersect(colnames(tx), c('chrom', 'start', 'end', 'strand', 'gene_id', 'gene_name', 'gene_biotype', 'gene_type', 'transcript_id', 'transcript_type'))]
write.table(tx, gsub('gtf$','transcripts.txt', basename(gtf_file)), sep = '\t', quote = F, row.names = F)



