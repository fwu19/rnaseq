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

## output genes.txt ####
exons <- gtf[gtf$type == 'exon']
gene_lens <- bind_rows(
    lapply(
    split(exons, exons$gene_id), 
    function(x){
        gr <- GenomicRanges::reduce(x)
        data.frame(
            gene_id = x$gene_id[1],
            gene_length = sum(width(gr))
        )
    })
)

genes <- as.data.frame(gtf[gtf$type == 'gene']) %>% 
    mutate(star = start - 1) %>% 
    dplyr::rename("chrom" = "seqnames") %>% 
    left_join(
        gene_lens, by = 'gene_id'
    )
genes <- genes[intersect(colnames(genes), c('chrom', 'start', 'end', 'strand', 'gene_id', 'gene_name', 'gene_biotype', 'gene_type', 'gene_length'))]
write.table(genes, 'genes.txt', sep = '\t', quote = F, row.names = F)

