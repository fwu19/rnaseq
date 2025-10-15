#!/usr/bin/env Rscript

options(stringsAsFactors = F)
options(scipen = 99)
library(dplyr)
library(ggplot2)
library(patchwork)

## functions ####
generate_count_matrix_star <- function(gene.txt, length.col, ss, count.files, count.col){

    ann <- read.delim(gene.txt)
    ann$gene_length <- ann[, length.col]
    
    cts <- do.call(cbind, lapply(
            count.files, 
            function(fname){
                df <- read.delim(fname, header = F, comment.char = '#', skip = 4)
                if(sum(!ann$gene_id %in% df[,1]) != 0){
                    warning(paste("Some genes in", gene.txt, "don't have counts!"))
                }
                v <- df[,count.col][match(ann$gene_id, df[,1])]
                v[is.na(v)] <- 0
                v
            }
        ))
    
    colnames(cts) <- gsub('.ReadsPerGene.*', '', basename(count.files))
    if (!setequal(colnames(cts), ss$id)){
        warning(paste("Some samples don't have counts!"))
    }
    
    return(list(ann = ann, cts = cts)) 
    
}

generate_count_matrix_featurecounts <- function(gene.txt, length.col, ss, count.files){
    
    ann <- read.delim(gene.txt)
    ann$gene_length <- ann[, length.col]
    
    cts <- do.call(cbind, lapply(
            count.files, 
            function(fname){
                df <- read.delim(fname, header = T, comment.char = '#')
                if(sum(!ann$gene_id %in% df$Geneid) != 0){
                    warning(paste("Some genes in", gene.txt, "don't have counts!"))
                }
                v <- df[,7][match(ann$gene_id, df$Geneid)]
                v[is.na(v)] <- 0
                v
            }
        ))
    
    colnames(cts) <- gsub('.exonic.*', '', basename(count.files))
    if (!setequal(colnames(cts), ss$id)){
        warning(paste("Some samples don't have counts!"))
    }
    return(list(ann = ann, cts = cts)) 
}

count2dgelist <- function(
    counts.tsv=NULL, feature.cols=1:8, pattern2remove="^X|.bam$", counts=NULL, features = NULL, samples = NULL, group.col = 'sample_group'
    ){
    options(stringsAsFactors = F)
    require(edgeR)
    
    if(is.null(counts)){
        counts <- read.delim(counts.tsv)
    }
    if(!is.null(pattern2remove)){
        colnames(counts) <- gsub(pattern2remove, '', colnames(counts))
    }
    if(is.null(features)){
        features <- counts[feature.cols]
        counts <- counts[-(feature.cols)]
    }
        
        
    y0 <- DGEList(counts=counts, genes=features, remove.zeros = T, samples = samples) 
    if(!is.null(samples) & group.col %in% colnames(samples)){
        y0$samples$group <- samples[,group.col]
    }
    y0 <- calcNormFactors(y0)
    
    saveRDS(y0, 'y0.rds')
    
    return(y0)
}


## read arguments ####
args <- as.vector(commandArgs(T)) 
lst <- strsplit(args, split = '=')
for (x in lst){
    assign(x[1],x[2])
} # read arguments: ss, count.dir, gene.txt, length.col

ss <- read.csv(input) %>% 
    mutate(id = factor(id)) %>% 
    relocate(fastq_1, fastq_2, .after = last_col()) %>% 
    unique.data.frame() # sample sheet

count.col <- as.integer(strand) + 2 # for strand=0,1,2


## generate count matrix ####
count.files <- list.files(count.dir, full.names = T)
if (grepl('ReadsPerGene.out.tab', count.files[1])){
    lst <- generate_count_matrix_star(gene.txt,length.col,ss,count.files,count.col)
}else if (grepl('exonic.*.txt', count.files[1])){
    lst <- generate_count_matrix_featurecounts(gene.txt,length.col,ss,count.files)
}else {
    stop (paste("Cannot parse files in", count.dir, "!"))
}
## write out raw counts ####
bind_cols(lst$ann, lst$cts) %>% 
    write.table('all_samples.gene_raw_counts.txt', sep = '\t', quote = F, row.names = F)


## create DGElist ####
y0 <- count2dgelist(
    counts = lst$cts, 
    features = lst$ann, 
    samples = ss %>% 
        arrange(factor(id, levels = colnames(lst$cts))),
    group = 'sample_group'
)
rpkm <- rpkm(y0, gene.length = "gene_length", normalized.lib.sizes = T, log = F)
cbind(y0$genes, rpkm) %>% 
    write.table('all_samples.gene_FPKM.txt', sep = '\t', quote = F, row.names = F)




