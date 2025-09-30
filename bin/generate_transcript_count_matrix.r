#!/usr/bin/env Rscript

options(stringsAsFactors = F)
options(scipen = 99)
library(dplyr)
library(ggplot2)
library(patchwork)

## functions ####
generate_count_matrix_salmon <- function(gene.txt, count.dirs){
    require(edgeR)
    
    ann <- read.delim(gene.txt) 
    
    catch <- catchSalmon(paths = count.dirs)
    scaled.counts <- catch$counts/catch$annotation$Overdispersion
    colnames(scaled.counts) <- basename(colnames(scaled.counts))
    lst <- list(
        features = data.frame(transcript_id = rownames(scaled.counts), catch$annotation) %>% 
            left_join(ann, by = 'transcript_id'), 
        counts = scaled.counts
    )
    return(lst) 
    
}

count2dgelist <- function(
    counts.tsv = NULL, return.counts = T, 
    pattern2remove="^X|.bam$", 
    counts = NULL, features = NULL,
    feature.cols=NULL, 
    samples = NULL, 
    group.col = 'sample_group',
    out.dir=NULL
){
    options(stringsAsFactors = F)
    require(edgeR)
    
    if(is.null(counts)){
        counts <- read.delim(counts.tsv)
    }
    if(!is.null(pattern2remove)){
        colnames(counts) <- gsub(pattern2remove, '', colnames(counts))
    }
    
    if (!is.null(feature.cols)){
        counts <- counts[-(feature.cols)]
        features <- cbind(features, counts[feature.cols])
    }
    
    y0 <- DGEList(
        counts = counts, 
        genes = features, 
        remove.zeros = T, 
        samples = samples
        ) 
    if(!is.null(samples) & group.col %in% colnames(samples)){
        y0$samples$group <- samples[,group.col]
    }
    y0 <- calcNormFactors(y0)
    
    if(is.null(out.dir)){
        if (is.null(counts.tsv)){
            out.dir <- './'
        }else{
            out.dir <- dirname(counts.tsv)
        }
    }
    if(!dir.exists(out.dir)){dir.create(out.dir, recursive = T)}
    saveRDS(y0, 'y0.rds')
    return(y0)
}

## compute normalized counts 
normalize_counts <- function(y, out.prefix, return = c('rpkm','cpm'), gene.length = "gene_length", log = F){
    require(edgeR)
    
    out.dir <- dirname(out.prefix)
    if (!dir.exists(out.dir)){
        dir.create(out.dir, recursive = T)
    }
    if(return[1] == 'cpm'){
        cpm <- cpm(y, normalized.lib.sizes = T, log = log)
        colnames(cpm) <- paste('CPM.TMMnormalized', colnames(cpm), sep = '.')
        df <- cbind(y$genes, cpm)
        out.suffix <- ifelse(log, 'log2CPM.txt', 'CPM.txt')
        write.table(df,paste(out.prefix, out.suffix, sep = '.'), sep = '\t',quote = F,row.names = F)
    }  
    if(return[1] == 'rpkm'){
        rpkm <- rpkm(y, gene.length = gene.length, normalized.lib.sizes = T, log = log)
        colnames(rpkm) <- paste('FPKM.TMMnormalized', colnames(rpkm), sep = '.')
        df <- cbind(y$genes, rpkm)
        out.suffix <- ifelse(log, 'log2FPKM.txt', 'FPKM.txt')
        write.table(df,paste(out.prefix, out.suffix, sep = '.'), sep = '\t',quote = F,row.names = F)
    }  
    
}

## add default value if a column is missing
add_colv <- function(df, colv, value){
    if (!colv %in% colnames(df)){df[,colv] <- value}
    return(df)
}

## read arguments ####
args <- as.vector(commandArgs(T)) 
lst <- strsplit(args, split = '=')
for (x in lst){
    assign(x[1],x[2])
} # read arguments: ss, comparison, count.dir, gene.txt, length.col

ss <- read.csv(input) %>% 
    relocate(fastq_1, fastq_2, .after = last_col()) %>% 
    unique.data.frame() # sample sheet

## generate count matrix ####
count.dirs <- list.dirs(count.dir, full.names = T, recursive = F)
lst <- generate_count_matrix_salmon(gene.txt, count.dirs)
lst  %>% 
    bind_cols() %>% 
    write.table('all_samples.transcript_raw_counts.txt', sep = '\t', quote = F, row.names = F)

## create DGElist ####
y0 <- count2dgelist(
    counts = lst$counts, features = lst$features,
    out.dir = NULL, 
    feature.cols = NULL, 
    samples = ss %>% 
        arrange(factor(id, levels = colnames(lst$counts))),
    group.col = 'sample_group'
)

