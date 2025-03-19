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
    cts <- cbind(data.frame(transcript_id = rownames(scaled.counts), catch$annotation) %>% 
            left_join(ann, by = 'transcript_id'), 
        scaled.counts
    )
    return(cts) 
    
}

count2dgelist <- function(counts.tsv=NULL, return.counts = T, pattern2remove="^X|.bam$", counts=NULL, out.dir=NULL, feature.cols=1:8, samples = NULL, group.col = 'sample_group'){
    options(stringsAsFactors = F)
    require(edgeR)
    
    if(is.null(counts)){
        counts <- read.delim(counts.tsv)
    }
    if(!is.null(pattern2remove)){
        colnames(counts) <- gsub(pattern2remove, '', colnames(counts))
    }
    
    y0 <- DGEList(counts=counts[-(feature.cols)], genes=counts[feature.cols], remove.zeros = T, samples = samples) 
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

## PCA
plot_pca <- function(y, out.prefix, var.genes = NULL, color = NULL, plot.title = '', sample.label = T, feature.length = 'gene_length', scree.plot = F){
    options(stringsAsFactors = F)
    require(ggrepel)
    require(edgeR)
    
    ## Create outdir if needed
    out.dir <- dirname(out.prefix)
    if(!dir.exists(out.dir)){dir.create(out.dir,recursive = T)}
    
    if('chrom' %in% colnames(y$genes)){
        y <- y[!y$genes$chrom %in% c('chrX', 'chrY', 'chrM'), ]
    }
    
    log2rpkm <- rpkm(y, gene.length = feature.length, normalized.lib.sizes = T, log = T)
    
    if (!is.null(var.genes)){
        keep <- rank(-apply(log2rpkm, 1, var)) < var.genes
        log2rpkm <- log2rpkm[keep,]
    }
    
    ## Run PCA
    pca <- prcomp(t(log2rpkm), center = T, scale = T)
    
    ## Scree plot
    pca.variance.prop <- (pca$sdev^2)/sum(pca$sdev^2)*100
    
    if (scree.plot){
        pdf(paste(out.dir,'pca.scree.plot.pdf',sep = '/'), width = 3, height = 3)
        barplot(
            pca.variance.prop[1:10],
            cex.names = 1,
            xlab = 'Principal component (PC), 1-10',
            ylab = 'Proportion of variance (%)',
            main = 'Scree plot',
            ylim = c(0,100)
        )
        
        points(
            cumsum(pca.variance.prop)[1:10], col = 'red', type = 'l'
        )
        dev.off()
    }
    
    ## PC1 vs PC2
    df <- cbind(pca$x[,1:2],data.frame(label=rownames(pca$x)))
    if(is.null(color)){
        df$color <- y$samples$group
    }else{
        df$color <- color
    }
    
    if(length(unique(df$color))>1){
        p <- ggplot(df,aes(x=PC1,y=PC2,label=label,color=color))
    }else{
        p <- ggplot(df,aes(x=PC1,y=PC2,label=label))
    }
    p <- p +
        geom_point(shape = 1)
    
    if(sample.label){
        p <- p +
            geom_text_repel(size = 2.4, color = 'black', position = 'jitter',max.overlaps = 80)
    }
    
    p <- p +
        labs(
            title = plot.title,
            color = '', 
            x = paste0('PC1 (',round(pca.variance.prop[1],1),'%)'),
            y = paste0('PC2 (',round(pca.variance.prop[2],1),'%)')
        )+
        theme_bw()+
        theme(
            legend.position = 'top'
        )
    ggsave(paste(out.prefix,'PCA.pdf',sep='.'),width = 5,height = 5)
    
    ## return data
    return(p)
    
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
cts <- generate_count_matrix_salmon(gene.txt, count.dirs)
cts %>% 
    write.table('all_samples.transcript_raw_counts.txt', sep = '\t', quote = F, row.names = F)

## create DGElist ####
y0 <- count2dgelist(
    counts = cts, 
    out.dir = NULL, 
    feature.cols = 1:11, 
    samples = ss %>% 
        arrange(factor(id, levels = colnames(cts)[12:ncol(cts)])),
    group.col = 'sample_group'
)

## plot PCA of all samples ####
p <- plot_pca(y0, out.prefix = 'all_samples', var.genes = 500, color = y0$samples$group, sample.label = T, feature.length = length.col)
saveRDS(p, 'pca.rds')


