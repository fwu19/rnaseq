#!/usr/bin/env Rscript

options(stringsAsFactors = F)
options(scipen = 99)
library(dplyr)
library(ggplot2)
library(patchwork)

## functions ####
generate_count_matrix_star <- function(gene.txt, length.col, ss, count.files, count.col){

    ann <- read.delim(gene.txt) 
    
    
    cts <- cbind(
        ann[1:7], 
        data.frame(gene_length = ann[, length.col]),
        do.call(cbind, lapply(
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
    )
    colnames(cts)[9:ncol(cts)] <- gsub('.ReadsPerGene.*', '', basename(count.files))
    if (!setequal(colnames(cts)[9:ncol(cts)], ss$id)){
        warning(paste("Some samples don't have counts!"))
    }
    
    return(cts) 
    
}

generate_count_matrix_featurecounts <- function(gene.txt, length.col, ss, count.files){
    
    ann <- read.delim(gene.txt) 
    
    cts <- cbind(
        ann[1:7], 
        data.frame(gene_length = ann[, length.col]),
        do.call(cbind, lapply(
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
    )
    colnames(cts)[9:ncol(cts)] <- gsub('.exonic.*', '', basename(count.files))
    if (!setequal(colnames(cts)[9:ncol(cts)], ss$id)){
        warning(paste("Some samples don't have counts!"))
    }
    return(cts)
}

count2dgelist <- function(counts.tsv=NULL, pattern2remove="^X|.bam$", counts=NULL, feature.cols=1:8, samples = NULL, group.col = 'sample_group'){
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
    
    saveRDS(y0, 'y0.rds')
    
    return(y0)
}

## for all
run_da <- function(
    y0, out.prefix, 
    control.group, test.group, group=NULL, 
    fdr=0.01, fc=2, fdr2=NULL, fc2=NULL, 
    report.cpm=F, report.rpkm=T, 
    TMM=T, method = 'QL', 
    rename.feature = NULL, feature.length = 'gene_length', 
    design.object = ~0+group,
    target = NULL
){
    require(edgeR)
    
    ## create output directory ####
    out.dir <- dirname(out.prefix)
    if(!dir.exists(out.dir)){dir.create(out.dir, recursive = T)}
    
    ## retrieve and process data ####
    if(!is.null(group)){y0$samples$group <- group}
    j <- y0$samples$group %in% c(control.group, test.group)
    
    y <- y0[,j]
    y$samples$group <- ifelse(y$samples$group %in% control.group, 'control', 'test')
    
    if(TMM){
        keep <- filterByExpr(y, group = y$samples$group, min.count=10, min.total.count = 15)
        
        y <- y[keep,,keep.lib.sizes=F]
        y<-calcNormFactors(y)
    }
    
    
    ## Create design matrix ####
    design <- model.matrix(design.object,data = y$samples)
    colnames(design) <- gsub('^group','',colnames(design))
    
    ## Make contrasts ####
    contrasts <- makeContrasts(
        cmp = test - control,
        levels = design
    )
    
    ## Run DP test ####
    y <- estimateDisp(y, design, robust = T)
    
    pdf(paste(out.dir,'bcv.pdf',sep = '/'),width = 4, height = 4); plotBCV(y); dev.off()
    
    if (method == 'QL'){
        fit<-glmQLFit(y, design=design, dispersion = y$trended.dispersion, robust = T)
        
        # pdf(paste(out.dir,'qldisp.pdf',sep = '/'),width = 4, height = 4); plotQLDisp(fit); dev.off()
        
        test <- glmQLFTest(fit, contrast = contrasts)
    }else{
        fit <- glmFit(y, design=design, dispersion = y$trended.dispersion, robust = T)
        test <- glmLRT(fit, contrast = contrasts)
    }
    
    ## prepare for plots
    is.sig <- decideTests(object = test,adjust.method = 'BH',p.value = fdr,lfc = log2(fc))
    df <- test$table
    df$FDR <- p.adjust(df$PValue, method = 'BH')
    df$is.sig <- as.vector(is.sig)
    
    ## convert y$group back ####
    y$samples$group <- ifelse(y$samples$group %in% 'control', control.group, test.group)
    
    ## write out results ####
    # saveRDS(list(y=y, design=design, test=test), paste(out.dir,'da.rds',sep = '/'))
    if(!is.null(rename.feature)){colnames(test$genes)[1] <- rename.feature}
    df <- cbind(test$genes,df[c('logFC','logCPM','PValue','FDR','is.sig')])
    if(!is.null(fdr2) & !is.null(fc2)){
        df$is.sig2 <- (df$FDR < fdr2) * sign(df$logFC) * (abs(df$logFC) > log2(fc2))
    }
    if(report.cpm){
        cpm <- cpm(y, normalized.lib.sizes = T, log = F)
        colnames(cpm) <- paste('CPM.TMMnormalized', colnames(cpm), sep = '.')
        df <- cbind(df, cpm)
    }  
    if(report.rpkm){
        rpkm <- rpkm(y, gene.length = feature.length, normalized.lib.sizes = T, log = F)
        colnames(rpkm) <- paste('FPKM.TMMnormalized', colnames(rpkm), sep = '.')
        df <- cbind(df, rpkm)
    }  
    write.table(df,paste(out.prefix,'txt',sep = '.'), sep = '\t',quote = F,row.names = F)
    
    ## return results ####
    df_sum <- data.frame(
        output.folder = gsub('.*differential_peaks/', '', out.dir), 
        control.group = paste(control.group, collapse = '+'),
        test.group = paste(test.group, collapse = '+'),
        control.samples = sum(y$samples$group %in% control.group),
        test.samples = sum(y$samples$group %in% test.group),
        features.tested = nrow(y),
        features.up = sum(is.sig %in% 1),
        features.down = sum(is.sig %in% -1),
        FDR.cutoff = fdr,
        FC.cutoff = fc
    )
    if(!is.null(target)){
        df_sum$target <- target  
    }
    
    if('is.sig2' %in% colnames(df)){
        df_sum <- cbind(
            df_sum,
            data.frame(
                features.up2 = sum(df$is.sig2 %in% 1),
                features.down2 = sum(df$is.sig2 %in% -1),
                FDR.cutoff2 = fdr2,
                FC.cutoff2 = fc2
            )
        )
    }
    
    return(list(summary = df_sum, y = y, df = df, design=design, test=test))
    
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
    df <- cbind(pca$x[,1:2],data.frame(label=factor(rownames(pca$x))))
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
} # read arguments: ss, count.dir, gene.txt, length.col

ss <- read.csv(input) %>% 
    mutate(id = factor(id)) %>% 
    relocate(fastq_1, fastq_2, .after = last_col()) %>% 
    unique.data.frame() # sample sheet

count.col <- as.integer(strand) + 2 # for strand=0,1,2


## generate count matrix ####
count.files <- list.files(count.dir, full.names = T)
if (grepl('ReadsPerGene.out.tab', count.files[1])){
    cts <- generate_count_matrix_star(gene.txt,length.col,ss,count.files,count.col)
}else if (grepl('exonic.*.txt', count.files[1])){
    cts <- generate_count_matrix_featurecounts(gene.txt,length.col,ss,count.files)
}else {
    stop (paste("Cannot parse files in", count.dir, "!"))
}


## create DGElist ####
y0 <- count2dgelist(
    counts = cts, 
    feature.cols = 1:8, 
    samples = ss %>% 
        arrange(factor(id, levels = colnames(cts)[9:ncol(cts)])),
    group = 'sample_group'
)

## write out raw counts ####
if(workflow %in% 'exome' & !length.col %in% 'full_exome'){
    ann <- read.delim(gene.txt) %>% 
        dplyr::select(gene_id, full_exome)
    cts <- cts %>% 
        left_join(
            ann, by = 'gene_id'
        ) %>% 
        relocate(full_exome, .after = 'gene_length')
    colnames(cts)[8:9] <- paste(c(length.col, 'full_exome'), 'length', sep = '_')
        
}
write.table(cts, paste('all_samples', 'raw_counts.txt', sep = '.'), sep = '\t', quote = F, row.names = F)

## plot PCA of all samples ####
p <- plot_pca(y0, out.prefix = 'all_samples', var.genes = 500, color = factor(y0$samples$group), sample.label = T, feature.length = "gene_length")
saveRDS(p, 'pca.rds')


