#!/usr/bin/env Rscript

options(stringsAsFactors = F)
options(scipen = 99)
library(dplyr)
library(ggplot2)
library(patchwork)

## functions ####
generate_count_matrix <- function(gene.txt, length.col, ss, count.dir, count.col){

    ann <- read.delim(gene.txt) 
    
    cts <- cbind(
        cbind(ann[1:7], data.frame(gene_length = ann[, length.col])),
        lapply(
            file.path(count.dir, ss$sample_id, 'ReadsPerGene.out.tab'), 
            function(fname){
                df <- read.delim(fname, header = F, comment.char = '#', skip = 4)
                stopifnot(sum(!ann$gene_id %in% df[,1]) == 0)
                df[,count.col][match(ann$gene_id, df[,1])]
        }
    ))
    colnames(cts)[9:ncol(cts)] <- ss$sample_id
    
    return(cts)
}

count2dgelist <- function(counts.tsv=NULL, return.counts = T, pattern2remove="^X|.bam$", counts=NULL, out.dir=NULL, feature.cols=1:7, samples = NULL, group.col = 'sample_group'){
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
        if (is.null(count.tsv)){
            out.dir <- './'
        }else{
            out.dir <- dirname(counts.tsv)
        }
    }
    if(!dir.exists(out.dir)){dir.create(out.dir, recursive = T)}
    saveRDS(y0, paste(out.dir, 'y0.rds', sep = '/'))
    
    if(return.counts){
        write.table(cbind(y0$genes, y0$counts), paste(out.dir, 'all_samples.raw_counts.txt', sep = '/'), quote = F, row.names = F)
    }
    return(y0)
}

## for all
pca_log2rpkm <- function(y, out.dir, prefix, var.genes = NULL, color = NULL, plot.title = '', sample.label = T, feature.length = 'gene_length'){
    options(stringsAsFactors = F)
    require(ggrepel)
    
    if('chrom' %in% colnames(y$genes)){
        y <- y[!y$chrom %in% c('chrX', 'chrY', 'chrM'), ]
    }
    
    log2rpkm <- rpkm(y, gene.length = feature.length, normalized.lib.sizes = T, log = T)
    
    if (!is.null(var.genes)){
        keep <- rank(-apply(log2rpkm, 1, var)) < var.genes
        log2rpkm <- log2rpkm[keep,]
    }
    
    ## Run PCA
    pca <- prcomp(t(log2rpkm), center = T, scale = T)
    
    ## Create outdir if needed
    if(!dir.exists(out.dir)){dir.create(out.dir,recursive = T)}
    
    ## Scree plot
    pca.variance.prop <- (pca$sdev^2)/sum(pca$sdev^2)*100
    
    pdf(paste(out.dir,'pca.scree.plot.pdf',sep = '/'))
    barplot(
        pca.variance.prop[1:50],
        cex.names = 1,
        xlab = 'Principal component (PC), 1-50',
        ylab = 'Proportion of variance (%)',
        main = 'Scree plot',
        ylim = c(0,80)
    )
    
    points(
        cumsum(pca.variance.prop)[1:50], col = 'red', type = 'l'
    )
    dev.off()
    
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
    
    p +
        labs(
            title = plot.title,
            color = '', 
            x = paste0('PC1 (',round(pca.variance.prop[1],1),'%)'),
            y = paste0('PC2 (',round(pca.variance.prop[2],1),'%)')
        )+
        theme_bw()
    ggsave(paste(out.dir,paste(prefix,'PCA.pdf',sep='.'),sep = '/'),width = 6,height = 5)
    
    ## return data
    return(pca)
    
}

## for all
run_da <- function(
    y0, out.prefix, 
    control.group, test.group, group=NULL, 
    fdr=0.01, fc=2, fdr2=NULL, fc2=NULL, 
    report.cpm=F, report.rpkm=T, 
    TMM=T, method = 'QL', 
    rename.feature = NULL, feature.length = 'length', 
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
        
        pdf(paste(out.dir,'qldisp.pdf',sep = '/'),width = 4, height = 4); plotQLDisp(fit); dev.off()
        
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
    saveRDS(list(y=y, design=design, fit=fit, test=test), paste(out.dir,'da.rds',sep = '/'))
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
    
    return(list(summary = df_sum, y = y, df = df))
    
}


## PCA
plot_pca <- function(y, out.prefix, var.genes = NULL, color = NULL, plot.title = '', sample.label = T, feature.length = 'gene_length'){
    options(stringsAsFactors = F)
    require(ggrepel)
    require(edgeR)
    
    ## Create outdir if needed
    out.dir <- dirname(out.prefix)
    if(!dir.exists(out.dir)){dir.create(out.dir,recursive = T)}
    
    
    log2rpkm <- rpkm(y, gene.length = feature.length, normalized.lib.sizes = T, log = T)
    
    if (!is.null(var.genes)){
        keep <- rank(-apply(log2rpkm, 1, var)) < var.genes
        log2rpkm <- log2rpkm[keep,]
    }
    
    ## Run PCA
    pca <- prcomp(t(log2rpkm), center = T, scale = T)
    
    ## Scree plot
    pca.variance.prop <- (pca$sdev^2)/sum(pca$sdev^2)*100
    
    pdf(paste(out.dir,'pca.scree.plot.pdf',sep = '/'))
    barplot(
        pca.variance.prop[1:50],
        cex.names = 1,
        xlab = 'Principal component (PC), 1-50',
        ylab = 'Proportion of variance (%)',
        main = 'Scree plot',
        ylim = c(0,80)
    )
    
    points(
        cumsum(pca.variance.prop)[1:50], col = 'red', type = 'l'
    )
    dev.off()
    
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

## Plot MD
plot_MD <- function(df, out.prefix, plot.title = ""){
    require(ggplot2)
    
    ## Create outdir if needed
    out.dir <- dirname(out.prefix)
    if(!dir.exists(out.dir)){dir.create(out.dir,recursive = T)}
    
    p <- ggplot(df, aes(x = logCPM, y = logFC, color=factor(is.sig)))+
        geom_hline(yintercept = 0)+
        geom_point(size = 0.2)+
        scale_color_manual(
            '',
            values = c("-1" = "blue", "0" = "gray", "1" = "red"),
            breaks = c(-1,0,1), labels = c('Down', 'No Sig.', 'Up')
        )+
        labs( 
            x = "Average log CPM", 
            y = "log-fold-of-change",
            title = plot.title
        )+
        theme_bw()+
        theme(
            text = element_text(size = 8),
            legend.position = 'top'
        )
    
    ggsave(paste(out.prefix,'MD.pdf',sep = '.'),width = 4,height = 5)
    return(p)
} 

## Plot volcano 
plot_volcano <- function(df, out.prefix, plot.title = ""){
    require(ggplot2)
    
    ## Create outdir if needed
    out.dir <- dirname(out.prefix)
    if(!dir.exists(out.dir)){dir.create(out.dir,recursive = T)}
    
    p <- ggplot(df,aes(x=logFC,y=-log10(FDR),color=factor(is.sig)))+
        geom_point(size = 0.2)+ 
        scale_color_manual(
            values = c('-1'='blue','0'='gray','1'='red'),
            breaks = c('-1','0','1'),
            labels = c('Down','No Sig.','Up'),
            drop = T
        )+
        labs(
            x='log2(fold change)',
            y='-log10(FDR)',
            color='',
            title = plot.title
        )+
        theme_bw()+
        theme(
            text = element_text(size = 8),
            legend.position = 'top'
        )
    
    ggsave(paste(out.prefix,'Volcano.pdf',sep = '.'),width = 4,height = 5)
    return(p)
}

## recompute is.sig2
recal_sig <- function(txt, col.sig, fdr, fc){
    de <- read.delim(txt)
    de[,col.sig] <- sign(de$logFC) * (abs(de$logFC) > log2(fc)) * (de$FDR < fdr)
    write.table(de, txt, sep = '\t', quote = F, row.names = F)
    return(
        data.frame(
            source.file = basename(txt),
            modify.col = col.sig,
            features.test = nrow(de),
            features.up = sum(de[,col.sig] %in% 1),
            features.down = sum(de[,col.sig] %in% -1),
            FDR.cutoff = fdr,
            FC.cutoff = fc
        )
    )
}

## wrapper
wrap_one_cmp <- function(icmp, ss, fdr = 0.05, fc = 1.5, fdr2 = 0.01, fc2 = 2){

    out.prefix <- icmp$out.prefix[1]
    control.group <- icmp$control.group[[1]]
    test.group <- icmp$test.group[[1]]
    plot.title <- icmp$plot.title[1]
    
    ## run DGE ####
    lst <- run_da(
        y0, 
        out.prefix = out.prefix, 
        control.group = control.group, 
        test.group = test.group, 
        group = y0$samples$sample_group, 
        feature.length = 'gene_length',
        fdr = fdr, fc = fc, fdr2 = fdr2, fc2 = fc2
    )
    
    y <- lst$y
    df <- lst$df
    lst$plots <- list(
        PCA = plot_pca(
            y, 
            out.prefix, 
            color = y$samples$sample_group, 
            sample.label = T, 
            plot.title = "", 
            var.genes = 500, 
            feature.length = "gene_length"),
        MD = plot_MD(
            df, out.prefix = out.prefix, 
            plot.title = plot.title
        ),
        volcano = plot_volcano(
            df, out.prefix = out.prefix, 
            plot.title = plot.title
        )
    )
    
    ##
    return(lst)


}


## read arguments ####
out.dir <- 'differential_expression' # path/to/output/directory
if (!dir.exists(out.dir)){dir.create(out.dir, recursive = T)}

args <- as.vector(commandArgs(T)) 
params <- rjson::fromJSON(file = args[1])

ss <- read.csv(params$input) %>% 
    mutate(sample_group = gsub('-','_', sample_group)) %>% 
    dplyr::select(-c(fastq_1, fastq_2)) %>% 
    unique.data.frame() # sample sheet

cmp <- readRDS(params$comparison)

count.dir <- file.path(params$outdir, 'STAR') # path/to/STAR_alignment
gene.txt <- params$geneTxt
length.col <- params$lengthCol # full_exome,TruSeq_RNA_exome
count.col <- params$strand + 2 # for strand=0,1,2
fdr <- params$fdr
fc <- params$fc
fdr2 <- params$fdr2
fc2 <- params$fc2

## generate count matrix ####
cts <- generate_count_matrix(gene.txt,length.col,ss,count.dir,count.col)

## create DGElist ####
y0 <- count2dgelist(
    counts = cts, 
    out.dir = out.dir, 
    feature.cols = 1:8, 
    samples = ss %>% 
        arrange(factor(sample_id, levels = colnames(cts)[9:ncol(cts)])),
    group = 'sample_group'
)
y0$samples$group <- y0$samples$sample_group

## plot PCA of all samples ####
plot_pca(y0, out.prefix = paste(out.dir, 'all_samples', sep = '/'), var.genes = 500, color = y0$samples$group, sample.label = T, feature.length = "gene_length")

## detect differential expression ####
de.list <- list()
for (i in 1:nrow(cmp)){
    de.list[[i]] <- wrap_one_cmp(cmp[i,], ss, fdr, fc, fdr2, fc2)
}
names(de.list) <- basename(cmp$out.prefix)
saveRDS(de.list, paste(out.dir, 'de.rds', sep = '/'))

de_sum <- bind_rows(lapply(de.list, function(de){de$summary})) %>% 
    dplyr::relocate(output.folder, control.group, test.group, features.tested:FC.cutoff, control.samples:test.samples)
de_sum %>% 
    write.table(paste(out.dir, 'DE_summary.txt', sep = '/'), sep = '\t', quote = F, row.names = F)