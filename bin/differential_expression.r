#!/usr/bin/env Rscript

options(stringsAsFactors = F)
options(scipen = 99)
library(dplyr)
library(ggplot2)
library(patchwork)

## functions ####

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
    
    j1 <- sum(y0$samples$group %in% control.group)
    j2 <- sum(y0$samples$group %in% test.group)
    if (j1 == 0 | j2 == 0){
        return(list(error = data.frame(
            control.group = paste(control.group, collapse = '+'),
            test.group = paste(test.group, collapse = '+'),
            control.samples = j1,
            test.samples = j2
        )))
    }else if (j1 == 1 & j2 == 1){
        return(list(error = data.frame(
            control.group = paste(control.group, collapse = '+'),
            test.group = paste(test.group, collapse = '+'),
            control.samples = j1,
            test.samples = j2
        )))        
    }
    
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
    
    ## revert y$group ####
    y$samples$group <- ifelse(y$samples$group %in% 'control', control.group, test.group)
    
    ## write out results ####
    
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
        keep <- rank(-apply(log2rpkm, 1, var)) <= var.genes
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

## wrapper
wrap_one_cmp <- function(y0, icmp, ss, fdr = 0.05, fc = 1.5, fdr2 = 0.01, fc2 = 2, genes2keep = NULL, gene_types = NULL){

    out.prefix <- icmp$out.prefix[1]
    control.group <- icmp$control.group[[1]]
    test.group <- icmp$test.group[[1]]
    plot.title <- icmp$plot.title[1]
    if (icmp$sample.group[1] %in% colnames(y0$samples)){
        group <- y0$samples[,icmp$sample.group[1]]
    }else{
        group <- NULL
    }
    
    ## filter gene if needed ####
    if (!is.null(genes2keep)){
        y0 <- y0[y0$genes$gene_id %in% genes2keep$gene_id, ]
    }
    if (!is.null(gene_types)){
        y0 <- y0[y0$genes$gene_type %in% gene_types, ]
    }
    
    ## run DGE ####
    lst <- run_da(
        y0, 
        out.prefix = file.path(out.prefix, out.prefix),
        control.group = control.group, 
        test.group = test.group, 
        group = group, 
        feature.length = 'gene_length',
        fdr = fdr, fc = fc, fdr2 = fdr2, fc2 = fc2
    )
    
    if ('error' %in% names(lst)){
        write.table(lst$error, file.path(out.prefix, paste0(out.prefix, '.error.txt')), sep = '\t', quote = F, row.names = F)
        return(lst)
        }
    
    y <- lst$y
    df <- lst$df
    lst$plots <- list(
        PCA = plot_pca(
            y, 
            file.path(out.prefix, out.prefix), 
            color = y$samples$group, 
            sample.label = T, 
            plot.title = "", 
            var.genes = 500, 
            feature.length = "gene_length"),
        MD = plot_MD(
            df, out.prefix = file.path(out.prefix, out.prefix),
            plot.title = plot.title
        ),
        volcano = plot_volcano(
            df, out.prefix = file.path(out.prefix, out.prefix),
            plot.title = plot.title
        )
    )
    
    ##
    return(lst)


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
} # read arguments: ss, comparison, rds, fdr, fc, fdr2, fc2

ss <- read.csv(input) %>% 
    relocate(fastq_1, fastq_2, .after = last_col()) %>% 
    unique.data.frame() # sample sheet

if(grepl('dummy', comparison)){
    cat(comparison, "is a dummy file! Provide --comparison path/to/comparison_file (a comparison table in csv, txt, tsv or rds format)!")
    quit()
}else if(grepl('.csv$', comparison)){
    cmp <- read.csv(comparison)
}else if (grepl('.rds$', comparison)){
    cmp <- readRDS(comparison)
}else {
    stop(paste(comparison, 'should be either csv or rds file!'))
}

y0 <- readRDS(rds)
fdr <- as.numeric(fdr)
fc <- as.numeric(fc)
fdr2 <- as.numeric(fdr2)
fc2 <- as.numeric(fc2)
if (file.exists(gene_txt)){
    genes2keep <- read.delim(gene_txt)
}else{
    genes2keep <- NULL
}

gene_types <- unlist(strsplit(gene_type, split = ','))
if(setequal(gene_types, 'all')){gene_types <- NULL}

## detect differential expression ####
cmp <- cmp %>% 
    add_colv('out.prefix', paste(cmp$test.group, cmp$control.group, sep = '_vs_')) %>% 
    add_colv('plot.title', paste(cmp$test.group, cmp$control.group, sep = ' vs ')) %>% 
    add_colv('sample.group', '')

de.list <- list()
for (i in 1:nrow(cmp)){
    de.list[[i]] <- wrap_one_cmp(y0, cmp[i,], ss, fdr, fc, fdr2, fc2, genes2keep, gene_types)
}
names(de.list) <- basename(cmp$out.prefix)
saveRDS(de.list, 'differential_genes.rds')
