#!/usr/bin/env Rscript

options(stringsAsFactors = F)
options(scipen = 99)
library(dplyr)
library(ggplot2)
library(patchwork)

## functions ####

## for all
run_da <- function(
    y0, out_prefix, 
    control_group, test_group, group=NULL, 
    exclude_samples = NULL, include_samples = NULL,
    fdr=0.01, fc=2, fdr2=NULL, fc2=NULL, 
    report_cpm=F, report_rpkm=T, 
    TMM=T, method = 'QL', 
    rename_feature = NULL, feature_length = 'gene_length', 
    design_object = ~0+group,
    target = NULL
){
    require(edgeR)
    
    ## create output directory ####
    out.dir <- dirname(out_prefix)
    if(!dir.exists(out.dir)){dir.create(out.dir, recursive = T)}
    
    ## retrieve and process data ####
    if(!is.null(group)){y0$samples$group <- group}
    
    j <- y0$samples$group %in% c(control_group, test_group)
    if (!is.null(exclude_samples)){
        j <- j & !rownames(y0$samples) %in% exclude_samples
    }
    
    if (!is.null(include_samples)){
        j <- j & rownames(y0$samples) %in% include_samples
    }
    
    if (sum(j)==0){
        return(list(error = data.frame(
            control_group = paste(control_group, collapse = ';'),
            test_group = paste(test_group, collapse = ';'),
            control_samples = 0,
            test_samples = 0
        )))   
    }
    
    y <- y0[,j]
    
    j1 <- sum(y$samples$group %in% control_group)
    j2 <- sum(y$samples$group %in% test_group)
    if (j1 == 0 | j2 == 0){
        return(list(error = data.frame(
            control_group = paste(control_group, collapse = ';'),
            test_group = paste(test_group, collapse = ';'),
            control_samples = j1,
            test_samples = j2
        )))
    }else if (j1 == 1 & j2 == 1){
        return(list(error = data.frame(
            control_group = paste(control_group, collapse = ';'),
            test_group = paste(test_group, collapse = ';'),
            control_samples = j1,
            test_samples = j2
        )))        
    }
    
    
    y$samples$group <- ifelse(y$samples$group %in% control_group, 'control', 'test')
    
    if(TMM){
        keep <- filterByExpr(y, group = y$samples$group, min.count=10, min.total.count = 15)
        
        y <- y[keep,,keep.lib.sizes=F]
        y<-calcNormFactors(y)
    }
    
    
    ## Create design matrix ####
    design <- model.matrix(design_object,data = y$samples)
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
    y$samples$group <- ifelse(y$samples$group %in% 'control', control_group, test_group)
    
    ## write out results ####
    
    if(!is.null(rename_feature)){colnames(test$genes)[1] <- rename_feature}
    df <- cbind(test$genes,df[c('logFC','logCPM','PValue','FDR','is.sig')])
    if(!is.null(fdr2) & !is.null(fc2)){
        df$is.sig2 <- (df$FDR < fdr2) * sign(df$logFC) * (abs(df$logFC) > log2(fc2))
    }
    if(report_cpm){
        cpm <- cpm(y, normalized.lib.sizes = T, log = F)
        colnames(cpm) <- paste('CPM.TMMnormalized', colnames(cpm), sep = '.')
        df <- cbind(df, cpm)
    }  
    if(report_rpkm){
        rpkm <- rpkm(y, gene.length = feature_length, normalized.lib.sizes = T, log = F)
        colnames(rpkm) <- paste('FPKM.TMMnormalized', colnames(rpkm), sep = '.')
        df <- cbind(df, rpkm)
    }  
    write.table(df,paste(out_prefix,'txt',sep = '.'), sep = '\t',quote = F,row.names = F)
    
    ## return results ####
    df_sum <- data.frame(
        output_folder = gsub('.*differential_peaks/', '', out.dir), 
        test_group = paste(test_group, collapse = '+'),
        control_group = paste(control_group, collapse = '+'),
        features_tested = nrow(y),
        features_up = sum(is.sig %in% 1),
        features_down = sum(is.sig %in% -1),
        FDR_cutoff = fdr,
        FC_cutoff = fc,
        test_samples = sum(y$samples$group %in% test_group),
        control_samples = sum(y$samples$group %in% control_group)
    )
    if(!is.null(target)){
        df_sum$target <- target  
    }
    
    if('is.sig2' %in% colnames(df)){
        df_sum <- cbind(
            df_sum,
            data.frame(
                features_up2 = sum(df$is.sig2 %in% 1),
                features_down2 = sum(df$is.sig2 %in% -1),
                FDR_cutoff2 = fdr2,
                FC_cutoff2 = fc2
            )
        )
    }
    
    return(list(summary = df_sum, y = y, df = df, design=design, test=test))
    
}


## PCA
plot_pca <- function(y, out_prefix, var.genes = NULL, color = NULL, plot_title = '', sample.label = T, feature_length = 'gene_length', scree.plot = F){
    options(stringsAsFactors = F)
    require(ggrepel)
    require(edgeR)
    
    ## Create outdir if needed
    out.dir <- dirname(out_prefix)
    if(!dir.exists(out.dir)){dir.create(out.dir,recursive = T)}
    
    if('chrom' %in% colnames(y$genes)){
        y <- y[!y$genes$chrom %in% c('chrX', 'chrY', 'chrM'), ]
    }
    
    log2rpkm <- rpkm(y, gene.length = feature_length, normalized.lib.sizes = T, log = T)
    
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
            title = plot_title,
            color = '', 
            x = paste0('PC1 (',round(pca.variance.prop[1],1),'%)'),
            y = paste0('PC2 (',round(pca.variance.prop[2],1),'%)')
        )+
        theme_bw()+
        theme(
            legend.position = 'top'
        )
    ggsave(paste(out_prefix,'PCA.pdf',sep='.'),width = 5,height = 5)
    
    ## return data
    return(p)
    
}

## Plot MD
plot_MD <- function(df, out_prefix, plot_title = ""){
    require(ggplot2)
    
    ## Create outdir if needed
    out.dir <- dirname(out_prefix)
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
            title = plot_title
        )+
        theme_bw()+
        theme(
            text = element_text(size = 8),
            legend.position = 'top'
        )
    
    ggsave(paste(out_prefix,'MD.pdf',sep = '.'),width = 4,height = 5)
    return(p)
} 

## Plot volcano 
plot_volcano <- function(df, out_prefix, plot_title = ""){
    require(ggplot2)
    
    ## Create outdir if needed
    out.dir <- dirname(out_prefix)
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
            title = plot_title
        )+
        theme_bw()+
        theme(
            text = element_text(size = 8),
            legend.position = 'top'
        )
    
    ggsave(paste(out_prefix,'Volcano.pdf',sep = '.'),width = 4,height = 5)
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
normalize_counts <- function(y, out_prefix, return = c('rpkm','cpm'), gene.length = "gene_length", log = F){
    require(edgeR)
    
    out.dir <- dirname(out_prefix)
    if (!dir.exists(out.dir)){
        dir.create(out.dir, recursive = T)
    }
    if(return[1] == 'cpm'){
        cpm <- cpm(y, normalized.lib.sizes = T, log = log)
        colnames(cpm) <- paste('CPM.TMMnormalized', colnames(cpm), sep = '.')
        df <- cbind(y$genes, cpm)
        out.suffix <- ifelse(log, 'log2CPM.txt', 'CPM.txt')
        write.table(df,paste(out_prefix, out.suffix, sep = '.'), sep = '\t',quote = F,row.names = F)
    }  
    if(return[1] == 'rpkm'){
        rpkm <- rpkm(y, gene.length = gene.length, normalized.lib.sizes = T, log = log)
        colnames(rpkm) <- paste('FPKM.TMMnormalized', colnames(rpkm), sep = '.')
        df <- cbind(y$genes, rpkm)
        out.suffix <- ifelse(log, 'log2FPKM.txt', 'FPKM.txt')
        write.table(df,paste(out_prefix, out.suffix, sep = '.'), sep = '\t',quote = F,row.names = F)
    }  
    
}

## wrapper
wrap_one_cmp <- function(y0, icmp, ss, fdr = 0.05, fc = 1.5, fdr2 = 0.01, fc2 = 2, outdir = './'){
    
    file_base <- icmp$file_base[1]
    control_group <- unlist(strsplit(as.character(icmp$control_group[[1]]), split = ';'))
    test_group <- unlist(strsplit(as.character(icmp$test_group[[1]]), split = ';'))
    plot_title <- icmp$plot_title[1]
    if (icmp$comparison_group[1] %in% colnames(ss)){
        group <- ss[,icmp$comparison_group[1]][match(rownames(y0$samples), ss$id)]
    }else if (icmp$comparison_group[1] %in% colnames(y0$samples)){
        group <- y0$samples[,icmp$comparison_group[1]]
    }else{
        group <- NULL
    }
    design_object <- as.formula(icmp$design_object[1])
    
    ## filter samples ####
    if (is.na(icmp$include_samples[1])){
        include_samples <- NULL
    }else{
        include_samples <- unlist(strsplit(icmp$include_samples[1], split = ';'))
    }
    
    if (is.na(icmp$exclude_samples[1])){
        exclude_samples <- NULL
    }else{
        exclude_samples <- unlist(strsplit(icmp$exclude_samples[1], split = ';'))
    }
    
    
    ## run DGE ####
    lst <- run_da(
        y0, 
        out_prefix = file_base,
        control_group = control_group, 
        test_group = test_group, 
        group = group, 
        include_samples = include_samples,
        exclude_samples = exclude_samples,
        feature_length = length_col,
        fdr = fdr, fc = fc, fdr2 = fdr2, fc2 = fc2,
        design_object = design_object
    )
    
    if ('error' %in% names(lst)){
        write.table(lst$error, paste0(file_base, '.error.txt'), sep = '\t', quote = F, row.names = F)
        return(lst)
    }
    
    y <- lst$y
    df <- lst$df
    lst$plots <- list(
        PCA = plot_pca(
            y, 
            file_base, 
            color = y$samples$group, 
            sample.label = T, 
            plot_title = "", 
            var.genes = 500, 
            feature_length = length_col),
        MD = plot_MD(
            df, out_prefix = file_base,
            plot_title = plot_title
        ),
        volcano = plot_volcano(
            df, out_prefix = file_base,
            plot_title = plot_title
        )
    )
    
    ##
    return(lst)
    
    
}

## add default value if a column is missing
fill_column <- function(df, colv, default.value, na.value = NULL, missing.value = NULL){
    if (!colv %in% colnames(df)){
        df[,colv] <- default.value
    }
    
    if (!is.null(na.value)){
        df[,colv] <- ifelse(is.na(df[,colv]), na.value, df[,colv])
    }
    
    if (!is.null(missing.value)){
        df[,colv] <- ifelse(df[,colv] == "", missing.value, df[,colv])
    }
    
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
if (!exists('length_col')){
    if (grepl('transcript', rds)){
        length_col <- 'EffectiveLength'
    }else{
        length_col <- 'gene_length'
    }
}
if (!exists('outdir')){
    outdir <- '.'
}

## parse comparison table ####
colnames(cmp) <- gsub(' +|\\.', '_', colnames(cmp))
cmp <- cmp %>% 
    mutate(
        test_group = gsub('-| +|&', '.', test_group),
        control_group = gsub('-| +|&', '.', control_group)
    ) %>% 
    fill_column('out_prefix', paste(gsub(';','-', cmp$test_group), gsub(';','-', cmp$control_group), sep = '_vs_')) %>% 
    fill_column('plot_title', paste(gsub(';','+', cmp$test_group), gsub(';','+', cmp$control_group), sep = ' vs ')) %>% 
    fill_column('comparison_group', '') %>% 
    fill_column('design_object', '~0+group', '~0+group', '~0+group' ) %>% 
    fill_column('include_samples', NA, NA, NA) %>% 
    fill_column('exclude_samples', NA, NA, NA) %>% 
    arrange(out_prefix, !is.na(exclude_samples), !is.na(include_samples)) %>% 
    group_by(out_prefix) %>% 
    mutate(
        n = 1:n()
    ) %>% 
    mutate(
        file_base = file.path(outdir, paste0(out_prefix, ifelse(n>1, paste0('_run',n), '')), out_prefix)
    )

## update sample sheet and DGElist ####
common_samples <- intersect(ss$id, y0$samples$id)
if (length(common_samples) == 0 ){
    stop(paste("No common samples are found between", input, "and", rds, "!"))
}
y0 <- y0[ ,y0$samples$id %in% common_samples]
ss <- ss[ss$id %in% common_samples, ]
new_meta <- setdiff(colnames(ss), colnames(y0$samples))
if (length(new_meta)>0){
    y0$samples[,new_meta] <- ss[match(y0$samples$id, ss$id),new_meta]
}


## filter gene if needed ####
if (file.exists(gene_txt)){
    genes2keep <- read.delim(gene_txt)
    y0 <- y0[y0$genes$gene_id %in% genes2keep$gene_id, ]
}

gene_types <- unlist(strsplit(gene_type, split = ','))
if(!setequal(gene_types, 'all')){
    y0 <- y0[y0$genes$gene_type %in% gene_types, ]
}

cbind(y0$genes, y0$counts) %>% 
    write.table(
        file.path(outdir, 'all_samples.genes_filtered_for_DE.raw_counts.txt'), 
        sep = '\t', quote = F, row.names = F
    )

rpkm <- rpkm(y0, gene.length = length_col, normalized.lib.sizes = T, log = F)
cbind(y0$genes, rpkm) %>% 
    write.table(
        file.path(outdir, 'all_samples.genes_filtered_for_DE.FPKM.txt'), 
        sep = '\t', quote = F, row.names = F
    )


## detect differential expression ####
de.list <- list()
for (i in 1:nrow(cmp)){
    de.list[[i]] <- wrap_one_cmp(
        y0, cmp[i,], ss, fdr, fc, fdr2, fc2, outdir = outdir
    )
}
names(de.list) <- basename(cmp$out_prefix)
saveRDS(de.list, ifelse(grepl('transcript', rds), 'differential_transcripts.rds', 'differential_genes.rds'))

