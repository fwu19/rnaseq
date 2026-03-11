#!/usr/bin/env Rscript

options(stringsAsFactors = F)
options(scipen = 99)
library(dplyr)
library(ggplot2)
library(patchwork)

## functions ####

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

## for all
run_da <- function(
    y0, file_base, 
    control_group, test_group, group=NULL, 
    exclude_samples = NULL, include_samples = NULL,
    fdr=0.01, fc=2, fdr2=NULL, fc2=NULL, 
    report_cpm=F, report_rpkm=T, 
    TMM=T, method = 'QL', 
    rename_feature = NULL, feature_length = 'gene_length', 
    model_formula = ~0+group,
    target = NULL
){
    require(edgeR)
    
    ## create output directory ####
    out.dir <- dirname(file_base)
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
    design <- model.matrix(model_formula,data = y$samples)
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
    write.table(df,paste(file_base,'txt',sep = '.'), sep = '\t',quote = F,row.names = F)
    
    ## return results ####
    df_sum <- data.frame(
        output_folder = basename(out.dir), 
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
    
    # if(sample.label){
    #     p <- p +
    #         geom_text_repel(size = 2.4, color = 'black', position = 'jitter',max.overlaps = 80)
    # }
    
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
    wrap_plots(
      p, 
      p +
        geom_text_repel(size = 2.4, color = 'black', position = 'jitter',max.overlaps = 80),
      nrow = 1
    )
    ggsave(paste(out_prefix,'PCA.pdf',sep='.'),width = 10,height = 5)
    
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
    model_formula <- as.formula(icmp$model_formula[1])
    
    ## filter samples ####
    if (is.na(icmp$include_samples[1])){
        include_samples <- NULL
    }else{
        include_samples <- trimws(unlist(strsplit(icmp$include_samples[1], split = ';')))
    }
    
    if (is.na(icmp$exclude_samples[1])){
        exclude_samples <- NULL
    }else{
        exclude_samples <- trimws(unlist(strsplit(icmp$exclude_samples[1], split = ';')))
    }
    
    
    ## run DGE ####
    lst <- run_da(
        y0, 
        file_base = file_base,
        control_group = control_group, 
        test_group = test_group, 
        group = group, 
        include_samples = include_samples,
        exclude_samples = exclude_samples,
        feature_length = length_col,
        fdr = fdr, fc = fc, fdr2 = fdr2, fc2 = fc2,
        model_formula = model_formula
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
# input=${samplesheet} comparison=$c{omparison} gene_txt=${gene_txt} count_file=${count_file} fdr=${params.fdr} fc=${params.fc} fdr2=${params.fdr2} fc2=${params.fc2} gene_type=${params.de_gene_type} length_col=${length_col}
args <- as.vector(commandArgs(T)) 
for (arg in strsplit(args, split = '=')){
  assign(trimws(arg[1]), trimws(arg[2]))
}
# for (arg in args){eval(parse(text = arg))}

ss <- read.csv(input, colClasses = 'character', comment.char = '#') %>% 
    unique.data.frame() %>% 
    mutate(
      sample_group = gsub('-| +|&', '.', sample_group)
  )

if(grepl('dummy', comparison)){
    cat(comparison, "is a dummy file! Provide --comparison path/to/comparison_file (a comparison table in csv, txt, tsv or rds format)!")
    quit()
}else if(grepl('.csv$', comparison)){
    cmp <- read.csv(comparison, colClasses = 'character', comment.char = '#')
}else if(grepl('.tsv$|.txt$', comparison)){
    cmp <- read.delim(comparison, colClasses = 'character', comment.char = '#')
}else if (grepl('.rds$', comparison)){
    cmp <- readRDS(comparison)
}else {
    stop('--comparison only takes csv, txt, tsv or rds file.')
}

fdr <- as.numeric(fdr)
fc <- as.numeric(fc)
fdr2 <- as.numeric(fdr2)
fc2 <- as.numeric(fc2)
if (!exists('length_col')){
    if (grepl('transcript', count_file)){
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
colnames(cmp) <- gsub('design.object|design_object', 'model_formula', colnames(cmp))
cmp <- cmp %>% 
    mutate(
        test_group = gsub('-| +|&', '.', test_group),
        control_group = gsub('-| +|&', '.', control_group)
    ) %>% 
    fill_column(
      'out_prefix', 
      paste(gsub(';','-', cmp$test_group), gsub(';','-', cmp$control_group), sep = '_vs_'),
      paste(gsub(';','-', cmp$test_group), gsub(';','-', cmp$control_group), sep = '_vs_'),
      paste(gsub(';','-', cmp$test_group), gsub(';','-', cmp$control_group), sep = '_vs_')
      ) %>% 
    fill_column(
      'plot_title', 
      paste(gsub(';','+', cmp$test_group), gsub(';','+', cmp$control_group), sep = ' vs '),
      paste(gsub(';','+', cmp$test_group), gsub(';','+', cmp$control_group), sep = ' vs '),
      paste(gsub(';','+', cmp$test_group), gsub(';','+', cmp$control_group), sep = ' vs ')
      ) %>% 
    fill_column('comparison_group', NA, NA, NA) %>% 
    fill_column('model_formula', '~0+group', '~0+group', '~0+group' ) %>% 
    fill_column('include_samples', NA, NA, NA) %>% 
    fill_column('exclude_samples', NA, NA, NA) %>% 
    fill_column('run_version', "", "", "") %>% 
    arrange(out_prefix, !is.na(exclude_samples), !is.na(include_samples)) %>% 
    group_by(out_prefix, run_version) %>% 
    mutate(
      n = 1:n(),
      run_version = ifelse(n == 1, run_version, ifelse(run_version == "", paste0('run', n), paste0(run_version, '.run', n))),
      file_base = file.path(outdir, ifelse(run_version == '', out_prefix, paste(out_prefix, run_version, sep = '.')), out_prefix),
      include_samples = gsub(' +', '.', include_samples),
      exclude_samples = gsub(' +', '.', exclude_samples)
    )


## create y0 if not available ####
if (grepl('rds$', count_file)){
  y0 <- readRDS(count_file)
}else{
  nsamples <- nrow(ss)
  gene_expr <- read.delim(count_file)
  cts <- gene_expr[,tail(1:ncol(gene_expr), nsamples)]
  mid <- gsub('-', '.', ifelse(grepl('^[0-9]', ss$id), paste0('X', ss$id), ss$id))
  colnames(cts) <- ss$id[match(colnames(cts), mid)] # revert "." to "-"
  
  y0 <- count2dgelist(
    counts = cts, 
    features = gene_expr[, 1:(ncol(gene_expr)-nsamples)], 
    samples = ss %>% arrange(factor(id, levels = colnames(cts))),
    group = 'sample_group'
  )
  rm(nsamples, mid, gene_expr, cts)
}

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
    if (!'gene_id' %in% colnames(genes2keep)){
      stop(paste('No gene_id column is found in', gene_txt))
    }
    y0 <- y0[y0$genes$gene_id %in% genes2keep$gene_id, ]
}

gene_types <- unlist(strsplit(gene_type, split = ','))
if(!setequal(gene_types, 'all') & 'gene_type' %in% colnames(y0$genes)){
    y0 <- y0[y0$genes$gene_type %in% gene_types, ]
}

cbind(y0$genes, y0$counts) %>% 
    write.table(
        file.path(outdir, 'all_samples.genes_filtered_for_DE.raw_counts.txt'), 
        sep = '\t', quote = F, row.names = F
    )

cpm <- cpm(y0, normalized.lib.sizes = T, log = F)
cbind(y0$genes, cpm) %>% 
  write.table(
    file.path(outdir, 'all_samples.genes_filtered_for_DE.CPM.txt'), 
    sep = '\t', quote = F, row.names = F
    )

rpkm <- rpkm(y0, gene.length = length_col, normalized.lib.sizes = T, log = F)
cbind(y0$genes, rpkm) %>% 
    write.table(
        file.path(outdir, 'all_samples.genes_filtered_for_DE.FPKM.txt'), 
        sep = '\t', quote = F, row.names = F
    )

tpm <- rpkm/matrix(colSums(rpkm), nrow = nrow(rpkm), ncol = ncol(rpkm), byrow = T)*1e6
cbind(y0$genes, tpm) %>% 
  write.table(
    file.path(outdir, 'all_samples.genes_filtered_for_DE.TPM.txt'), 
    sep = '\t', quote = F, row.names = F
    )


## detect differential expression ####
de.list <- list()
for (i in 1:nrow(cmp)){
    de.list[[i]] <- wrap_one_cmp(
        y0, cmp[i,], ss, fdr, fc, fdr2, fc2, outdir = outdir
    )
}
names(de.list) <- basename(dirname(cmp$file_base))
saveRDS(de.list, ifelse(grepl('transcript', count_file), 'differential_transcripts.rds', 'differential_genes.rds'))

## write a copy of comparison table ####
cmp %>% 
  dplyr::select(-file_base) %>%
  relocate(test_group, control_group, out_prefix, run_version, model_formula, comparison_group, exclude_samples, include_samples, plot_title) %>% 
  write.table(
    ifelse(grepl('transcript', count_file), 'comparisons.differential_transcripts.csv', 'comparisons.differential_genes.csv'),
    sep = ',', quote = F, row.names = F
    )

## write out DE summary ####
de_sum <- bind_rows(lapply(de.list, function(ide){ide$summary})) %>% 
  dplyr::relocate(test_samples, control_samples, .after = last_col()) %>% 
  mutate(output_folder = basename(output_folder))
out.txt <- '_README.txt'
headers <- c(
  '# For each comparison, differentially expressed genes or transcripts are defined using two sets of FDR and fold change cutoffs, which correspond to columns is.sig and is.sig2 in [test group]_vs_[control group]/[test group]_vs_[control group].txt.', 
  '# 1 = upregulation in test group, -1 = downregulation in test group, 0 = no significance.', 
  '# Plots are based on the first cutoffs, i.e. column is.sig.')
writeLines(headers, '_README.txt')
de_sum %>% 
  write.table(
    out.txt, append = T,
    sep = '\t', quote = F, row.names = F)