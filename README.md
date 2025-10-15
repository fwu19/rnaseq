# rnaseq

This pipeline can analyze RNAseq data derived from mRNA-seq, total RNA-seq and RNA exome libraries as well as data from Patient-derive Xenofgrapt (PDX) samples. 

## Pipeline summary

- Trim adapters. 

- Split or concatenate fastq files. 

- Align reads to the reference genome. 

- For PDX samples, align reads to both graft and host genomes, respectively, and filter graft alignment by removing reads of host origin. 

- Run QC on raw fastq files, trimmed fastq files and read alignments. 

- Quantify expression at the gene level. 

- Quantify expression at the transcript level. 

- Perform differential expression analysis at the gene or transcript level. 

- Identify gene fusions.

- Generate an analysis report including QC metrics, differential expression, description of deliverables and analytical methods. 

