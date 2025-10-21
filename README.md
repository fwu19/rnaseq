# rnaseq

This pipeline can analyze RNAseq data derived from mRNA-seq, total RNA-seq and RNA exome libraries as well as data from Patient-derive Xenofgrapt (PDX) samples. 

## Pipeline summary

-   Build aligner index and prepare other reference files as necessary.

- Trim adapters. 

-   Split or concatenate fastq files as needed.

-   Align reads to the reference genome.

-   For PDX samples, align reads to both graft and host genomes, respectively, and filter graft alignment by removing reads of host origin.
    Save filtered fastq files if needed.

-   Run QC on raw fastq files, trimmed fastq files and read alignments.

-   Quantify expression at the gene and/or transcript level.

-   Perform differential expression analysis at the gene or transcript level.

-   Identify gene fusions.

-   Generate an analysis report including QC metrics, differential expression, description of deliverables and analytical methods.


