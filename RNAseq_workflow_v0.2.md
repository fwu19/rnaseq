$$\\\\\[0.2in\]$$

# 

## Introduction

This pipeline can analyze RNAseq data derived from mRNA-seq, total
RNA-seq and RNA exome libraries as well as data from Patient-derived
Xenofgrapt (PDX) samples.

### Pipeline summary

-   Trim adapters.

-   Split or concatenate fastq files as needed.

-   Align reads to the reference genome.

-   For PDX samples, align reads to both graft and host genomes,
    respectively, and filter graft alignment by removing reads of host
    origin. Save filtered fastq files if needed.

-   Run QC on raw fastq files, trimmed fastq files and read alignments.

-   Quantify expression at the gene level.

-   Quantify expression at the transcript level.

-   Perform differential expression analysis at the gene or transcript
    level.

-   Identify gene fusions.

-   Generate an analysis report including QC metrics, differential
    expression, description of deliverables and analytical methods.

## Usage

To launch the pipeline, set up a run script using the template below:

``` bash
#!/usr/bin/env bash
module purge all; 
module load Nextflow/23.04.2; 
module load Singularity; export PATH=$SINGULARITYROOT/bin/:$PATH; 
repo="fwu19/rnaseq -r v0.2.0"
config=path/to/nextflow.config # see examples in tab "Resouces for rhino/gizmo users"
params=path/to/params.json # see examples in tab "Resouces for rhino/gizmo users"
report=nextflow_report.html

sh="nextflow run $repo -c $config -params-file $params -with-report $report $@" # $@ for taking command arguments
sbatch --wrap="$sh"
```

$$\\\\\[0.2in\]$$

Save the above script as `run.sh` and run one of the following command
lines:

``` bash
bash run.sh \
  --input path/to/sample_sheet.csv \
  --metadata path/to/metadata.csv \
  --comparison path/to/comparisons.csv

# --input_dir takes a path to fastq files
bash run.sh \
  --input_dir path/to/fastq/directory/ \
  --metadata path/to/metadata.csv \
  --comparison path/to/comparisons.csv

# if multiple flow cells are analyzed together, can point to a text file listing paths to flowcells, one path per line
bash run.sh \
  --input_dir path/to/flowcell_list.txt \
  --metadata path/to/metadata.csv \
  --comparison path/to/comparisons.csv
```

Alternatively, command line arguments (starting with “--”) can also be
specified in params.json.  

## Parameters

### Input/Output options

#### `--input` and `--input_dir`

Either one is required and only one should be supplied.

  

#### `--input_dir path/to/fastq/directory` 

This retrieves fastq files from a directory and create a sample sheet
with columns id, fastq_1 and fastq_2. id is obtained from fastq file
names, e.g. **H_1**\_S50_L004_R1_001.fastq.gz. The underlying codes
assume that fastq file names follow the pattern
`id_S[0-9]+(_L[0-9]+)?_R[12]`  
- If fastq files are located in different folders, can point to a text
file listing all paths, one path per line  

  

#### `--input path/to/input.csv`

A comma-delimited file with **.csv** extension . It requires 3 columns:
id, fastq_1 and fastq_2.  
**id** = sample ID, use for naming output files.  
**fastq_1** = path/to/read1.fastq.gz  
**fastq_2** = path/to/read2.fastq.gz.  

For samples sequenced across multiple lanes, provide each pair of fastq
files in one row.  
Additional columns will either be ignored or taken as metadata if column
names are recognized (see `--metadata` below), but will be overwritten
if the same column names also exit in `--metadata`.  

An example below:  
  

<div class="datatables html-widget html-fill-item-overflow-hidden html-fill-item" id="htmlwidget-4c4a03d31cd919ee27be" style="width:100%;height:auto;"></div>
<script type="application/json" data-for="htmlwidget-4c4a03d31cd919ee27be">{"x":{"filter":"none","vertical":false,"extensions":["Buttons"],"data":[["H1048_EV_1","H1048_EV_2","H1048_BMP_1","H1048_BMP_2"],["H1048_EV_1_S115_L006_R1_001.fastq.gz","H1048_EV_2_S116_L005_R1_001.fastq.gz","H1048_BMP_1_S118_L005_R1_001.fastq.gz","H1048_BMP_2_S119_L005_R1_001.fastq.gz"],["H1048_EV_1_S115_L006_R2_001.fastq.gz","H1048_EV_2_S116_L005_R2_001.fastq.gz","H1048_BMP_1_S118_L005_R2_001.fastq.gz","H1048_BMP_2_S119_L005_R2_001.fastq.gz"],["group1","group1","group2","group2"]],"container":"<table class=\"cell-border stripe\">\n  <thead>\n    <tr>\n      <th>id<\/th>\n      <th>fastq_1<\/th>\n      <th>fastq_2<\/th>\n      <th>sample_group<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"dom":"Bfrtip","buttons":["csv"],"columnDefs":[],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>

  

#### `--metadata path/to/metadata.csv`

A comma-delimted file with **.csv** extension . Each row is one sample.
It requires a minimum of 1 column: id. Additional columns are optional.
If column *sample_group* exists, it will be used to group samples in QC
and differential expression analysis.  
**sample_group** = condition, cell line, etc.  

An example below:  
  

<div class="datatables html-widget html-fill-item-overflow-hidden html-fill-item" id="htmlwidget-3b22640f360a5df6b038" style="width:100%;height:auto;"></div>
<script type="application/json" data-for="htmlwidget-3b22640f360a5df6b038">{"x":{"filter":"none","vertical":false,"extensions":["Buttons"],"data":[["H1048_BMP_1","H1048_BMP_2","H1048_EV_1","H1048_EV_2"],["group2","group2","group1","group1"]],"container":"<table class=\"cell-border stripe\">\n  <thead>\n    <tr>\n      <th>id<\/th>\n      <th>sample_group<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"dom":"Bfrtip","buttons":["csv"],"columnDefs":[],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>

-   When `--metadata` is supplied, only samples (by column id) in both
    `--metadata` and `--input`/`--input_dir` will be analyzed. This is
    useful for analyzing only a subset of samples.  

-   If necessary metadata information is included in `--input`, there is
    no need to provide `--metadata`.  

  

#### `--comparison path/to/comparison_table.{csv,txt,tsv,rds}` 

A table listing comparisons of interest. It requires 2 columns:
control.group and test.group, Each row corresponds to a comparison.
Values in both columns should be listed under column *sample_group* in
–metadata or –input. The acceptable file formats are comma-delimited
(.csv), tab-delimited (.txt, .tsv) or a data.frame in R (.rds).  

An example below:  
  

<div class="datatables html-widget html-fill-item-overflow-hidden html-fill-item" id="htmlwidget-f6d67475e6d448c9b7e8" style="width:100%;height:auto;"></div>
<script type="application/json" data-for="htmlwidget-f6d67475e6d448c9b7e8">{"x":{"filter":"none","vertical":false,"extensions":["Buttons"],"data":[["group1"],["group2"]],"container":"<table class=\"cell-border stripe\">\n  <thead>\n    <tr>\n      <th>control.group<\/th>\n      <th>test.group<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"dom":"Bfrtip","buttons":["csv"],"columnDefs":[],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>

-   It is OK to list all possible pairs of control and test groups. The
    workflow performs a comparison only when data are available from
    both control and test groups.  

-   By default, output folders will be named as
    *test.group_vs_control.group/* and plot title as *“test.group vs
    control.group (referene)”*. If needed, can specify them by adding
    columns *out.prefix* and *plot.title* into the table.  

$$\\\\\[0.1in\]$$

### Pipeline options

#### `--workflow`

Possible values: regular, exome, pdx. Default: regular.

  

#### `--run_arriba`

Run Arriba. Default: false

  

#### `--run_cat_fastq`

Concatenate multiple fastq files from the same library Default: false

  

#### `--run_cut_adapt`

Run adapter. Default: true

  

#### `--run_de`

Run differential gene expression. Default: true

  

#### `--run_dt`

Run differential transcript expression. Default: false

  

#### `--run_featurecounts`

Run FeatureCounts. Default: true if –workflow pdx, and false otherwise.

  

#### `--run_multiqc`

Run MultiQC. Default: true

  

#### `--run_report`

Generate analysis report. Default: true

  

#### `--run_salmon`

Run Salmon. Default: false

  

#### `--run_split_fastq`

Split a fastq file into multiple smaller files. Default: false

  

#### `--split_size`

The number of reads to split a fastq file by. Default: 50000000

  

$$\\\\\[0.1in\]$$

### QC options

#### `--run_qc`

If false, do not run any QC tools. Default: true

  

#### `--run_hs_metrics`

Run GATK function CollectHsMetrics. Default: true if –workflow exome,
and false otherwise.

  

#### `--run_fastqc`

Run FastqC. Default: true

  

#### `--run_rnaseqc`

Run RNA-SeQC. Default: true

  

#### `--run_rseqc`

Run RSeQC. Default: true

  

$$\\\\\[0.1in\]$$

### Flow switching options

#### `--only_alignment`

Only run until STAR alignment. Default: false.

  

#### `--only_filter_fastq`

Only run until removing reads of host origin if –workflow pdx. Default:
false.

  

#### `--only_input`

Only run until generating the sample sheet. Default: false

  

#### `--only_merge_fastq`

Only run until merging fastq files if –run_cat_fastq true. Default:
false.

  

#### `--only_split_fastq`

Only run until split fastq files if –run_split_fastq true. Default:
false.

  

## Resources for rhino/gizmo users

### Example params.json and nextflow.config

\`/fh/fast/\_SR/Genomics/proj/fwu/nextflow/rnaseq/assets/params_config/’

`nextflow.config` - set up working environment. env.R_LIBS and
params.local_assets are shared paths. Other paths may point to users’
space.  
`regular/` - used for mRNAseq or total RNAseq.  
`exome/` - used for RNA exome. It additionally runs Picard function
CollectHsMetrics to examine metrics such as off-target rates.  
`pdx/` - used for PDX samples. It additionally aligns reads to the host
genome and filters the alignment to the graft genome by removing reads
of host origin.  

$$\\\\\[0.1in\]$$

### Demo analyses

`/fh/fast/_SR/Genomics/proj/fwu/nextflow/rnaseq/demo/` has multiple
subdirectories, each demonstrating a different workflow:  
`regular/` - mRNAseq or total RNAseq.  
`exome/` - RNA exome.  
`pdx/` - PDX samples.  

  

#### Deliverables

-   `params.outdir` (Default: Analysis/) saves the analysis results.  

-   `params.outdir/00_RNAseq_analysis_report.html` provides a summary of
    results, description of deliverables, methods and references.  

$$\\\\\[0.1in\]$$

#### Workflow reports and debugging

-   `./nextflow_report.html` - nextflow report (upon workflow completion
    or failure) in the launch directory. This file contains a
    comprehensive summary of the workflow run.  

-   `${params.outdir}/pipeline_info/samplesheet.valid.csv` - the
    formatted samplesheet (merging `--input` and `--metadata`). This
    file can be used for a rerun,
    e.g. `bash run.sh --input samplesheet.valid.csv`, in case original
    files for –input and –metadata are not available.  

-   `${params.outdir}/pipeline_info/execution_trace.txt` - This file
    saves the cached locations of each analysis.  

-   `./saved_data/` - This directory (saved under the launch directory)
    contains a R markdown file and necessary data files to regenerate
    `${params.outdir}/00_RNAseq_analysis_report.html` in case a
    modification is needed, e.g. changing figure heights and widths, or
    adding/deleting a section.  

$$\\\\\[0.2in\]$$
