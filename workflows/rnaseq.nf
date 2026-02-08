#!/usr/bin/env nextflow 

nextflow.enable.dsl=2

include { CHECK_REFERENCE } from '../subworkflows/check_reference.nf'
include { CHECK_INPUT } from '../subworkflows/check_input.nf'
include { PROCESS_FASTQ } from '../subworkflows/process_fastq.nf'
include { ALIGN_FASTQ } from '../subworkflows/align_fastq.nf'
include { CHECK_EXPERIMENT } from '../subworkflows/check_experiment.nf'
include { QUANT_GENES } from '../subworkflows/quant_genes.nf'
include { QUANT_TRANSCRIPTS } from '../subworkflows/quant_transcripts.nf'
include { QC_FASTQ } from '../subworkflows/qc_fastq.nf'
include { QC_ALIGNMENT } from '../subworkflows/qc_alignment.nf'
include { GENERATE_REPORT } from '../subworkflows/generate_report.nf'
include { WRITE_OUTPUT_CSV } from '../subworkflows/write_output_csv.nf'
include { UPDATE_OUTPUT_CSV } from '../subworkflows/update_output_csv.nf'
include { WRITE_PARAMS } from '../subworkflows/write_params.nf'

include { ARRIBA  } from '../modules/arriba.nf'

//Validate required parameters
if ( params.genome_fa == null && params.run_build_index ){
    exit 1, "Need to provide a valid path with --genome_fa path/to/genome/fasta."
}

if ( params.gtf == null && (params.run_gene_count || params.run_tx_count)){
    exit 1, "Need to provide a valid path with --gtf path/to/genes/gtf."
}

// Define step options
def step_options = ["mapping","expression_quantification","differential_expression"]
if (!(params.step in step_options)){ 
    exit 1, "Invalid option for --step. Available options: ${step_options.join(', ')}"
}

// Define update options
def update_options = [null, "qc_fastq", "qc_alignment", "gene_expression", "differential_genes", "transcript_expression", "differential_transcripts", "multiqc", "report", "csv"]
if(!(params.update in update_options )){ 
    exit 1, "Invalid option for --update. Available options: ${update_options.findAll{it != null}.join(', ')}"
}


// Define aligner options
def aligner_options = ["star", "bwa-mem"]
if(!params.aligner in aligner_options){
    exit 1, "Invalid option for --aligner. Available options: ${aligner_options.join(', ')}"
}

// Define tools to run and handle skip_tools parameter
def tools = "cutadapt,fastp,fastqc,rnaseqc,rseqc,star,bwa-mem,salmon,arriba"
def all_tools = tools.split(',')*.trim() 
if ( params.skip_tools != null ) {
    def skip_tools = params.skip_tools.split(',')*.trim()

    def use_tools = all_tools.intersect(skip_tools)

    if (!common) {
        exit 1, "Invalid option for --skip_tools. Available options: ${all_tools.split(',').join(', ')}"
    }
}

// Define a dummy file path
dummy_file = file("$projectDir/assets/dummy_file.csv")


workflow RNASEQ {

    def outdir  = new File("${params.outdir}").absolutePath
    def srcdir  = params.srcdir ? new File("${params.srcdir}").absolutePath : outdir

    ch_software_versions = Channel.empty()
    /*
    * Run reference check
    */
    gene_txt = Channel.empty()
    tx_bed = Channel.empty()
    tx_txt = Channel.empty()
    tx_fa = Channel.empty()
    collapsed_gtf = Channel.empty()

    if (params.run_reference_check) {
        CHECK_REFERENCE(
            params.genome_fa,
            params.gtf,
            params.aligner,
            params.aligner_index ?: params.aligner,
            srcdir
        )
        index_dir = CHECK_REFERENCE.out.index_dir // match params.aligner
        gene_txt = CHECK_REFERENCE.out.gene_txt
        tx_txt = CHECK_REFERENCE.out.tx_txt
        tx_bed = CHECK_REFERENCE.out.tx_bed
        tx_fa = CHECK_REFERENCE.out.tx_fa
        collapsed_gtf = CHECK_REFERENCE.out.collapsed_gtf
        ch_software_versions = ch_software_versions.mix(CHECK_REFERENCE.out.versions)
    }


    /*
    * Run input check
    */
    ch_input = Channel.fromPath( params.input, checkIfExists: true )
    if (params.run_input_check){
        CHECK_INPUT(
            params.input ? file(params.input) : dummy_file,
            params.input_dir ? file(params.input_dir) : dummy_file,
            params.metadata ? file(params.metadata) : dummy_file
        )
        samplesheet = CHECK_INPUT.out.samplesheet
        fq = CHECK_INPUT.out.fq
    ch_software_versions = ch_software_versions.mix(CHECK_INPUT.out.versions)
    }


    /*
    * cat/split fastq files
    */
    ch_reads = Channel.empty()
    ch_reads_trimmed = Channel.empty()
    ch_cutadapt_js = Channel.empty()
    ch_fastp_js = Channel.empty()
    if (params.run_process_fastq){
        PROCESS_FASTQ(
            samplesheet,
            fq,
            params.cat_fastq,
            params.trimmer,
            srcdir
        )
        ch_reads = PROCESS_FASTQ.out.reads
        ch_reads_trimmed = PROCESS_FASTQ.out.reads_trimmed
        ch_cutadapt_js = PROCESS_FASTQ.out.cutadapt_js
        ch_fastp_js = PROCESS_FASTQ.out.fastp_js
        ch_software_versions = ch_software_versions.mix(PROCESS_FASTQ.out.versions)
    }

    /*
    * Infer strandedness and read type
    */
    infer_experiment = Channel.empty()
    if (params.run_infer_experiment ){
        CHECK_EXPERIMENT(
            fq,
            index_dir,
            tx_bed,
            srcdir
        )
        infer_experiment = CHECK_EXPERIMENT.out.csv
        ch_software_versions = ch_software_versions.mix(CHECK_EXPERIMENT.out.versions)
        
    }  

    /*
    * Gene-level analysis
    */
    ch_bam_bai = Channel.empty()
    ch_star_counts = Channel.empty()
    ch_fc_counts = Channel.empty()
    ch_tx_bam = Channel.empty() // STAR's transcript bam, used for downstream quantification with Salmon
    ch_star_log = Channel.empty()
    ch_bam_bai_host = Channel.empty() 
    ch_star_log_host = Channel.empty() 
    ch_bam_bai_xeno = Channel.empty()
    ch_graft_reads = Channel.empty()
    ch_gene_expr = Channel.empty()
    ch_de = Channel.empty()
    de_csv = Channel.empty()

    if (params.gene_level){
        /*
        * run alignment 
        * for pdx workflow, align to both graft and host genomes 
        * for pdx workflow, run XenofilteR to remove reads with host origin
        */
        if (params.run_alignment){
            ALIGN_FASTQ(
                params.run_cut_adapt ? ch_reads_trimmed : ch_reads,
                params.run_split_fastq, // if true, dont run xenofilteR
                index_dir,
            )
        
            ch_bam_bai = ALIGN_FASTQ.out.bam_bai
            ch_star_counts = ALIGN_FASTQ.out.counts
            ch_tx_bam = ALIGN_FASTQ.out.tx_bam
            ch_star_log = ALIGN_FASTQ.out.star_log
            ch_bam_bai_host = ALIGN_FASTQ.out.bam_bai_host
            ch_star_log_host = ALIGN_FASTQ.out.star_log_host
            ch_bam_bai_xeno = ALIGN_FASTQ.out.bam_bai_xeno 
            ch_graft_reads = ALIGN_FASTQ.out.graft_reads
            ch_software_versions = ch_software_versions.mix(ALIGN_FASTQ.out.versions)
    
        }
    

        /*
        *   collect gene-level count matrix and call differential expression
        */
        if ( params.run_quant_genes ){
            QUANT_GENES(
                samplesheet,
                params.workflow == 'pdx' ? ch_bam_bai_xeno.ifEmpty([]) : ch_bam_bai.ifEmpty([]),
                ch_star_counts.ifEmpty([]),
                gene_txt,
                infer_experiment.ifEmpty([]),
                srcdir
            )
        
            ch_fc_counts = QUANT_GENES.out.fc_counts
            ch_gene_expr = QUANT_GENES.out.expr 
            ch_de = QUANT_GENES.out.de
            de_csv = QUANT_GENES.out.de_csv
            ch_software_versions = ch_software_versions.mix(QUANT_GENES.out.versions)

        }
    }


    /*
    * Transcript-level analysis
    */
    ch_salmon = Channel.empty()
    ch_tx_expr = Channel.empty()
    ch_dt = Channel.empty()
    dt_csv = Channel.empty()
    if (params.transcript_level ){
        /* Map to transcripts */
        // not available yet

        /*  Collect transcript-level counts and call differential transcripts */
        if (params.run_quant_transcripts){

            QUANT_TRANSCRIPTS(
                samplesheet,
                ch_tx_bam.ifEmpty([]),
                tx_fa,
                tx_txt,
                gene_txt,
                srcdir
            )
            ch_salmon = QUANT_TRANSCRIPTS.out.salmon
            ch_tx_expr = QUANT_TRANSCRIPTS.out.tx_expr
            ch_dt = QUANT_TRANSCRIPTS.out.dt
            dt_csv = QUANT_TRANSCRIPTS.out.dt_csv
            ch_software_versions = ch_software_versions.mix(QUANT_TRANSCRIPTS.out.versions)
        }
    }


    /*
    * Identify gene fusions
    */
    if (params.gene_fusion){
        // can build a subworkflow for gene_fusion analysis
        if (params.run_arriba){
            ARRIBA(
            params.workflow == 'pdx' ? ch_graft_reads : (params.run_cut_adapt ? ch_reads_trimmed : ch_reads), 
            params.genome, 
            index_dir, 
            params.gtf,
            params.genome_fa,
            params.blacklist,
            params.known_fusions,
            params.protein_domains
            )
            ch_software_versions = ch_software_versions.mix(ARRIBA.out.versions)
        }
    }


    /*
    * QC fastq
    */
    ch_fastqc = Channel.empty()
    ch_fastqc_trimmed = Channel.empty()
    if (params.run_qc_fastq){
        QC_FASTQ(
            ch_reads,
            ch_reads_trimmed
        )
        ch_fastqc = QC_FASTQ.out.fastqc.map{it[1]}.flatten().collect()
        ch_fastqc_trimmed = QC_FASTQ.out.fastqc_trimmed.map{it[1]}.flatten().collect()
        ch_software_versions = ch_software_versions.mix(QC_FASTQ.out.versions)
    }

    /*
    * QC alignment
    */
    ch_rnaseqc = Channel.empty()
    ch_rseqc = Channel.empty()
    ch_hs_metrics = Channel.empty()
    ch_bam_stat = Channel.empty()
    ch_bam_stat_host = Channel.empty()
    ch_bam_stat_xeno = Channel.empty()

    if (params.run_qc_alignment){
        // need to add a parser to align_fastq.csv
        QC_ALIGNMENT(
            ch_bam_bai,
            ch_bam_bai_host,
            ch_bam_bai_xeno,
            collapsed_gtf,
            tx_bed,
            infer_experiment.ifEmpty([]),
            srcdir
        )
        ch_rnaseqc = QC_ALIGNMENT.out.rnaseqc.flatten().collect()
        ch_rseqc = QC_ALIGNMENT.out.rseqc.flatten().collect()
        ch_hs_metrics = QC_ALIGNMENT.out.hs_metrics.flatten().collect()
        ch_bam_stat = QC_ALIGNMENT.out.bam_stat.flatten().collect()
        ch_bam_stat_host = QC_ALIGNMENT.out.bam_stat_host.flatten().collect()
        ch_bam_stat_xeno = QC_ALIGNMENT.out.bam_stat_xeno.flatten().collect()
        ch_software_versions = ch_software_versions.mix(QC_ALIGNMENT.out.versions)
        
    }


    /*
    * generate multiqc report and analysis report
    */
    report_rmd = Channel.empty()
    GENERATE_REPORT(
        samplesheet, 
        ch_star_counts.map{it[2]}.flatten().collect().ifEmpty([]), 
        ch_fc_counts.map{it[2]}.flatten().collect().ifEmpty([]), 
        ch_salmon.map{it[2]}.flatten().collect().ifEmpty([]), 
        ch_fastqc.flatten().collect().ifEmpty([]), 
        ch_cutadapt_js.flatten().collect().ifEmpty([]),
        ch_fastp_js.flatten().collect().ifEmpty([]),
        ch_fastqc_trimmed.flatten().collect().ifEmpty([]),
        ch_star_log.flatten().collect().ifEmpty([]),
        ch_rnaseqc.flatten().collect().ifEmpty([]), 
        ch_rseqc.flatten().collect().ifEmpty([]), 
        ch_bam_stat.flatten().collect().ifEmpty([]), 
        ch_bam_stat_host.flatten().collect().ifEmpty([]), 
        ch_bam_stat_xeno.flatten().collect().ifEmpty([]), 
        ch_hs_metrics.flatten().collect().ifEmpty([]), 
        ch_gene_expr.ifEmpty([]), 
        ch_de.ifEmpty([]), 
        ch_tx_expr.ifEmpty([]), 
        ch_dt.ifEmpty([]),
        srcdir
    )
    report_rmd = GENERATE_REPORT.out.rmd
    ch_software_versions = ch_software_versions.mix(GENERATE_REPORT.out.versions)

    /*
    * Write output file paths
    */
    if (params.update == 'csv'){
        UPDATE_OUTPUT_CSV(
            samplesheet
        )
    }else{
        WRITE_OUTPUT_CSV(
            ch_bam_bai.ifEmpty([]),
            ch_bam_bai_host.ifEmpty([]),
            ch_bam_bai_xeno.ifEmpty([]),
            ch_star_counts.ifEmpty([]),
            ch_fc_counts.ifEmpty([]),
            ch_salmon.ifEmpty([])
        )
    }
       

    /*
    * Write software versions to a yml file and overridden params to a json file
    */
    WRITE_PARAMS(
        ch_software_versions
    )
}

////////////////////////////////////////////////////
/* --              COMPLETION            -- */
////////////////////////////////////////////////////

