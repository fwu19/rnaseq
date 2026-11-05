#!/usr/bin/env nextflow 

nextflow.enable.dsl=2

include { BUILD_INDEX } from '../subworkflows/build_index.nf'
include { CHECK_REFERENCE } from '../subworkflows/check_reference.nf'
include { CHECK_REFERENCE as CHECK_REFERENCE_HOST } from '../subworkflows/check_reference.nf'
include { CHECK_INPUT } from '../subworkflows/check_input.nf'
include { PROCESS_FASTQ } from '../subworkflows/process_fastq.nf'
include { CHECK_EXPERIMENT } from '../subworkflows/check_experiment.nf'
include { ALIGN_FASTQ } from '../subworkflows/align_fastq.nf'
include { ALIGN_FASTQ as ALIGN_FASTQ_HOST } from '../subworkflows/align_fastq.nf'
include { PDX_XENOFILTER } from '../subworkflows/pdx_xenofilter.nf'
include { QUANT_GENES } from '../subworkflows/quant_genes.nf'
include { QUANT_TRANSCRIPTS } from '../subworkflows/quant_transcripts.nf'
include { QC_FASTQ } from '../subworkflows/qc_fastq.nf'
include { QC_ALIGNMENT } from '../subworkflows/qc_alignment.nf'
include { DETECT_FUSIONS } from '../subworkflows/detect_fusions.nf'
include { GENERATE_REPORT } from '../subworkflows/generate_report.nf'
include { WRITE_OUTPUT_CSV } from '../subworkflows/write_output_csv.nf'
include { WRITE_PARAMS } from '../subworkflows/write_params.nf'

include { XENGSORT_CLASSIFY } from '../modules/xengsort_classify.nf'

//Validate required parameters
if ( params.genome_fa == null && params.run_build_index ){
    exit 1, "Need to provide a valid path with --genome_fa path/to/genome/fasta."
}

if ( params.gtf == null && (params.run_gene_count || params.run_tx_count)){
    exit 1, "Need to provide a valid path with --gtf path/to/genes/gtf."
}

// Define step options
def step_options = ["mapping", "expression_quantification", "differential_expression", "gene_fusion", "multiqc" ]
if (!(params.step in step_options)){ 
    exit 1, "Invalid option for --step. Available options: ${step_options.join(', ')}"
}

// Define update options
def update_options = [null, "qc_alignment", "gene_level", "gene_expression", "differential_genes", "transcript_level", "transcript_expression", "differential_transcripts", "gene_fusion", "multiqc", "report", "output"]
if(!(params.update in update_options )){ 
    exit 1, "Invalid option for --update. Available options: ${update_options.findAll{it != null}.join(', ')}"
}


// Define aligner options
def aligner_options = ["star", "bwa-mem"]
if(!params.aligner in aligner_options){
    exit 1, "Invalid option for --aligner. Available options: ${aligner_options.join(', ')}"
}

// Define tools to run and handle skip_tools parameter
def tools = "cutadapt,fastp,fastqc,rnaseqc,rseqc,star,bwa-mem,salmon,arriba,xengsort,xenofilter"
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

    /*
    * Write software versions to a yaml file and overridden params to a json file
    */
    WRITE_PARAMS()

    def outdir  = new File("${params.outdir}").absolutePath
    def srcdir  = params.srcdir ? new File("${params.srcdir}").absolutePath : outdir

    ch_software_versions = Channel.empty()

    /*
    * Build index only
    */
    if (params.run_build_index){
        BUILD_INDEX(
            params.genome_fa,
            params.gtf,
            params.genome_fa_host,
            params.gtf_host,
            params.reference_index ?: params.aligner
        )
        ch_software_versions = ch_software_versions.mix(BUILD_INDEX.out.versions)
    }

    /*
    * Run reference check
    */
    index_dir = Channel.empty()
    gene_txt = Channel.empty()
    tx_bed = Channel.empty()
    tx_txt = Channel.empty()
    tx_fa = Channel.empty()
    collapsed_gtf = Channel.empty()
    index_host_dir = Channel.empty()
    gene_txt_host = Channel.empty()
    tx_bed_host = Channel.empty()
    tx_txt_host = Channel.empty()
    tx_fa_host = Channel.empty()
    collapsed_gtf_host = Channel.empty()
    xengsort_dir = Channel.empty()

    if (params.run_reference_check) {
        CHECK_REFERENCE(
            params.genome,
            params.genome_fa,
            params.gtf,
            params.aligner,
            params.star,
            params.bwa,
            params.rnaseqc_gtf,
            params.rseqc_bed,
            params.gene_txt,
            params.tx_bed,
            params.tx_fa,
            params.tx_txt,
            srcdir
        )
        index_dir = CHECK_REFERENCE.out.index_dir // match params.aligner
        gene_txt = CHECK_REFERENCE.out.gene_txt
        tx_txt = CHECK_REFERENCE.out.tx_txt
        tx_bed = CHECK_REFERENCE.out.tx_bed
        tx_fa = CHECK_REFERENCE.out.tx_fa
        collapsed_gtf = CHECK_REFERENCE.out.collapsed_gtf
        ch_software_versions = ch_software_versions.mix(CHECK_REFERENCE.out.versions)

        if (params.workflow == 'pdx'){
            CHECK_REFERENCE_HOST(
                params.genome_host,
                params.genome_fa_host,
                params.gtf_host,
                params.aligner,
                params.star_host,
                params.bwa_host,
                params.rnaseqc_gtf_host,
                params.rseqc_bed_host,
                params.gene_txt_host,
                params.tx_bed_host,
                params.tx_fa_host,
                params.tx_txt_host,
                srcdir
            )
            index_host_dir = CHECK_REFERENCE_HOST.out.index_dir
            gene_txt_host = CHECK_REFERENCE_HOST.out.gene_txt
            tx_txt_host = CHECK_REFERENCE_HOST.out.tx_txt
            tx_bed_host = CHECK_REFERENCE_HOST.out.tx_bed
            tx_fa_host = CHECK_REFERENCE_HOST.out.tx_fa
            collapsed_gtf_host = CHECK_REFERENCE_HOST.out.collapsed_gtf
            xengsort_dir = CHECK_REFERENCE_HOST.out.xengsort_dir
            ch_software_versions = ch_software_versions.mix(CHECK_REFERENCE_HOST.out.versions)

        }
    }


    /*
    * Run input check
    */
    samplesheet = Channel.empty()
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
    ch_trimgalore_js = Channel.empty()
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
        ch_trimgalore_js = PROCESS_FASTQ.out.trimgalore_js
        ch_cutadapt_js = PROCESS_FASTQ.out.cutadapt_js
        ch_fastp_js = PROCESS_FASTQ.out.fastp_js
        ch_software_versions = ch_software_versions.mix(PROCESS_FASTQ.out.versions)
    }

    /*
    * For PDX samples, classify graft and host reads
    */
    ch_reads_xengsort = Channel.empty()
    ch_reads_xengsort_host = Channel.empty()
    ch_xengsort_log = Channel.empty()
    if (params.run_classify_fastq && params.run_xengsort){
        XENGSORT_CLASSIFY(
            params.run_cut_adapt ? ch_reads_trimmed : ch_reads,
            xengsort_dir
        )
        ch_reads_xengsort = XENGSORT_CLASSIFY.out.graft_fq
        ch_reads_xengsort_host = XENGSORT_CLASSIFY.out.host_fq
        ch_xengsort_log = XENGSORT_CLASSIFY.out.log
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
    * run alignment 
    * for pdx workflow, align to both graft and host genomes and run XenofilteR to remove reads with host origin
    */
    ch_bam_bai = Channel.empty()
    ch_star_counts = Channel.empty()
    ch_tx_bam = Channel.empty() // STAR's transcript bam, used for downstream quantification with Salmon
    ch_star_log = Channel.empty()

    ch_bam_bai_host = Channel.empty() 
    ch_bam_bai_xeno = Channel.empty()
    ch_star_counts_host = Channel.empty()
    ch_tx_bam_host = Channel.empty()
    ch_star_log_host = Channel.empty() 
    ch_graft_reads = Channel.empty()

    if (params.run_alignment){
            ALIGN_FASTQ(
                params.workflow == 'pdx' && params.pdx_tool == 'xengsort' ? ch_reads_xengsort : (params.run_cut_adapt ? ch_reads_trimmed : ch_reads),
                index_dir,
                params.gtf,
                params.workflow == 'pdx' && params.pdx_tool == 'xengsort' ? ch_reads_xengsort_host : (params.run_cut_adapt ? ch_reads_trimmed : ch_reads),
                index_host_dir,
                params.gtf_host
            )
        
            ch_bam_bai = ALIGN_FASTQ.out.bam_bai
            ch_star_counts = ALIGN_FASTQ.out.counts
            ch_tx_bam = ALIGN_FASTQ.out.tx_bam
            ch_star_log = ALIGN_FASTQ.out.star_log
            ch_bam_bai_host = ALIGN_FASTQ.out.bam_bai_host
            ch_star_counts_host = ALIGN_FASTQ.out.counts_host
            ch_tx_bam_host = ALIGN_FASTQ.out.tx_bam_host
            ch_star_log_host = ALIGN_FASTQ.out.star_log_host
            ch_software_versions = ch_software_versions.mix(ALIGN_FASTQ.out.versions)

            if (params.workflow == 'pdx' && params.run_filter_fastq && params.run_xenofilter){
                PDX_XENOFILTER(
                        params.run_cut_adapt ? ch_reads_trimmed : ch_reads,
                        index_dir,
                        index_host_dir,
                        ch_bam_bai,
                        ch_bam_bai_host
                )
                ch_bam_bai_xeno = PDX_XENOFILTER.out.bam_bai_xeno
                ch_graft_reads = PDX_XENOFILTER.out.graft_reads
                ch_software_versions = ch_software_versions.mix(PDX_XENOFILTER.out.versions)

            }
    
    }

    /*
    *   collect gene-level count matrix and call differential expression
    */
    ch_fc_counts = Channel.empty()
    ch_fc_summary = Channel.empty()
    ch_gene_expr = Channel.empty()
    ch_de = Channel.empty()
    de_csv = Channel.empty()

    ch_fc_counts_host = Channel.empty()
    ch_fc_summary_host = Channel.empty()
    ch_gene_expr_host = Channel.empty()
    ch_de_host = Channel.empty()
    de_csv_host = Channel.empty()

    if (params.run_quant_genes){
        QUANT_GENES(
                samplesheet,
                params.run_xenofilter ? ch_bam_bai_xeno : ch_bam_bai,
                ch_star_counts,
                params.de_gene_txt ? file(params.de_gene_txt, checkIfExists: true) : gene_txt,
                ch_bam_bai_host,
                ch_star_counts_host,
                params.de_gene_txt_host ? file(params.de_gene_txt_host, checkIfExists: true) : gene_txt_host,
                srcdir,
                infer_experiment
        )
        
        ch_fc_counts = QUANT_GENES.out.fc_counts
        ch_fc_summary = QUANT_GENES.out.fc_summary
        ch_gene_expr = QUANT_GENES.out.expr 
        ch_de = QUANT_GENES.out.de
        de_csv = QUANT_GENES.out.de_csv
        ch_fc_counts_host = QUANT_GENES.out.fc_counts_host
        ch_fc_summary_host = QUANT_GENES.out.fc_summary_host
        ch_gene_expr_host = QUANT_GENES.out.expr_host
        ch_de_host = QUANT_GENES.out.de_host
        de_csv_host = QUANT_GENES.out.de_csv_host
        ch_software_versions = ch_software_versions.mix(QUANT_GENES.out.versions)        

    }


    /*
    *  Collect transcript-level counts and call differential transcripts
    */
    ch_salmon = Channel.empty()
    ch_tx_expr = Channel.empty()
    ch_dt = Channel.empty()
    dt_csv = Channel.empty()

    ch_salmon_host = Channel.empty()
    ch_tx_expr_host = Channel.empty()
    ch_dt_host = Channel.empty()
    dt_csv_host = Channel.empty()

    if (params.run_quant_transcripts){
        QUANT_TRANSCRIPTS(
                samplesheet,
                ch_tx_bam,
                tx_fa,
                tx_txt,
                params.de_gene_txt ? file(params.de_gene_txt, checkIfExists: true) : gene_txt,
                ch_tx_bam_host,
                tx_fa_host,
                tx_txt_host,
                params.de_gene_txt_host ? file(params.de_gene_txt_host, checkIfExists: true) : gene_txt_host,
                srcdir
        )
        ch_salmon = QUANT_TRANSCRIPTS.out.salmon
        ch_tx_expr = QUANT_TRANSCRIPTS.out.tx_expr
        ch_dt = QUANT_TRANSCRIPTS.out.dt
        dt_csv = QUANT_TRANSCRIPTS.out.dt_csv
        ch_salmon_host = QUANT_TRANSCRIPTS.out.salmon_host
        ch_tx_expr_host = QUANT_TRANSCRIPTS.out.tx_expr_host
        ch_dt_host = QUANT_TRANSCRIPTS.out.dt_host
        dt_csv_host = QUANT_TRANSCRIPTS.out.dt_csv_host
        ch_software_versions = ch_software_versions.mix(QUANT_TRANSCRIPTS.out.versions)

    }

    /*
    * Identify gene fusions
    */
    ch_arriba = Channel.empty()
    if (params.run_fusion){
        DETECT_FUSIONS(
            params.run_xenofilter ? ch_graft_reads : (params.run_xengsort ? ch_reads_xengsort : (params.run_cut_adapt ? ch_reads_trimmed : ch_reads)), 
            params.genome, 
            index_dir, 
            params.gtf,
            params.genome_fa,
            params.blacklist,
            params.known_fusions,
            params.protein_domains
        )
        ch_arriba = DETECT_FUSIONS.out.arriba
        ch_software_versions = ch_software_versions.mix(DETECT_FUSIONS.out.versions)
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
        ch_fastqc = QC_FASTQ.out.fastqc
        ch_fastqc_trimmed = QC_FASTQ.out.fastqc_trimmed
        ch_software_versions = ch_software_versions.mix(QC_FASTQ.out.versions)
    }

    /*
    * QC alignment
    */
    ch_rnaseqc = Channel.empty()
    ch_rseqc = Channel.empty()
    ch_hs_metrics = Channel.empty()
    ch_bam_stat = Channel.empty()

    ch_rnaseqc_host = Channel.empty()
    ch_rseqc_host = Channel.empty()
    ch_hs_metrics_host = Channel.empty()
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
            collapsed_gtf_host,
            tx_bed_host,
            infer_experiment,
            srcdir
        )
        ch_rnaseqc = QC_ALIGNMENT.out.rnaseqc
        ch_rseqc = QC_ALIGNMENT.out.rseqc
        ch_hs_metrics = QC_ALIGNMENT.out.hs_metrics
        ch_bam_stat = QC_ALIGNMENT.out.bam_stat
        ch_rnaseqc_host = QC_ALIGNMENT.out.rnaseqc_host
        ch_rseqc_host = QC_ALIGNMENT.out.rseqc_host
        ch_hs_metrics_host = QC_ALIGNMENT.out.hs_metrics_host
        ch_bam_stat_host = QC_ALIGNMENT.out.bam_stat_host
        ch_bam_stat_xeno = QC_ALIGNMENT.out.bam_stat_xeno
        ch_software_versions = ch_software_versions.mix(QC_ALIGNMENT.out.versions)


        
    }


    /*
    * generate multiqc report and analysis report
    */
    report_rmd = Channel.empty()
    GENERATE_REPORT(
        samplesheet, 
        ch_fastqc, 
        ch_cutadapt_js,
        ch_fastp_js,
        ch_trimgalore_js,
        ch_fastqc_trimmed,
        ch_xengsort_log,
        ch_star_log,
        ch_star_log_host,
        ch_rnaseqc,
        ch_rnaseqc_host,
        ch_rseqc,
        ch_rseqc_host,
        ch_bam_stat,
        ch_bam_stat_host,
        ch_bam_stat_xeno,
        ch_hs_metrics, 
        ch_star_counts, 
        ch_star_counts_host,
        ch_fc_summary, 
        ch_fc_summary_host,
        ch_gene_expr, 
        ch_gene_expr_host,
        ch_de, 
        ch_de_host,
        ch_salmon, 
        ch_salmon_host,
        ch_tx_expr, 
        ch_tx_expr_host,
        ch_dt,
        ch_dt_host,
        ch_arriba,
        srcdir
    )
    report_rmd = GENERATE_REPORT.out.rmd
    ch_software_versions = ch_software_versions.mix(GENERATE_REPORT.out.versions)

    /*
    * Write output file paths
    */
    WRITE_OUTPUT_CSV(
            samplesheet,
            srcdir,
            ch_bam_bai,
            ch_bam_bai_host,
            ch_bam_bai_xeno,
            ch_tx_bam,
            ch_tx_bam_host,
            ch_star_counts,
            ch_star_counts_host,
            ch_fc_counts,
            ch_fc_counts_host,
            ch_salmon,
            ch_salmon_host,
            ch_arriba
        )

    
    /* Collect software versions */
    ch_software_versions
        .collectFile(storeDir: "${params.outdir}/pipeline_info/${params.info_dir}/", name: 'software_versions.yml', sort: false, newLine: true)

}

////////////////////////////////////////////////////
/* --              COMPLETION            -- */
////////////////////////////////////////////////////

