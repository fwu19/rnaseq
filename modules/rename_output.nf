process RENAME_OUTPUT {

    label 'process_single'

    tag "Rename output files from old pipelines to match v0.4."

    input:
    val (srcdir)
    val (outdir)

    script:
    """
    rename_output.sh ${params.workflow} $srcdir $outdir ${params.star_dir} ${params.star_host_dir} ${params.star_xeno_dir} ${params.bwa_dir} ${params.bwa_host_dir} ${params.bwa_xeno_dir} ${params.featurecounts_dir} ${params.salmon_dir} ${params.arriba_dir} ${params.fastqc_dir} ${params.fastqc_trimmed_dir} ${params.cutadapt_dir} ${params.fastp_dir} ${params.star_log_dir} ${params.star_log_host_dir} ${params.rnaseqc_dir} ${params.rseqc_dir} ${params.gatk_dir} ${params.samtools_dir} ${params.samtools_host_dir} ${params.samtools_xeno_dir} ${params.multiqc_dir} ${params.gene_expr_dir} ${params.tx_expr_dir} ${params.de_dir} ${params.dt_dir} 
    """
}
