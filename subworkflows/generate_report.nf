/*
* Generate analysis report 
*/

include { RENDER_REPORT } from '../modules/render_report.nf'

workflow GENERATE_REPORT {
    take:
    samplesheet
    multiqc
    gene_expr
    de
    tx_expr
    dt
    hs_metrics
    report_dir

    main: 
    ch_versions  = Channel.empty()

    /*
    indir = "${params.outdir}/MultiQC/multiqc_data/"
    if (multiqc.isEmpty() && file(indir).exists() ){ multiqc = Channel.fromPath(indir, type: 'dir').ifEmpty([]) }

    infile = "${params.outdir}/expression_quantification/all_samples.gene_raw_counts.txt"
    if (gene_expr.isEmpty() && file(infile).exists()){ gene_expr = Channel.fromPath(infile).ifEmpty([]) }

    indir = "${params.outdir}/differential_genes/"
    if ( de.isEmpty() && file(indir, type: 'dir').exists()){ de = Channel.fromPath(indir, type: 'dir').ifEmpty([]) }

    infile = "${params.outdir}/expression_quantification/all_samples.transcript_raw_counts.txt"
    if (tx_expr.isEmpty() && file(infile).exists()){ tx_expr = Channel.fromPath(infile).ifEmpty([]) }

    indir = "${params.outdir}/differential_transcripts/"
    if (dt.isEmpty() && file(indir, type: 'dir').exists()){ dt = Channel.fromPath(indir, type: 'dir').ifEmpty([])}

    indir = "${params.outdir}/QC/gatk/"
    if (hs_metrics.isEmpty() && file(indir).exists()){ hs_metrics = Channel.fromPath("${indir}/*", type: 'dir').collect().ifEmpty([]) }
    */

    RENDER_REPORT(
        samplesheet, 
        multiqc, 
        gene_expr, 
        de, 
        tx_expr, 
        dt, 
        hs_metrics, 
        report_dir
    )
    ch_versions = RENDER_REPORT.out.versions

    emit:
    versions = ch_versions

}
