process HS_METRICS {
    label "process_medium"
    module = ['GATK/4.4.0.0-GCCcore-12.2.0-Java-17', 'SAMtools/1.17-GCC-12.2.0']


    tag "collect hs_metrics on ${meta.id}"

    publishDir "${params.outdir}/QC/gatk/", pattern: "*.{tsv,gct}", mode: 'copy'

    input:
    tuple val(meta), val(out_prefix), path(bam), path(bai)
    val(genome_fa)
    path(target_region)

    output:
    tuple val(meta), path("${out_prefix}*.hs_metrics.txt"), emit: qc
    tuple val(meta), path("${out_prefix}*.hs_metrics.txt")

    script:
    """
    collect_hs_metrics.sh ${out_prefix} $bam $genome_fa $target_region

    """
}


	