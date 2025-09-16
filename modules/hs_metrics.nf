process HS_METRICS {
    time = '1d'
    cpus = 1
    memory = '24G'
    //module = ['GATK/4.4.0.0-GCCcore-12.2.0-Java-17', 'SAMtools/1.17-GCC-12.2.0']
    container "quay.io/biocontainers/mulled-v2-a4c30dc1a2dfc3f31070c6a8acc1c627f7a22916:da999000e91310fca6d5021998dab12999a6ad0c-0" // gatk 4.1.4 samtools 1.9

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
    collect_hs_metrics.sh ${meta.id} $bam $genome_fa $target_region

    """
}


	