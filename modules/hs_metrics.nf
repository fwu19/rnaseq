process HS_METRICS {
    time = '1d'
    cpus = 1
    memory = '24G'
    module = ['GATK/4.1.8.1-GCCcore-8.3.0-Java-1.8.0_181', 'SAMtools/1.10-GCCcore-8.3.0']


    tag "collect hs_metrics on ${meta.id}"

    publishDir "${params.outdir}/QC/gatk/", pattern: "*.{tsv,gct}", mode: 'copy'

    input:
    tuple val(meta), path("*.{bam,bai}", stageAs: "input/*")
    val(genome_fa)
    path(target_region)

    output:
    tuple val(meta), path("*.hs_metrics.txt"), emit: qc

    script:
    """
    collect_hs_metrics.sh ${meta.id} input/Aligned.sortedByCoord.out.bam $genome_fa $target_region

    """
}


	