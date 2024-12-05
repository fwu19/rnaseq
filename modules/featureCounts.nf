process FEATURECOUNTS {
    time = '1d'
    cpus = 6
    memory = '36G'
    module = ['Subread/2.0.0-GCC-8.3.0']


    tag "featureCounts on ${meta.id}"

    publishDir "${params.outdir}/featureCounts/", pattern: "${meta.id}*.{txt,summary}", mode: 'copy'

    input:
    tuple val(meta), val(out_prefix), path( bam )
    path (gtf)
    val (read_type)
    val (strand)
    
    output:
    tuple val(meta), path("${out_prefix}*.txt"), emit: counts
    path("*.{txt,summary}")

    script:
    """
    featureCounts.sh ${out_prefix} ${gtf} ${bam} ${read_type} ${strand} ${task.cpus}
    """
}