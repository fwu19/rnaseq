process FEATURECOUNTS {
    time = '1d'
    cpus = 6
    memory = '36G'
    module = ['Subread/2.0.0-GCC-8.3.0']


    tag "featureCounts on ${meta.id}"

    publishDir "${params.outdir}/featureCounts/", pattern: "${meta.id}*.{txt,summary}", mode: 'copy'

    input:
    tuple val(meta), path( bam )
    path (gtf)
    val (read_type)
    val (strand)
    
    output:
    tuple val(meta), path("${meta.id}*.txt"), emit: counts
    tuple val(meta), path("*")

    script:
    """
    featureCounts.sh ${meta.id} $gtf *.bam ${read_type} ${strand} ${task.cpus}
    """
}