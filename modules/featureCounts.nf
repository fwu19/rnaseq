process FEATURECOUNTS {
    time = '1d'
    cpus = 6
    memory = '36G'
    module = ['Subread/2.0.0-GCC-8.3.0']


    tag "featureCounts on $sample_id"

    publishDir "${params.outdir}/featureCounts/", pattern: "${sample_id}*.{txt,summary}", mode: 'copy'

    input:
    tuple val(sample_id), path(bam)
    path (gtf)
    
    output:
    tuple val(sample_id), path("*")
    tuple path("${sample_id}*.txt"), emit: counts

    script:
    """
    featureCounts.sh $sample_id $gtf ${sample_id}.bam ${params.read_type} ${params.strand} ${task.cpus}
    """
}