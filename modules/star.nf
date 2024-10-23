
process STAR {
    time = '1d'
    cpus = 6
    memory = '36G'
    module = ['STAR/2.7.7a-GCC-10.2.0', 'SAMtools/1.11-GCC-10.2.0']


    tag "STAR on $sample_id"

    publishDir "${params.outdir}/STAR/$genome", mode: 'copy'

    input:
    tuple val(sample_id), path("*.fastq.gz", stageAs: "fastq/*")
    val (genome)
    path (star)
    path (gtf)
    
    output:
    tuple path("${sample_id}/"), emit: star
    tuple val(sample_id), path("${sample_id}/${sample_id}.{bam,bam.bai}"), emit: bam 
    
    script:
    """
    star.sh ${sample_id} ${star} ${gtf} ${task.cpus} 
    """
}
