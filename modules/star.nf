
process STAR {
    time = '1d'
    cpus = 6
    memory = '36G'
    module = ['STAR/2.7.7a-GCC-10.2.0', 'SAMtools/1.11-GCC-10.2.0']


    tag "STAR on $sample_id"

    publishDir "${params.outdir}/STAR/$genome/$sample_id", mode: 'copy'

    input:
    tuple val(sample_id), path(read1), path(read2)
    val (genome)
    path (star)
    path (gtf)
    
    output:
    tuple path("*"), emit: star
    tuple val(sample_id), path("${sample_id}.{bam,bam.bai}"), emit: bam 
    tuple path("${sample_id}.ReadsPerGene.out.tab"), emit: counts
    
    script:
    """
    star.sh ${sample_id} ${star} ${gtf} ${task.cpus} $read1 $read2
    mv ReadsPerGene.out.tab ${sample_id}.ReadsPerGene.out.tab
    rm $read1 $read2
    """
}
