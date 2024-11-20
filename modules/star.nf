
process STAR {
    time = '1d'
    cpus = 6
    memory = '36G'
    module = ['STAR/2.7.7a-GCC-10.2.0', 'SAMtools/1.11-GCC-10.2.0']


    tag "STAR on ${meta.id}"

    publishDir "${params.outdir}/STAR/$genome/${meta.id}", mode: 'copy'

    input:
    tuple val(meta), path(read1), path(read2)
    val (genome)
    path (star)
    path (gtf)
    
    output:
    tuple val(meta), path( "${meta.id}/Log.final.out" ), emit: log
    tuple val(meta), path("${meta.id}/Aligned.sortedByCoord.out.{bam,bam.bai}"), emit: bam 
    tuple val(meta), path("${meta.id}.ReadsPerGene.out.tab"), emit: counts
    path("${meta.id}", type: 'dir')

    script:
    """
    star.sh ${meta.id} ${star} ${gtf} ${task.cpus} $read1 $read2
    cp ${meta.id}/ReadsPerGene.out.tab ${meta.id}.ReadsPerGene.out.tab
    """
}
