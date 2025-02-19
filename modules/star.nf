
process STAR {
    time = '1d'
    cpus = 8
    memory = '48G'
    module = ['STAR/2.7.7a-GCC-10.2.0', 'SAMtools/1.11-GCC-10.2.0']


    tag "STAR on ${out_prefix}"

    publishDir "${params.outdir}/STAR/$genome/${meta.id}", mode: 'copy'

    input:
    tuple val(meta), val(out_prefix), path(read1), path(read2)
    val (genome)
    path (star_index)
    path (gtf)
    
    output:
    tuple val(meta), val(out_prefix), path( "output/${out_prefix}.bam" ), emit: bam 
    tuple val(meta), val(out_prefix), path( "output/${out_prefix}.bam.bai" ), emit: bai
    tuple val(meta), val(out_prefix), path( "output/${out_prefix}.ReadsPerGene.out.tab" ), emit: counts
    tuple val(meta), val(out_prefix), path( "output/${out_prefix}.Log.final.out" ), emit: log
    path("${out_prefix}/", type: 'dir')

    script:
    def args = task.ext.args ?: ""
    """
    star.sh ${out_prefix} ${star_index} ${gtf} ${task.cpus} $read1 $read2 "$args"
    mkdir output
    cd output
    ln -s ../${out_prefix}/Aligned.sortedByCoord.out.bam ${out_prefix}.bam
    ln -s ../${out_prefix}/Aligned.sortedByCoord.out.bam.bai ${out_prefix}.bam.bai
    ln -s ../${out_prefix}/ReadsPerGene.out.tab ${out_prefix}.ReadsPerGene.out.tab
    ln -s ../${out_prefix}/Log.final.out ${out_prefix}.Log.final.out
    """
}
