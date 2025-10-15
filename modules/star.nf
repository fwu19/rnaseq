
process STAR {
    label "process_high"
    module = ['STAR/2.7.10b-GCC-12.2.0', 'SAMtools/1.17-GCC-12.2.0']


    tag "STAR on ${out_prefix}"

    publishDir "${params.outdir}/STAR/$genome/${meta.id}", mode: 'copy'

    input:
    tuple val(meta), val(out_prefix), path(read1), path(read2)
    val (genome)
    path (star_index)
    path (gtf)
    
    output:
    tuple val(meta), val(out_prefix), path( "${out_prefix}.bam" ), emit: bam 
    tuple val(meta), val(out_prefix), path( "${out_prefix}.bam.bai" ), emit: bai
    tuple val(meta), val(out_prefix), path( "output/${out_prefix}.toTranscriptome.out.bam" ),emit: tx_bam 
    tuple val(meta), val(out_prefix), path( "output/${out_prefix}.ReadsPerGene.out.tab" ), emit: counts
    tuple val(meta), val(out_prefix), path( "output/${out_prefix}.Log.final.out" ), emit: log
    path("${out_prefix}/", type: 'dir')
    path("${out_prefix}.{bam,bam.bai}")

    script:
    def args = task.ext.args ?: ""
    def rg = task.ext.rg ?: "none"
    """
    star.sh ${out_prefix} ${star_index} ${gtf} ${task.cpus} $read1 $read2 "$rg" "$args"

    mv ${out_prefix}/Aligned.sortedByCoord.out.bam ${out_prefix}.bam
    mv ${out_prefix}/Aligned.sortedByCoord.out.bam.bai ${out_prefix}.bam.bai

    mkdir output
    cd output
    ln -s ../${out_prefix}/Aligned.toTranscriptome.out.bam ${out_prefix}.toTranscriptome.out.bam
    ln -s ../${out_prefix}/ReadsPerGene.out.tab ${out_prefix}.ReadsPerGene.out.tab
    ln -s ../${out_prefix}/Log.final.out ${out_prefix}.Log.final.out
    """
}
