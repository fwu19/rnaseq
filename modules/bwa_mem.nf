
process BWA_MEM {
    label "process_high"
    module = ['BWA/0.7.17-GCCcore-12.2.0', 'SAMtools/1.17-GCC-12.2.0']


    tag "BWA_MEM on ${out_prefix}"

    publishDir "${params.outdir}/BWA/$genome/${meta.id}", mode: 'copy'

    input:
    tuple val(meta), val(out_prefix), path(read1), path(read2)
    val (genome)
    path (bwa_index)
    
    output:
    tuple val(meta), val(out_prefix), path( "${out_prefix}.bam" ), emit: bam 
    tuple val(meta), val(out_prefix), path( "${out_prefix}.bam.bai" ), emit: bai
    path("${out_prefix}.{bam,bam.bai}")

    script:
    def args = task.ext.args ?: ""
    """
    bwa_mem.sh ${task.cpus} ${bwa_index} ${read1} ${read2} ${out_prefix}
    
    """
}
