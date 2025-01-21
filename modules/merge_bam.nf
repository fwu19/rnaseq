
process MERGE_BAM {
    time = '1d'
    cpus = 8
    memory = '48G'
    module = ['SAMtools/1.11-GCC-10.2.0']


    tag "merge bam files for ${out_prefix}"

    publishDir "${params.outdir}/merged_bam/", mode: 'copy'

    input:
    tuple val(meta), path("input/*")
    
    output:
    tuple val(meta), val(meta.id), path("*.bam"), emit: bam 
    tuple val(meta), val(meta.id), path("*.bai"), emit: bai
    path("*.{bam,bai}")

    script:
    def suffix = task.ext.suffix ?: ""
    def args = task.ext.args ?: ""
    """
    samtools merge ${args} -@ ${task.cpus} -f ${meta.id}${suffix}.bam input/*.bam
    samtools index ${meta.id}${suffix}.bam

    """
}
