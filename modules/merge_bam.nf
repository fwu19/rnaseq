
process MERGE_BAM {
    label "process_high"
    module = ['SAMtools/1.17-GCC-12.2.0']


    tag "merge bam files for ${out_prefix}"

    publishDir "${params.outdir}/merged_bam/", mode: 'copy'

    input:
    tuple val(meta), path("input/*")
    
    output:
    tuple val(meta), val(meta.id), path("${meta.id}*.bam"), emit: bam 
    tuple val(meta), val(meta.id), path("${meta.id}*.bam.bai"), emit: bai
    path("${meta.id}*.{bam,bam.bai}")

    script:
    def suffix = task.ext.suffix ?: ""
    def args = task.ext.args ?: ""
    """
    samtools merge ${args} -@ ${task.cpus} -f -o ${meta.id}${suffix}.bam input/*.bam
    samtools index ${meta.id}${suffix}.bam

    """
}
