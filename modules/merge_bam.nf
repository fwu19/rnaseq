
process MERGE_BAM {
    time = '1d'
    cpus = 8
    memory = '48G'
    //module = ['SAMtools/1.17-GCC-12.2.0']
    container "quay.io/biocontainers/mulled-v2-a4c30dc1a2dfc3f31070c6a8acc1c627f7a22916:da999000e91310fca6d5021998dab12999a6ad0c-0" // gatk 4.1.4 samtools 1.9

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
