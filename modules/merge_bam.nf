
process MERGE_BAM {
    label "process_high"

    container = 'quay.io/biocontainers/samtools:1.17--hd87286a_1'
    //module = ['SAMtools/1.17-GCC-12.2.0']

    tag "merge bam files for ${out_prefix}"

    publishDir "${params.outdir}/merged_bam/", mode: 'copy'

    input:
    tuple val(meta), path("input/*")
    
    output:
    tuple val(meta), val(meta.id), path("${meta.id}*.bam"), emit: bam 
    tuple val(meta), val(meta.id), path("${meta.id}*.bam.bai"), emit: bai
    path("${meta.id}*.{bam,bam.bai}")
    path ('versions.yml'), emit: versions

    script:
    def suffix = task.ext.suffix ?: ""
    def args = task.ext.args ?: ""
    """
    samtools merge ${args} -@ ${task.cpus} -f -o ${meta.id}${suffix}.bam input/*.bam
    samtools index ${meta.id}${suffix}.bam
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$( samtools --version | head -n 1 | sed -e "s/.* //g" )
    END_VERSIONS

    """
}
