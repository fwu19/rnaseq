
process SAMTOOLS_VIEW {
    time = '1d'
    cpus = 8
    memory = '48G'
    //module = ['SAMtools/1.17-GCC-12.2.0']
    //container "quay.io/biocontainers/mulled-v2-227a1cb61b41a4b207e98ffab745211e900486fc:0d43472fc20149be3538449002b5af50bc71883b-0" // samtools 1.9
    container "quay.io/biocontainers/mulled-v2-31457b56c9a9c57b2aa348290faca1623bd406a7:4aff614678cdb4148839deb33790b3f304a862bd-0" // samtools v1.21

    tag "samtools on ${out_prefix}"

    publishDir "${params.outdir}/BAM/", mode: 'copy'

    input:
    tuple val(meta), val(out_prefix), path(bam), path(bai)
    
    output:
    tuple val(meta), path("${out_prefix}.{stat,flagstat}"), emit: data
    path("${out_prefix}.{stat,flagstat}")

    script:
    def suffix = task.ext.suffix ?: ""
    def args = task.ext.args ?: ""
    """
    samtools stats -@ ${task.cpus} $bam >${out_prefix}.stat
    samtools flagstat -@ ${task.cpus} -O tsv $bam >${out_prefix}.flagstat
    """
}
