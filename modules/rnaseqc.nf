process RNASEQC {
    time = '1d'
    cpus = 6
    memory = '36G'
    //module = ['RNA-SeQC/2.4.2-foss-2021b', 'SAMtools/1.17-GCC-12.2.0']
    container "quay.io/biocontainers/rna-seqc:2.4.2--h8492097_0" // RNA-SeQC 2.4.2

    tag "RNA-SeQC on ${meta.id}"

    publishDir "${params.outdir}/QC/rnaseqc/", pattern: "*.{tsv,gct}", mode: 'copy'

    input:
    tuple val(meta), val(out_prefix), path(bam), path(bai)
    path(gtf)
    val(strand)
    val(read_type)

    output:
    tuple val(meta), path("${out_prefix}*.{tsv,gct}"), emit: qc
    path("${out_prefix}*.{tsv,gct}")
    
    script:
    """
    rnaseqc.sh ${out_prefix} ${gtf} ${strand} ${read_type} ${bam} 

    """
}


	