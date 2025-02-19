process RNASEQC {
    time = '1d'
    cpus = 6
    memory = '36G'
    module = ['RNA-SeQC/2.3.4-foss-2019b', 'SAMtools/1.11-GCC-10.2.0']


    tag "RNA-SeQC on ${meta.id}"

    publishDir "${params.outdir}/QC/rnaseqc/", pattern: "*.{tsv,gct}", mode: 'copy'

    input:
    tuple val(meta), val(out_prefix), path(bam)
    tuple val(meta), val(out_prefix), path(bai)
    path(gtf)
    val(strand)
    val(read_type)

    output:
    tuple val(meta), path("*.{tsv,gct}"), emit: qc
    path("*.{tsv,gct}")
    
    script:
    """
    rnaseqc.sh ${out_prefix} ${gtf} ${strand} ${read_type} ${bam} 

    """
}


	