process RSEQC {
    time = '1d'
    cpus = 6
    memory = '36G'
    module = ['R/4.1.2-foss-2020b', 'SAMtools/1.11-GCC-10.2.0']


    tag "RSeQC on ${meta.id}"

    publishDir "${params.outdir}/QC/rseqc/", pattern: "*.{txt,pdf,r,bed,log,xls,xlsx}", mode: 'copy'

    input:
    tuple val(meta), path("*.{bam,bai}", stageAs: "input/*")
    path(rseqc_bed)
    path(tx_bed)
    path(gene_bed)

    output:
    tuple val(meta), path("*.{txt,pdf,r,bed,log,xls,xlsx}"), emit: qc
    path("*.{txt,pdf,r,bed,log,xls,xlsx}")

    script:
    """
    rseqc.sh ${meta.id} ${rseqc_bed} ${tx_bed} ${gene_bed} input/*.bam 

    """
}

