process RSEQC {
    time = '1d'
    cpus = 6
    memory = '36G'
    module = ['RSeQC/5.0.1-foss-2021b', 'SAMtools/1.17-GCC-12.2.0', 'R/4.1.2-foss-2021b']


    tag "RSeQC on ${meta.id}"

    publishDir "${params.outdir}/QC/rseqc/", pattern: "*.{txt,pdf,r,bed,log,xls,xlsx}", mode: 'copy'

    input:
    tuple val(meta), val(out_prefix), path(bam)
    tuple val(meta), val(out_prefix), path(bai)
    path(rseqc_bed)
    path(tx_bed)
    path(gene_bed)

    output:
    tuple val(meta), path("${out_prefix}*.{txt,pdf,r,bed,log,xls,xlsx}"), emit: qc
    path("${out_prefix}*.{txt,pdf,r,bed,log,xls,xlsx}")

    script:
    """
    rseqc.sh ${out_prefix} ${rseqc_bed} ${tx_bed} ${gene_bed} $bam 

    """
}

