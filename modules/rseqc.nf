process RSEQC {
    label "process_high"
    module = ['RSeQC/5.0.1-foss-2021b', 'SAMtools/1.17-GCC-12.2.0', 'R/4.1.2-foss-2021b']


    tag "RSeQC on ${meta.id}"

    publishDir "${params.outdir}/QC/rseqc/", pattern: "*.{txt,pdf,r,bed,log,xls,xlsx}", mode: 'copy'

    input:
    tuple val(meta), val(out_prefix), path(bam), path(bai)
    path(tx_bed)

    output:
    tuple val(meta), path("${out_prefix}*.{txt,pdf,r,bed,log,xls,xlsx}"), emit: qc
    path("${out_prefix}*.{txt,pdf,r,bed,log,xls,xlsx}")
    path ('versions.yml'), emit: versions

    script:
    """
    rseqc.sh ${out_prefix} ${bam} ${tx_bed} 

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$( samtools --version | head -n 1 | sed -e "s/.* //g" )
        R: \$(R --version | head -n 1)
        Python: \$(python --version | head -n 1)
    END_VERSIONS


    """
}

