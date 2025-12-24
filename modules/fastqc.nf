process FASTQC {
    label "process_single"
    module = 'FastQC/0.12.1-Java-11'


    tag "FASTQC on ${meta.id}"

    publishDir "${params.outdir}/QC/fastqc/", pattern: '*.{html,zip}', mode: 'copy'

    input:
    tuple val(meta), val(out_prefix), path("${out_prefix}_R1.fastq.gz"), path("${out_prefix}_R2.fastq.gz")
    

    output:
    tuple val(meta), path("*.{html,zip}"), emit: qc

    script:
    """
    fastqc -o ./ --casava ${out_prefix}_R1.fastq.gz ${out_prefix}_R2.fastq.gz
    """
}
