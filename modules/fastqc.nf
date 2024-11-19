process FASTQC {
    time = '1d'
    cpus = 1
    memory = '12G'
    module = 'FastQC/0.11.9-Java-11'


    tag "FASTQC on ${meta.id}"

    publishDir "${params.outdir}/QC/fastqc/", pattern: '*.{html,zip}', mode: 'copy'

    input:
    tuple val(meta), path(read1), path(read2)
    

    output:
    tuple val(meta), path("*.{html,zip}"), emit: qc

    script:
    """
    fastqc -o ./ --casava $read1 $read2
    """
}
