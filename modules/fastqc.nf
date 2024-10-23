process FASTQC {
    time = '1d'
    cpus = 1
    memory = '12G'
    module = 'FastQC/0.11.9-Java-11'


    tag "FASTQC on $sample_id"

    publishDir "${params.outdir}/QC/fastqc/", pattern: '*.{html,zip}', mode: 'copy'

    input:
    tuple val(sample_id), path(read1), path(read2)
    

    output:
    tuple path("*.{html,zip}"), emit: qc

    script:
    """
    fastqc -o ./ --casava $read1 $read2
    """
}
