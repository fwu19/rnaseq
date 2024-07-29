params.outdir = "Analysis"

process FASTQC {
    time = '1d'
    cpus = 1
    memory = '12G'
    module = 'FastQC/0.11.9-Java-11'


    tag "FASTQC on $sample_id"

    publishDir "${outdir}/FastQC/", pattern: '*.{html,zip}', mode: 'copy'

    input:
    tuple val(sample_id), path(reads1), path(reads2)
    path(outdir)

    output:
    tuple path("*.{html,zip}")

    script:
    """
    fastqc -o ./ --casava ${reads1} ${reads2}
     """
}
