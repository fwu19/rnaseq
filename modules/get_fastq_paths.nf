process GET_FASTQ_PATHS {
    module = ['fhR/4.1.2-foss-2021b']

    label 'process_single'

    tag "Get paths to fastq files."

    input:
    path ( "fastq/" )

    output:
    path 'input.csv'        , emit: csv
    path  "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    get_fastq_paths.r fastq/*

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: \$(R --version | head -n 1)
    END_VERSIONS
    """
}
