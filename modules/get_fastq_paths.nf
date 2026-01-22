process GET_FASTQ_PATHS {
    container "docker://fwu19/r-libs:4.1.2" 

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
    def args = task.ext.args ?: ''
    """
    get_fastq_paths.r $args r1_pattern="${params.r1_pattern}" r2_pattern="${params.r2_pattern}"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: \$(R --version 2>&1 | head -n 1 | sed -e "s/R version //g; s/ .*//g" )
    END_VERSIONS
    """
}
