process CHECK_INPUT {
    container "docker://fwu19/r-libs:4.1.2" 

    label 'process_single'

    tag "Generate $samplesheet"

    input:
    path ( samplesheet )
    path ( metadata )

    output:
    path ('samplesheet.valid.csv'), emit: csv
    path ('fq.csv'), emit: fq
    path  "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    check_input.r $samplesheet $metadata

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: \$(R --version | head -n 1)
    END_VERSIONS
    """
}
