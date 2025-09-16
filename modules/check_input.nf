process CHECK_INPUT {

    label 'process_single'

    tag "Generate $samplesheet"

    container "docker://fwu19/r-libs:4.1.2" 

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
    R -e "sessionInfo()" > session_info.txt
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: \$(R --version | head -n 1)
    END_VERSIONS
    """
}
