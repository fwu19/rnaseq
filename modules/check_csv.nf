process CHECK_CSV {
    container "docker://fwu19/r-libs:4.1.2" 

    label 'process_single'

    tag "Generate $samplesheet"

    input:
    path ( samplesheet )
    path ( outdir )

    output:
    path ('samplesheet.checked.csv'), emit: csv
    path ('versions.yml'), emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def workflow = params.workflow?: "regular"
    def step = params.step?: "mapping"
    def salmon = params.run_salmon?: "false"
    def featurecounts = params.run_featurecounts?: "false"
    """
    check_csv.r csv=$samplesheet outdir=$outdir workflow=$workflow step=$step salmon=$salmon featurecounts=$featurecounts

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: \$(R --version | head -n 1)
    END_VERSIONS
    """
}
