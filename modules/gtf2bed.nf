process GTF2BED {
    container "docker://quay.io/biocontainers/mulled-v2-5c3b70e8ce9caf9b6b1999ab1fa0390894d4b503:2a7ee357fe425dbd45d1b771a8dcadedc6b7242a-0"
    label 'process_single'

    tag "Convert gtf to bed12."

    input:
    path ( gtf )

    output:
    path '*.bed'        , emit: bed
    path  "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    gtf2bed $gtf > transcripts.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        Perl: \$(perl --version | head -n 1)
    END_VERSIONS
    """

    stub:
    """
    touch transcripts.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        Perl: \$(perl --version | head -n 1)
    END_VERSIONS
    """

}
