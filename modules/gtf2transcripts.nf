process GTF2TRANSCRIPTS {
    container "docker://fwu19/r-libs:4.1.2" 

    label 'process_single'

    tag "Process a GTF file."

    input:
    path ( gtf )

    output:
    path '*.txt', emit: txt
    path  "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    gtf2transcripts.r $gtf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: \$(R --version 2>&1 | head -n 1 | sed -e "s/R version //g; s/ .*//g" )
    END_VERSIONS
    """

    stub:
    """
    touch genes.bed genes.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: \$(R --version | head -n 1)
    END_VERSIONS
    """

}
