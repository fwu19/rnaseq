process GTF2GENES {
    container "docker://fwu19/r-libs:4.1.2" 

    label 'process_single'

    tag "Retrieve genes from a GTF file."

    input:
    path ( gtf )

    output:
    path '*.txt', emit: txt
    path  "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    gtf2genes.r $gtf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: \$(R --version | head -n 1)
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
