process COLLAPSE_GTF {
    module = [ "fhPython/3.11.3-foss-2023a", "bx-python/0.13.0-foss-2023b"]

    label 'process_single'

    tag "Process a GTF file."

    input:
    path ( gtf )

    output:
    path 'genes.collapsed.gtf', emit: gtf
    path  "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    collapse_annotation.py $gtf genes.collapsed.gtf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        Python: \$(python --version | head -n 1)
    END_VERSIONS
    """

    stub:
    """
    touch genes.collapsed.gtf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        Python: \$(python --version | head -n 1)
    END_VERSIONS
    """

}
