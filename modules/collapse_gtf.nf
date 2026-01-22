process COLLAPSE_GTF {
    label 'process_single'

//    container = 'quay.io/biocontainers/mulled-v2-b09526e9169de63faa8bb47ddeb97043262b8aaf:7ff5bd4607a0f73068087f8a88f5d13031c7f58d-0'
    module = [ "fhPython/3.11.3-foss-2023a", "bx-python/0.13.0-foss-2023b"]

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
        Python: \$(python --version 2>&1 | head -n 1 | sed -e "s/.* //g" )
    END_VERSIONS
    """

}
