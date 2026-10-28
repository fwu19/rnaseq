process COLLAPSE_GTF {
    label 'process_single'

    container 'docker://broadinstitute/gtex_rnaseq:V10'
    //container 'quay.io/biocontainers/bx-python:0.15.1--py313haa19386_0' 
    //module = [ "fhPython/3.11.3-foss-2023a", "SciPy-bundle/2025.07-gfbf-2024a-numpy-2", "bx-python/0.13.0-foss-2023b"]

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
    unset PYTHONPATH
    export PYTHONNOUSERSITE=1

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
