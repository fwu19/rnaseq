process GTF2FASTA {
    label 'process_single'

    container = 'quay.io/biocontainers/gffread:0.12.7--h077b44d_6'
    //module = [ 'gffread/0.12.7-GCCcore-12.3.0' ]

    tag "Retrieve transcript sequences."

    input:
    path ( genome_fa )
    path ( gtf )

    output:
    path 'transcripts.fa', emit: fa
    path  "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    gffread -w transcripts.fa -g $genome_fa $gtf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gffread: \$(gffread --version 2>&1 | head -n 1)
    END_VERSIONS
    """

    stub:
    """
    touch transcripts.fa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gffread: \$(gffread --version | head -n 1)
    END_VERSIONS
    """

}
