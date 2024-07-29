params.outdir = "Analysis"

process TEST {
    time = '1d'
    cpus = 1
    memory = '5G'


    tag "test on $sample_id"

    publishDir "${outdir}/TEST/${sample_id}", overwrite: true

    input:
    tuple val(sample_id), path(output)
    path(outdir)

    output:
    tuple val(sample_id), path("out.txt")

    script:
    """
    echo $sample_id, $output > out.txt
     """
}