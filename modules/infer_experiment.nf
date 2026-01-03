process INFER_EXPERIMENT {
    label 'process_single'
    module = ['RSeQC/5.0.1-foss-2021b']

    tag "Infer strandedness and read type."

    input:
    tuple val(meta), val(out_prefix), path(bam), path(bai)
    path(tx_bed)

    output:
    path ('strand.txt', emit: strand)
    path ('read_type.txt', emit: read_type)
    path ('versions.yml'), emit: versions

    script:
    def args = task.ext.args ?: "-s 200000 -q 30"
    """
    infer_experiment.sh $bam $tx_bed $args 
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        Python: \$(python --version | head -n 1)
    END_VERSIONS

    """
    
}

