process INFER_STRAND {
    label 'process_single'
    module = ['RSeQC/5.0.1-foss-2021b']

    tag "Infer strandedness."

    input:
    tuple val(meta), val(out_prefix), path(bam), path(bai)
    path(tx_bed)

    output:
    path ('strand.txt', emit: strand)

    script:
    def args = task.ext.args ?: "-s 200000 -q 30"
    """
    infer_experiment.sh $bam $tx_bed $args > strand.txt
    
    """
    
}

