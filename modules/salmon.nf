
process SALMON {
    time = '1d'
    cpus = 8
    memory = '48G'
    //module = ['Salmon/1.10.1-GCC-12.2.0', 'SAMtools/1.17-GCC-12.2.0']
    container "quay.io/biocontainers/salmon:1.10.1--h7e5ed60_1" // salmon 1.10.1

    tag "Salmon on ${out_prefix}"

    publishDir "${params.outdir}/Salmon/${params.genome}/", mode: 'copy'

    input:
    tuple val(meta), val(out_prefix), path(bam)
    path (tx_fa)
    
    output:
    tuple val(meta), val(out_prefix), path( "${out_prefix}/" ), emit: sf
    path( "${out_prefix}/" )

    script:

    """
    salmon.sh ${task.cpus} ${bam} ${tx_fa} ${out_prefix}
    """
}
