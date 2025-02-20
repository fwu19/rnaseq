
process ARRIBA {
    time = '1d'
    cpus = 8
    memory = '48G'
    module = []


    tag "Arriba on ${out_prefix}"

    publishDir "${params.outdir}/Arriba/${params.genome}/", mode: 'copy'

    input:
    tuple val(meta), val(out_prefix), path(bam)
    
    output:
    tuple val(meta), val(out_prefix), path( "${out_prefix}.*" ), emit: out 
    path( "${out_prefix}.*" )

    script:
    def args = task.ext.args ?: ""
    """
    arriba.sh ${bam} ${out_prefix}
    
    """
}
