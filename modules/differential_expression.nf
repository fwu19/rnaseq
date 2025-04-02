
process DIFFERENTIAL_EXPRESSION {
    time = '1d'
    cpus = 1
    memory = '24G'
    module = ['fhR/4.1.2-foss-2021b']


    tag "Differential expression"

    publishDir "${params.outdir}/differential_expression/", mode: 'copy'
    
    input:
    path(samplesheet)
    path(comparison)
    path(rds)
    val(fdr)
    val(fc)
    val(fdr2)
    val(fc2)

    output:
    path ("differential_genes.rds"), emit: rds
    path ("*")
    
    script:
    """
    differential_expression.r input=$samplesheet comparison=$comparison rds=${rds} fdr=${fdr} fc=${fc} fdr2=${fdr2} fc2=${fc2}
    """
}
