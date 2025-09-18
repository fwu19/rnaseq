
process DIFFERENTIAL_EXPRESSION {
    time = '1d'
    cpus = 1
    memory = '24G'
    container "docker://fwu19/r-libs:4.1.2" 


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
    path(gene_txt)
    val(gene_type)

    output:
    path ("differential_genes.rds"), emit: rds, optional: true
    path ("*"), optional: true
    
    script:
    """
    differential_expression.r input=$samplesheet comparison=$comparison rds=${rds} fdr=${fdr} fc=${fc} fdr2=${fdr2} fc2=${fc2} gene_txt=$gene_txt gene_type=$gene_type
    """
}
