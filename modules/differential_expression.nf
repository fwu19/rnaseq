
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
    path("*", stageAs: "counts/*")
    path(gene_txt)
    val(length_col)
    val(strand)
    val(fdr)
    val(fc)
    val(fdr2)
    val(fc2)

    output:
    path ("*.rds"), emit: rds
    path ("*")
    
    script:
    """
    differential_expression.r input=$samplesheet comparison=$comparison count.dir=counts gene.txt=$gene_txt length.col=${length_col} strand=${strand} fdr=${fdr} fc=${fc} fdr2=${fdr2} fc2=${fc2}
    rm -r counts/
    """
}
