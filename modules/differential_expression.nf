
process DIFFERENTIAL_EXPRESSION {
    time = '1d'
    cpus = 1
    memory = '24G'
    module = ['fhR/4.1.2-foss-2021b']


    tag "Differential expression"

    publishDir "${params.outdir}/differential_expression/", mode: 'copy'
    
    input:
    path ("*.txt", stageAs: "counts/")

    output:
    tuple path ("*")
    
    script:
    """
    Rscript differential_expression.r ss=${params.input} comparison=${params.comparison} count.dir=counts gene.txt=${params.geneTxt} length.col=${params.lengthCol} strand=${params.strand} fdr=${params.fdr} fc=${params.fc} fdr2=${params.fdr2} fc2=${params.fc2}
    """
}
