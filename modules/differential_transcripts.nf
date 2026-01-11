
process DIFFERENTIAL_TRANSCRIPTS {
    label "process_single"
    container "docker://fwu19/r-libs:4.1.2" 


    tag "Differential transcripts"

    publishDir "${params.outdir}/differential_transcripts/", mode: 'copy'
    
    input:
    path(samplesheet)
    path(comparison)
    path(count_file)
    val(length_col)
    path(gene_txt)

    output:
    path ("differential_transcripts.rds"), emit: rds, optional: true
    path ('versions.yml'), emit: versions
    path ("*"), optional: true
    
    script:
    """
    differential_expression.r input="\'$samplesheet\'" comparison="\'$comparison\'" gene_txt="\'${gene_txt}\'" count_file="\'${count_file}\'" fdr=${params.fdr} fc=${params.fc} fdr2=${params.fdr2} fc2=${params.fc2} gene_type="\'${params.de_gene_type}\'" length_col="\'$length_col\'"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: \$(R --version | head -n 1)
    END_VERSIONS

    """
}
