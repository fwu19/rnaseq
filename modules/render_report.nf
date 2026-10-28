process RENDER_REPORT {
    container "docker://fwu19/r-libs:4.1.2" 

    label "process_single"

    tag "Render a report"

    input:
    path ( "sample_sheet.csv" )
    path ( "multiqc_data" )
    path ( "gene_counts/*" )
    path ( "gene_counts_host/*" )
    path ( "de/*" )
    path ( "de_host/*" )
    path ( "tx_counts/*" )
    path ( "tx_counts_host/*" )
    path ( "dt/*" )
    path ( "dt_host/*" )
    path ( "hs_metrics/" )
    path ( "arriba/" )
    path ( "_report.Rmd" )

    output:
    path( "*.{Rmd,rmd}" ), emit: rmd
    path( "*.{rds,html,Rmd,rmd}" )
    path ('versions.yml'), emit: versions

    script:
    def args = task.ext.args ?: ''
    """
    render_report.r report_rmd="_report.Rmd" report_html=${params.report_html} workflow=${params.workflow} genome=${params.genome} genome_host=${params.genome_host} pdx_tool=${params.pdx_tool} fdr=${params.fdr} fc=${params.fc} fdr2=${params.fdr2} fc2=${params.fc2} $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: \$(R --version 2>&1 | head -n 1 | sed -e "s/R version //g; s/ .*//g" )
    END_VERSIONS
    
    """
}
