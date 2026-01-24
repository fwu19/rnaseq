process RENDER_REPORT {
    container "docker://fwu19/r-libs:4.1.2" 

    label "process_single"

    tag "Render a report"

    input:
    path ( "sample_sheet.csv" )
    path ( "multiqc_data" )
    path ( gene_counts )
    path ( de )
    path ( txt_counts )
    path ( dt )
    path ( "hs_metrics/" )
    path ( "report" )

    output:
    path( "*.{rds,html,Rmd}" )
    path ('versions.yml'), emit: versions

    script:
    """
    render_report.r workflow=${params.workflow} fdr=${params.fdr} fc=${params.fc} fdr2=${params.fdr2} fc2=${params.fc2}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: \$(R --version 2>&1 | head -n 1 | sed -e "s/R version //g; s/ .*//g" )
    END_VERSIONS
    
    """
}
