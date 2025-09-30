process GENERATE_REPORT {
    container "docker://fwu19/r-libs:4.1.2" 

    label "process_single"

    tag "Generate an analysis report "

    input:
    val( workflow )
    path ( samplesheet, stageAs: "sample_sheet.csv" )
    path ( "multiqc_data/" )
    path ( gene_rds )
    path ( de_rds )
    path ( txt_rds )
    path ( dt_rds )
    path ( "hs_metrics/" )
    path ( report_dir )

    output:
    path( "*.{rds,html,Rmd}" )

    script:
    """
    prepare_report.r $workflow ${report_dir}
    render_report.r *.Rmd

    """
}
