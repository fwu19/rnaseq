process PREPARE_REPORT {
    container "docker://fwu19/r-libs:4.1.2" 

    label "process_single"

    tag "Prepare an analysis report "

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
    path( "*.{rds,Rmd}", emit: report )

    script:
    """
    prepare_report.r $workflow ${report_dir}
    
    """
}
