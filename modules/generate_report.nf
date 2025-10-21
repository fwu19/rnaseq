process GENERATE_REPORT {
    container "docker://fwu19/r-libs:4.1.2" 

    label "process_single"

    tag "Make plots of read and peak metrics "

    input:
    val( workflow )
    path ( samplesheet, stageAs: "sample_sheet.csv" )
    path ( "multiqc_data/" )
    path ( gene_rds )
    path ( de_rds )
    path ( txt_rds )
    path ( dt_rds )
    path ( "hs_metrics/" )
    path ( report_dir, stageAs: "report" )

    output:
    path( "*.{rds,html,Rmd}" )

    script:
    """
    generate_report.r $workflow
    
    """
}
