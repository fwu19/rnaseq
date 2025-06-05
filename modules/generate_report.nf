process GENERATE_REPORT {
    module = ['fhR/4.1.2-foss-2021b']

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
    path ( "*" )

    output:
    path( "*.{rds,html,Rmd}" )

    script:
    """
    prepare_report.r $workflow
    render_report.r *.Rmd

    """
}
