process GENERATE_REPORT {
    module = ['fhR/4.1.2-foss-2021b']

    label "process_single"

    tag "Make plots of read and peak metrics "

    input:
    path( samplesheet, stageAs: "sample_sheet.csv" )
    //path( read_metrics )
    path( "*" )
    path( "*" )

    output:
    tuple path( "*.{rds,html,Rmd}" )

    script:
    """
    prepare_report.r
    render_report.r _Analysis_report.Rmd

    """
}
