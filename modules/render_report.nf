process RENDER_REPORT {
    container "docker://fwu19/r-libs:4.1.2" 

    label "process_single"

    tag "Render an analysis report "

    input:
    path ( "*" )

    output:
    path( "*.{html}" )

    script:
    """
    render_report.r *.Rmd

    """
}
