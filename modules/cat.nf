
process CAT {
    label 'process_single'

    tag "CAT_FASTQ on $id"

    //publishDir "${params.outdir}/merged_fastq/", pattern: '*.fastq.gz', mode: 'copy'

    input:
    path( "input*/*" )
    val( out_file )
    val( type )
    
    output:
    path( out_file, emit: file)

    script:
    """
    #!/bin/bash
    
    if [[ $type == "gtf" ]]; then
        cat input*/* | egrep -v "^#" > $out_file
    else
        cat input*/* >$out_file
    fi


    """

}
