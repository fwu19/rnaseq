#!/usr/bin/env nextflow 
// https://github.com/nextflow-io/nextflow/discussions/2892

import groovy.json.JsonOutput

process OUTPUT_PARAMS {

    publishDir "${params.outdir}/pipeline_info/", pattern: '*', mode: 'copy'
    
    input:
        path (out.dir)
        
    output:
        path ('params.json')

    script:
    """
        "echo '${JsonOutput.toJson(params)}' > params.json"
    """
}

