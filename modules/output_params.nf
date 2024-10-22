#!/usr/bin/env nextflow 
// https://github.com/nextflow-io/nextflow/discussions/2892

import groovy.json.JsonOutput

process OUTPUT_PARAMS {
    label "process_single"
    
    input:
        path (outdir)

    output:
        path ('params.json')

    script:
    """
        echo '${JsonOutput.toJson(params)}' > params.json
    """
}

