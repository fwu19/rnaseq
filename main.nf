#!/usr/bin/env nextflow 

/*
* include workflows
*/
include { RNASEQ_REGULAR } from './workflows/rnaseq_regular'
include { RNASEQ_PDX } from './workflows/rnaseq_pdx'

/*
* named workflows
*/
workflow RNASEQ {  
    /*
    * Call named workflow
    */
    if (params.workflow_regular) {
        RNASEQ_REGULAR()
    }
    else if (params.workflow_pdx) {
        RNASEQ_PDX()
    }
    else {
        err "Invalid workflow! Possible options: [workflow_regular, workflow_pdx, workflow_exome]"
    }
}


/*
* call the main workflow
*/
workflow {
    RNASEQ()
}

