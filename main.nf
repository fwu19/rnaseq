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
    if (params.workflow == "regular" ) {
        RNASEQ_REGULAR()
    }
    else if (params.workflow == "pdx" ) {
        RNASEQ_PDX()
    }
    else if (params.workflow == "exome" ) {
        RNASEQ_REGULAR()
    }
    else {
        err "Invalid workflow! Possible options: [regular, pdx, exome]"
    }
}


/*
* call the main workflow
*/
workflow {
    RNASEQ()
}

