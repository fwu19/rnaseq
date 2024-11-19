#!/usr/bin/env nextflow 

/*
* include workflows
*/
include { RNASEQ_REGULAR } from './workflows/rnaseq_regular'
include { RNASEQ_PDX } from './workflows/rnaseq_pdx'
include { RNASEQ_EXOME } from './workflows/rnaseq_exome'

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
        RNASEQ_EXOME()
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

