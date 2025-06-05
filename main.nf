#!/usr/bin/env nextflow 

nextflow.enable.dsl = 2

/*
* include workflows
*/
include { RNASEQ } from './workflows/rnaseq'
//include { RNASEQ_PDX } from './workflows/rnaseq_pdx'

/*
* named workflows
*/
workflow RUN_RNASEQ {  
    /*
    * Call named workflow
    */
    if (params.workflow == "regular" ) {
        RNASEQ()
    }else if (params.workflow == "pdx" ) {
        RNASEQ()
    }else if (params.workflow == "exome" ) {
        RNASEQ()
    }else if (params.workflow == "pdx_exome" ) {
        RNASEQ()
    }else {
        err "Invalid workflow! Possible options: [regular, pdx, exome, pdx_exome]"
    }
}


/*
* call the main workflow
*/
workflow {
    RUN_RNASEQ()
}

