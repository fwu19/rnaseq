nextflow.enable.dsl=2

/*
* Writes out csv files for v0.3
*/

include { MAKE_SYMLINKS } from '../modules/make_symlinks.nf'
include { RENAME_OUTPUT } from '../modules/rename_output.nf'
include { WRITE_CSV as WRITE_CSV_ALIGN_FASTQ } from '../modules/write_csv.nf'
include { WRITE_CSV as WRITE_CSV_FC_COUNTS } from '../modules/write_csv.nf'
include { WRITE_CSV as WRITE_CSV_STAR_COUNTS } from '../modules/write_csv.nf'
include { WRITE_CSV as WRITE_CSV_MAP_TX } from '../modules/write_csv.nf'
include { WRITE_CSV as WRITE_CSV_SALMON } from '../modules/write_csv.nf'
include { WRITE_CSV as WRITE_CSV_ARRIBA } from '../modules/write_csv.nf'

workflow UPDATE_OUTPUT {
    take:
    samplesheet
    srcdir
    outdir

    main: 

    
    RENAME_OUTPUT(
        srcdir,
        outdir
    )

    samplesheet
        .splitCsv( header: true )
        .map { row -> [ id: row.id, sample_group: row.sample_group ] }
        .set { ch_metadata }

    /* regular or exome workflow, map2genome.csv */
    if (params.workflow != 'pdx' ){
        star_dir = file("${srcdir}/${params.star_dir}/")
        if (star_dir.exists()){
            WRITE_CSV_ALIGN_FASTQ(
                ch_metadata
                .map { 
                    it -> [id: it.id] + [sample_group: it.sample_group] + [bam: "${params.star_dir}/${it.id}.bam"] + [bai: "${params.star_dir}/${it.id}.bam.bai"]
                }
                .collect(),
                "map2genome.star.csv"        
            )
        
            /* map2transcriptome.star.csv */
            WRITE_CSV_MAP_TX(
                ch_metadata
                .map { 
                    it -> [id: it.id] + [sample_group: it.sample_group] + [tx_bam: "${params.star_dir}/_work/${it.id}/${it.id}.Aligned.toTranscriptome.out.bam"]
                }
                .collect(),
                "map2transcriptome.star.csv"        
            )

            /* gene_counts.star.csv */
            WRITE_CSV_STAR_COUNTS(
                ch_metadata
                .map { 
                    it -> [id: it.id] + [sample_group: it.sample_group] + [gene_count: "${params.star_dir}/_work/${it.id}/${it.id}.ReadsPerGene.out.tab" ] 
                }
                .collect(),
            "gene_counts.star.csv"        
            )

        }

        bwa_dir = file("${srcdir}/${params.bwa_dir}/")
        if (bwa_dir.exists()){
            WRITE_CSV_ALIGN_FASTQ(
                ch_metadata
                .map { 
                    it -> [id: it.id] + [sample_group: it.sample_group] + [bam: "${params.star_dir}/${it.id}.bam"] + [bai: "${params.bwa_dir}/${it.id}.bam.bai"]
                }
                .collect(),
                "map2genome.bwa-mem.csv"        
            )
        
        }
    }

    /* pdx workflow, map2genome.csv */
    if(params.workflow == 'pdx'){
        star_dir = file("${srcdir}/${params.star_dir}")
        if(star_dir.exists()){
            WRITE_CSV_ALIGN_FASTQ(
                ch_metadata
                .map { 
                    it -> [id: it.id] + [sample_group: it.sample_group] + [graft_bam: "${params.star_dir}/${it.id}.bam" ] + [graft_bai: "${params.star_dir}/${it.id}.bam.bai" ] + [host_bam: "${params.star_host_dir}/${it.id}.bam" ] + [host_bai: "${params.star_host_dir}/${it.id}.bam.bai" ] + [bam: "${params.star_xeno_dir}/${it.id}.bam" ] + [bai: "${params.star_xeno_dir}/${it.id}.bam.bai" ]
                    }
                .collect(),
            "map2genome.star.csv"        
            )
        }     

        bwa_dir = file("${srcdir}/${params.bwa_dir}")
        if(bwa_dir.exists()){
            WRITE_CSV_ALIGN_FASTQ(
                ch_metadata
                .map { 
                    it -> [id: it.id] + [sample_group: it.sample_group] + [graft_bam: "${params.bwa_dir}/${it.id}.bam" ] + [graft_bai: "${params.bwa_dir}/${it.id}.bam.bai" ] + [host_bam: "${params.bwa_host_dir}/${it.id}.bam" ] + [host_bai: "${params.bwa_host_dir}/${it.id}.bam.bai" ] + [bam: "${params.bwa_xeno_dir}/${it.id}.bam" ] + [bai: "${params.bwa_xeno_dir}/${it.id}.bam.bai" ]
                    }
                .collect(),
            "map2genome.bwa-mem.csv"        
            )
        }     
    }

    /* gene_counts.featurecounts.csv */
    fc_dir = file("${srcdir}/${params.featurecounts_dir}")
    if (fc_dir.exists()){
        WRITE_CSV_FC_COUNTS(
            ch_metadata
                .map { 
                    it -> [id: it.id] + [sample_group: it.sample_group] + [gene_count: "${params.featurecounts_dir}/${it.id}.exonicFragmentCounts.txt" ] 
                }
                .collect(),
            "gene_counts.featurecounts.csv"        
        )

    }

    /* salmon.csv */
    salmon_dir = file("${srcdir}/${params.salmon_dir}")  
    if (salmon_dir.exists()){
        WRITE_CSV_SALMON(
            ch_metadata
                .map { 
                    it -> [id: it.id] + [sample_group: it.sample_group] + [salmon: "${params.salmon_dir}/${it.id}/"]
                }
                .collect(),
        "salmon.csv"        
        )
    }

    /* arriba.csv */
    arriba_dir = file("${srcdir}/${params.arriba_dir}")
    if (arriba_dir.exists()){
        WRITE_CSV_ARRIBA(
            ch_metadata
                .map { 
                    it -> [id: it.id] + [sample_group: it.sample_group] + [arriba: "${params.arriba_dir}/${it.id}.fusions.tsv"]
                }
                .collect(),
        "arriba.csv"        
        )
    }

    
}
