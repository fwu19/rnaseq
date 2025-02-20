#!/usr/bin/env bash

# get arguments
BASE_DIR="/fh/fast/_SR/Genomics/proj/dstirewa_NPM1/arriba_v2.4.0/"
# STAR_INDEX_DIR="${BASE_DIR}/STAR_index_hg38_GENCODE38"
ANNOTATION_GTF="${BASE_DIR}/GENCODE38.gtf"
ASSEMBLY_FA="${BASE_DIR}/hg38.fa"
BLACKLIST_TSV="${BASE_DIR}/database/blacklist_hg38_GRCh38_v2.4.0.tsv.gz"
KNOWN_FUSIONS_TSV="${BASE_DIR}/database/known_fusions_hg38_GRCh38_v2.4.0.tsv.gz"
TAGS_TSV="$KNOWN_FUSIONS_TSV" # different files can be used for filtering and tagging, but the provided one can be used for both
PROTEIN_DOMAINS_GFF3="${BASE_DIR}/database/protein_domains_hg38_GRCh38_v2.4.0.gff3"

bam=$1; shift
out_prefix=$1; shift


"$BASE_DIR/arriba" \
        -x $bam \
        -o ${out_prefix}.fusions.tsv -O ${out_prefix}.fusions.discarded.tsv \
        -a "$ASSEMBLY_FA" -g "$ANNOTATION_GTF" -b "$BLACKLIST_TSV" -k "$KNOWN_FUSIONS_TSV" -t "$TAGS_TSV" -p "$PROTEIN_DOMAINS_GFF3" | tee ${out_prefix}.arriba.log
