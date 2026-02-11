#!/usr/bin/env bash

# rename_output.sh ${params.workflow} $srcdir ${params.star_dir} ${params.star_host_dir} ${params.star_log_dir} ${params.star_log_host_dir}
# ${params.workflow} ${params.genome} ${params.genme_host} $srcdir $outdir ${params.star_dir} ${params.star_host_dir} ${params.star_xeno_dir} ${params.bwa_dir} ${params.bwa_host_dir} ${params.bwa_xeno_dir} ${params.featurecounts_dir} ${params.salmon_dir} ${params.arriba_dir} ${params.fastqc_dir} ${params.fastqc_trimmed_dir} ${params.cutadapt_dir} ${params.fastp_dir} ${params.star_log_dir} ${params.star_log_host_dir} ${params.rnaseqc_dir} ${params.rseqc_dir} ${params.gatk_dir} ${params.samtools_dir} ${params.samtools_host_dir} ${params.samtools_xeno_dir} ${params.multiqc_dir} ${params.gene_expr_dir} ${params.tx_expr_dir} ${params.de_dir} ${params.dt_dir} 

workflow=$1; shift
#genome=$1; shift
#genome_host=$1; shift
srcdir=$1; shift
outdir=$1; shift
star_dir=$1; shift
star_host_dir=$1; shift
star_xeno_dir=$1; shift
bwa_dir=$1; shift
bwa_host_dir=$1; shift
bwa_xeno_dir=$1; shift
featurecounts_dir=$1; shift
salmon_dir=$1; shift
arriba_dir=$1; shift
fastqc_dir=$1; shift
fastqc_trimmed_dir=$1; shift
cutadapt_dir=$1; shift
fastp_dir=$1; shift
star_log_dir=$1; shift
star_log_host_dir=$1; shift
rnaseqc_dir=$1; shift
rseqc_dir=$1; shift
gatk_dir=$1; shift
samtools_dir=$1; shift
samtools_host_dir=$1; shift
samtools_xeno_dir=$1; shift
multiqc_dir=$1; shift
gene_expr_dir=$1; shift
tx_expr_dir=$1; shift
de_dir=$1; shift
dt_dir=$1; shift

if [[ $srcdir != $outdir ]]; then

subdir=${star_xeno_dir}
if [[ -d ${srcdir}/${subdir} ]]; then
    mkdir -p ${outdir}/${subdir}
    ln -sf ${srcdir}/${subdir}/* ${outdir}/${subdir}/
fi

subdir=${bwa_dir}
if [[ -d ${srcdir}/${subdir} ]]; then
    mkdir -p ${outdir}/${subdir}
    ln -sf ${srcdir}/${subdir}/* ${outdir}/${subdir}/
fi

subdir=${bwa_host_dir}
if [[ -d ${srcdir}/${subdir} ]]; then
    mkdir -p ${outdir}/${subdir}
    ln -sf ${srcdir}/${subdir}/* ${outdir}/${subdir}/
fi

subdir=${bwa_xeno_dir}
if [[ -d ${srcdir}/${subdir} ]]; then
    mkdir -p ${outdir}/${subdir}
    ln -sf ${srcdir}/${subdir}/* ${outdir}/${subdir}/
fi

subdir=${featurecounts_dir}
if [[ -d ${srcdir}/${subdir} ]]; then
    mkdir -p ${outdir}/${subdir}
    ln -sf ${srcdir}/${subdir}/* ${outdir}/${subdir}/
fi

subdir=${salmon_dir}
if [[ -d ${srcdir}/${subdir} ]]; then
    mkdir -p ${outdir}/${subdir}
    ln -sf ${srcdir}/${subdir}/* ${outdir}/${subdir}/
fi

subdir=${arriba_dir}
if [[ -d ${srcdir}/${subdir} ]]; then
    mkdir -p ${outdir}/${subdir}
    ln -sf ${srcdir}/${subdir}/* ${outdir}/${subdir}/
fi

subdir=${fastqc_dir}
if [[ -d ${srcdir}/${subdir} ]]; then
    mkdir -p ${outdir}/${subdir}
    ln -sf ${srcdir}/${subdir}/* ${outdir}/${subdir}/
fi

subdir=${fastqc_trimmed_dir}
if [[ -d ${srcdir}/${subdir} ]]; then
    mkdir -p ${outdir}/${subdir}
    ln -sf ${srcdir}/${subdir}/* ${outdir}/${subdir}/
fi

subdir=${cutadapt_dir}
if [[ -d ${srcdir}/${subdir} ]]; then
    mkdir -p ${outdir}/${subdir}
    ln -sf ${srcdir}/${subdir}/* ${outdir}/${subdir}/
fi

subdir=${fastp_dir}
if [[ -d ${srcdir}/${subdir} ]]; then
    mkdir -p ${outdir}/${subdir}
    ln -sf ${srcdir}/${subdir}/* ${outdir}/${subdir}/
fi

subdir=${rnaseqc_dir}
if [[ -d ${srcdir}/${subdir} ]]; then
    mkdir -p ${outdir}/${subdir}
    ln -sf ${srcdir}/${subdir}/* ${outdir}/${subdir}/
fi

subdir=${rseqc_dir}
if [[ -d ${srcdir}/${subdir} ]]; then
    mkdir -p ${outdir}/${subdir}
    ln -sf ${srcdir}/${subdir}/* ${outdir}/${subdir}/
fi

subdir=${gatk_dir}
if [[ -d ${srcdir}/${subdir} ]]; then
    mkdir -p ${outdir}/${subdir}
    ln -sf ${srcdir}/${subdir}/* ${outdir}/${subdir}/
fi

subdir=${samtools_dir}
if [[ -d ${srcdir}/${subdir} ]]; then
    mkdir -p ${outdir}/${subdir}
    ln -sf ${srcdir}/${subdir}/* ${outdir}/${subdir}/
fi

subdir=${samtools_host_dir}
if [[ -d ${srcdir}/${subdir} ]]; then
    mkdir -p ${outdir}/${subdir}
    ln -sf ${srcdir}/${subdir}/* ${outdir}/${subdir}/
fi

subdir=${samtools_xeno_dir}
if [[ -d ${srcdir}/${subdir} ]]; then
    mkdir -p ${outdir}/${subdir}
    ln -sf ${srcdir}/${subdir}/* ${outdir}/${subdir}/
fi

subdir=${multiqc_dir}
if [[ -d ${srcdir}/${subdir} ]]; then
    mkdir -p ${outdir}/${subdir}
    ln -sf ${srcdir}/${subdir}/* ${outdir}/${subdir}/
fi

subdir=${gene_expr_dir}
if [[ -d ${srcdir}/${subdir} ]]; then
    mkdir -p ${outdir}/${subdir}
    ln -sf ${srcdir}/${subdir}/* ${outdir}/${subdir}/
fi

subdir=${tx_expr_dir}
if [[ -d ${srcdir}/${subdir} ]]; then
    mkdir -p ${outdir}/${subdir}
    ln -sf ${srcdir}/${subdir}/* ${outdir}/${subdir}/
fi

subdir=${de_dir}
if [[ -d ${srcdir}/${subdir} ]]; then
    mkdir -p ${outdir}/${subdir}
    ln -sf ${srcdir}/${subdir}/* ${outdir}/${subdir}/
fi

subdir=${dt_dir}
if [[ -d ${srcdir}/${subdir} ]]; then
    mkdir -p ${outdir}/${subdir}
    ln -sf ${srcdir}/${subdir}/* ${outdir}/${subdir}/
fi

fi


## extra work for STAR output, to rename some files and organize them in subdirectories
if [[ -d ${srcdir}/${star_dir} ]]; then
    mkdir -p ${outdir}/${star_dir}/_work
    ln -sf ${srcdir}/${star_dir}/*.bam ${outdir}/${star_dir}/
    ln -sf ${srcdir}/${star_dir}/*.bam.bai ${outdir}/${star_dir}/

    for i in $(ls ${srcdir}/${star_dir}/_work/); do
        mkdir -p ${outdir}/${star_dir}/_work/${i}
        ln -sf ${srcdir}/${star_dir}/_work/${i}/* ${outdir}/${star_dir}/_work/${i}/

        srcfile="${srcdir}/${star_dir}/_work/${i}/ReadsPerGene.out.tab"
        dstfile="${outdir}/${star_dir}/_work/${i}/${i}.ReadsPerGene.out.tab"
        [[ -f "${srcfile}"  ]] && ln -sf "${srcfile}" "${dstfile}"

        srcfile="${srcdir}/${star_dir}/_work/${i}/Aligned.toTranscriptome.out.bam "
        dstfile="${outdir}/${star_dir}/_work/${i}/${i}.Aligned.toTranscriptome.out.bam "
        [[ -f "${srcfile}"  ]] && ln -sf "${srcfile}" "${dstfile}"

        srcfile="${srcdir}/${star_dir}/_work/${i}/Log.final.out"
        dstfile="${outdir}/${star_dir}/_work/${i}/${i}.Log.final.out"
        [[ -f "${srcfile}"  ]] && ln -sf "${srcfile}" "${dstfile}"

    done

    if [[ -d ${srcdir}/${star_log_dir} ]]; then
        mkdir -p ${outdir}/${star_log_dir}
        ln -sf ${srcdir}/${star_log_dir}/* ${outdir}/${star_log_dir}/
    else
        mkdir -p ${outdir}/${star_log_dir}
        for i in $(ls ${srcdir}/${star_dir}/_work/); do
                srcfile="${outdir}/${star_dir}/_work/${i}/${i}.Log.final.out"
                dstfile="${outdir}/${star_log_dir}/${i}.Log.final.out"
                [[ -f "${srcfile}"  ]] && ln -sf "${srcfile}" "${dstfile}"
        done
    fi
fi

if [[ -d ${srcdir}/${star_host_dir} ]]; then
    mkdir -p ${outdir}/${star_host_dir}/_work
    ln -sf ${srcdir}/${star_host_dir}/*.bam ${outdir}/${star_host_dir}/
    ln -sf ${srcdir}/${star_host_dir}/*.bam.bai ${outdir}/${star_host_dir}/

    for i in $(ls ${srcdir}/${star_host_dir}/_work/); do
        mkdir -p ${outdir}/${star_host_dir}/_work/${i}
        ln -sf ${srcdir}/${star_host_dir}/_work/${i}/* ${outdir}/${star_host_dir}/_work/${i}/

        srcfile="${srcdir}/${star_host_dir}/_work/${i}/ReadsPerGene.out.tab"
        dstfile="${outdir}/${star_host_dir}/_work/${i}/${i}.ReadsPerGene.out.tab"
        [[ -f "${srcfile}"  ]] && ln -sf "${srcfile}" "${dstfile}"

        srcfile="${srcdir}/${star_host_dir}/_work/${i}/Aligned.toTranscriptome.out.bam "
        dstfile="${outdir}/${star_host_dir}/_work/${i}/${i}.Aligned.toTranscriptome.out.bam "
        [[ -f "${srcfile}"  ]] && ln -sf "${srcfile}" "${dstfile}"

        srcfile="${srcdir}/${star_host_dir}/_work/${i}/Log.final.out"
        dstfile="${outdir}/${star_host_dir}/_work/${i}/${i}.Log.final.out"
        [[ -f "${srcfile}"  ]] && ln -sf "${srcfile}" "${dstfile}"
    done

    if [[ -d ${srcdir}/${star_log_host_dir} ]]; then
        mkdir -p ${outdir}/${star_log_host_dir}
        ln -sf ${srcdir}/${star_log_host_dir}/* ${outdir}/${star_log_dir}/
    else
        mkdir -p ${outdir}/${star_log_host_dir}
        for i in $(ls ${srcdir}/${star_host_dir}/_work/); do
                srcfile="${outdir}/${star_host_dir}/_work/${i}/${i}.Log.final.out"
                dstfile="${outdir}/${star_log_host_dir}/${i}.Log.final.out"
                [[ -f "${srcfile}"  ]] && ln -sf "${srcfile}" "${dstfile}"
        done
    fi
fi

## make links for analysis report
for i in $(find ${srcdir} -maxdepth 1 -type f  | egrep "*.html"); do
    [[ -f $i ]] && ln -sf $i ${outdir}/
done

