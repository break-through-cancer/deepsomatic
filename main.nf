#!/usr/bin/env nextflow

/*
 * DeepSomatic Somatic Variant Calling Pipeline
 *
 * GPU-accelerated somatic variant calling using NVIDIA Parabricks DeepSomatic.
 * Input: Paired tumor-normal BAMs (preprocessed with BQSR applied)
 * Output: VCF files with somatic variant calls
 */

nextflow.enable.dsl = 2

// ---------------------------------------------------------------------------
// PARAMS
// ---------------------------------------------------------------------------

if (!params.containsKey('outdir'))              params.outdir = 'results'
if (!params.containsKey('ref_fasta'))           params.ref_fasta = null
if (!params.containsKey('ref_fai'))             params.ref_fai = null
if (!params.containsKey('deepsomatic_runs'))    params.deepsomatic_runs = null

// DeepSomatic-specific options
if (!params.containsKey('num_streams_per_gpu')) params.num_streams_per_gpu = 3
if (!params.containsKey('enable_gvcf'))         params.enable_gvcf = false

// ---------------------------------------------------------------------------
// PROCESSES
// ---------------------------------------------------------------------------

process DEEPSOMATIC {
    tag "${sampleName}"
    label 'process_high'
    container "nvcr.io/nvidia/clara/clara-parabricks:4.6.0"
    publishDir "${params.outdir}/${sampleName}", mode: 'copy'
    errorStrategy 'retry'
    maxRetries 2

    input:
    tuple val(sampleName), path(tumor_bam), path(tumor_bai), path(normal_bam), path(normal_bai)
    path ref_fasta
    path ref_fai

    output:
    tuple val(sampleName), path("${sampleName}.vcf.gz"), path("${sampleName}.vcf.gz.tbi"), emit: vcf
    path "${sampleName}.log", emit: log

    shell:
    '''
    set -euxo pipefail

    # Run DeepSomatic
    parabricks deepsomatic \
        --in-tumor-bam !{tumor_bam} \
        --in-normal-bam !{normal_bam} \
        --ref !{ref_fasta} \
        --out-variants !{sampleName}.vcf.gz \
        --num-streams-per-gpu !{params.num_streams_per_gpu} \
        !{params.enable_gvcf ? '--run-deepsomatic-gvcf' : ''} \
        > !{sampleName}.log 2>&1

    # Verify output
    test -s !{sampleName}.vcf.gz
    test -s !{sampleName}.vcf.gz.tbi
    '''
}

// ---------------------------------------------------------------------------
// WORKFLOW
// ---------------------------------------------------------------------------

workflow {

    if (!params.deepsomatic_runs) {
        error "params.deepsomatic_runs is empty -- provide tumor/normal BAM pairs"
    }
    if (!params.ref_fasta || !params.ref_fai) {
        error "Missing required reference params: ref_fasta, ref_fai"
    }

    ref_fasta = file(params.ref_fasta, checkIfExists: true)
    ref_fai = file(params.ref_fai, checkIfExists: true)

    runs_ch = Channel
        .fromList(params.deepsomatic_runs)
        .map { run ->
            tuple(
                run.sample_name,
                file(run.tumor_bam, checkIfExists: true),
                file(run.tumor_bam_index, checkIfExists: true),
                file(run.normal_bam, checkIfExists: true),
                file(run.normal_bam_index, checkIfExists: true)
            )
        }

    DEEPSOMATIC(
        runs_ch,
        ref_fasta,
        ref_fai
    )

    // Collect VCFs
    vcf_files = DEEPSOMATIC.out.vcf
        .map { sampleName, vcf, tbi -> vcf }
        .collect()

}
