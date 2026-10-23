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
if (!params.containsKey('num_streams_per_gpu')) params.num_streams_per_gpu = 1
if (!params.containsKey('enable_gvcf'))         params.enable_gvcf = false
if (!params.containsKey('interval_file'))       params.interval_file = null   // optional BED to restrict calling

// Empty list stands in for "no file" so the optional path input stages nothing
def NO_FILE = []

// ---------------------------------------------------------------------------
// PROCESSES
// ---------------------------------------------------------------------------

process DEEPSOMATIC {
    tag "${sampleName}"
    label 'process_high'
    container "nvcr.io/nvidia/clara/clara-parabricks:4.6.0-1"
    // Published files carry the caller name (<sample>.deepsomatic.vcf.gz) so downstream
    // tools such as consensus-calling can tell callers apart by filename. saveAs only
    // changes the published copy; it is not part of the task cache key.
    publishDir "${params.outdir}/${sampleName}", mode: 'copy',
        saveAs: { fn -> fn.replace("${sampleName}.", "${sampleName}.deepsomatic.") }
    errorStrategy 'retry'
    maxRetries 2

    input:
    tuple val(sampleName), path(tumor_bam), path(tumor_bai), path(normal_bam), path(normal_bai)
    path ref_fasta
    path ref_fai
    path interval_file

    output:
    tuple val(sampleName), path("${sampleName}.vcf.gz"), emit: vcf
    path "${sampleName}.vcf.gz.tbi", optional: true, emit: tbi
    path "${sampleName}.log", emit: log

    shell:
    '''
    set -euxo pipefail

    # Optionally restrict calling to the regions in a BED file (genome-wide when unset)
    interval_arg=""
    if [ -n "!{interval_file}" ]; then
        interval_arg="--interval-file !{interval_file}"
    fi

    # Run DeepSomatic
    pbrun deepsomatic \
        --in-tumor-bam !{tumor_bam} \
        --in-normal-bam !{normal_bam} \
        --ref !{ref_fasta} \
        --out-variants !{sampleName}.vcf.gz \
        --num-streams-per-gpu !{params.num_streams_per_gpu} \
        !{params.enable_gvcf ? '--run-deepsomatic-gvcf' : ''} \
        ${interval_arg} \
        > !{sampleName}.log 2>&1

    # Verify output. pbrun may not write a .tbi, so index it if we can, but don't fail the
    # task over a missing index (the VCF itself can take hours to produce).
    test -s !{sampleName}.vcf.gz
    if [ ! -s !{sampleName}.vcf.gz.tbi ] && command -v tabix >/dev/null 2>&1; then
        tabix -p vcf !{sampleName}.vcf.gz || true
    fi
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

    interval_file = params.interval_file ? file(params.interval_file, checkIfExists: true) : NO_FILE

    DEEPSOMATIC(
        runs_ch,
        ref_fasta,
        ref_fai,
        interval_file
    )

    // Collect VCFs
    vcf_files = DEEPSOMATIC.out.vcf
        .map { sampleName, vcf -> vcf }
        .collect()

}
