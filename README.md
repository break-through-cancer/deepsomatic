# DeepSomatic Pipeline

GPU-accelerated somatic variant calling using NVIDIA Parabricks DeepSomatic.

## Requirements

- **GPU**: Minimum 40GB VRAM (H100, A100, A10, or T4)
- **Input**: Paired tumor-normal BAMs (preprocessed with BQSR applied)
- **Output**: VCF files with somatic variant calls

## Input Format

Provide tumor-normal BAM pairs via JSON config:

```json
{
  "deepsomatic_runs": [
    {
      "sample_name": "SAMPLE1",
      "tumor_bam": "s3://bucket/SAMPLE1.tumor.bam",
      "tumor_bam_index": "s3://bucket/SAMPLE1.tumor.bam.bai",
      "normal_bam": "s3://bucket/SAMPLE1.normal.bam",
      "normal_bam_index": "s3://bucket/SAMPLE1.normal.bam.bai"
    }
  ],
  "ref_fasta": "s3://broad-references/hg38/v0/Homo_sapiens_assembly38.fasta",
  "ref_fai": "s3://broad-references/hg38/v0/Homo_sapiens_assembly38.fasta.fai"
}
```

## Parameters

- `num_streams_per_gpu`: Number of parallel streams per GPU (1-6). Default: 3
- `enable_gvcf`: Generate gVCF output. Default: false

## Cost Estimation

**Per-sample (WGS 30x):**
- T4: ~$2-3 (4-6 hours)
- A10: ~$3-5 (2-3 hours)
- H100 Spot: ~$1-2 (30 min, cheapest)

## Usage

```bash
nextflow run main.nf -params-file config.json
nextflow run main.nf -params-file config.json --num_streams_per_gpu 5
nextflow run main.nf -params-file config.json --enable_gvcf true
```

## Notes

- BAMs must be preprocessed with BQSR
- GPU automatically parallelizes within sample (no interval scattering)
- Edit `.cirro/process-compute.config` to select GPU type
