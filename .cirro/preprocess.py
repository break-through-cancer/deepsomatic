#!/usr/bin/env python3
"""
preprocess.py -- Cirro dataset preprocessing hook for DeepSomatic pipeline.
Runs BEFORE Nextflow starts, inspects the dataset's file manifest, and
auto-discovers tumor-normal BAM pairs, injecting them as `deepsomatic_runs` param.

Shared-normal model: every tumor BAM is paired with the SAME single normal BAM.
Classification: sample names containing "PBMC" are treated as normal;
all others are treated as tumors.
"""
import json
from cirro.helpers.preprocess_dataset import PreprocessDataset

NORMAL_NAME_MARKERS = ["PBMC"]  # adjust if normals use different naming


def extract_bams(ds):
    """Extract BAM/BAI file pairs from dataset manifest, grouped by sample."""
    df = ds.files.copy()
    df["file"] = df["file"].astype(str)

    # Filter to only BAM and BAI files
    df = df[df["file"].str.endswith(".bam") | df["file"].str.endswith(".bam.bai")]

    bam_map = {}

    for sample, group in df.groupby("sample"):
        bam = ""
        bai = ""

        for f in group["file"]:
            if f.endswith(".bam") and not f.endswith(".bam.bai"):
                bam = f
            elif f.endswith(".bam.bai"):
                bai = f

        if bam:
            bam_map[str(sample)] = {"bam": bam, "bai": bai}

    if not bam_map:
        raise ValueError("No BAMs found in dataset")

    return bam_map


def is_normal(sample_name):
    """Check if sample name indicates normal/germline sample."""
    upper = sample_name.upper()
    return any(marker in upper for marker in NORMAL_NAME_MARKERS)


def main():
    """Extract BAM pairs and create deepsomatic_runs param."""
    ds = PreprocessDataset.from_running()

    print("=== Dataset files (preview) ===")
    print(ds.files.head(20).to_string(index=False))

    bam_map = extract_bams(ds)
    print(f"\n=== Found {len(bam_map)} BAM file(s) ===")

    # Auto-discover normal sample
    normal_bam = None
    normal_bai = None
    normal_sample_name = None
    tumor_samples = []

    for sample, files in bam_map.items():
        sample_name = str(sample)

        if is_normal(sample_name):
            if normal_bam is not None:
                raise ValueError(
                    f"Multiple normal-labeled samples found ({normal_sample_name!r} and "
                    f"{sample_name!r}) -- expects exactly one shared normal per dataset. "
                    f"Adjust NORMAL_NAME_MARKERS or dataset contents."
                )
            normal_bam = files["bam"]
            normal_bai = files["bai"]
            normal_sample_name = sample_name
            print(f"  Normal: {sample_name}")
        else:
            tumor_samples.append({
                "sample_name": sample_name,
                "tumor_bam": files["bam"],
                "tumor_bam_index": files["bai"],
            })
            print(f"  Tumor:  {sample_name}")

    if not tumor_samples:
        raise ValueError("No tumor BAM found")

    if normal_bam is None:
        raise ValueError(
            "No normal-labeled sample found -- expects a sample name containing "
            f"one of {NORMAL_NAME_MARKERS} (see NORMAL_NAME_MARKERS)"
        )

    # Create deepsomatic_runs with tumor-normal pairings
    deepsomatic_runs = []
    for tumor in tumor_samples:
        deepsomatic_runs.append({
            "sample_name": tumor["sample_name"],
            "tumor_bam": tumor["tumor_bam"],
            "tumor_bam_index": tumor["tumor_bam_index"],
            "normal_bam": normal_bam,
            "normal_bam_index": normal_bai,
        })

    ds.add_param("deepsomatic_runs", deepsomatic_runs)

    print("\n=== Final parameters ===")
    print(json.dumps(ds.params, indent=2, default=str))


if __name__ == "__main__":
    main()
