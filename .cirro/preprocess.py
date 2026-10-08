#!/usr/bin/env python3
"""
Cirro preprocess hook for DeepSomatic pipeline.

Converts Cirro input format (samples array) to Nextflow params format.
Validates that all BAM files exist and are properly indexed.
"""

import json
import sys
from pathlib import Path


def main():
    """Process Cirro inputs and generate Nextflow parameters."""

    # Read input from Cirro
    input_file = Path("/cirro/input/input.json")
    if not input_file.exists():
        print("ERROR: /cirro/input/input.json not found", file=sys.stderr)
        sys.exit(1)

    with open(input_file) as f:
        cirro_input = json.load(f)

    samples = cirro_input.get("samples", [])
    if not samples:
        print("ERROR: No samples provided", file=sys.stderr)
        sys.exit(1)

    # Validate and build deepsomatic_runs
    deepsomatic_runs = []
    for sample in samples:
        sample_name = sample.get("sample_name")
        tumor_bam = sample.get("tumor_bam")
        tumor_bam_index = sample.get("tumor_bam_index")
        normal_bam = sample.get("normal_bam")
        normal_bam_index = sample.get("normal_bam_index")

        # Validate required fields
        if not all([sample_name, tumor_bam, tumor_bam_index, normal_bam, normal_bam_index]):
            print(f"ERROR: Missing required fields for sample {sample_name}", file=sys.stderr)
            sys.exit(1)

        # Check if files exist (if local paths)
        for file_path in [tumor_bam, tumor_bam_index, normal_bam, normal_bam_index]:
            if not file_path.startswith("s3://"):
                path = Path(file_path)
                if not path.exists():
                    print(f"ERROR: File not found: {file_path}", file=sys.stderr)
                    sys.exit(1)

        deepsomatic_runs.append({
            "sample_name": sample_name,
            "tumor_bam": tumor_bam,
            "tumor_bam_index": tumor_bam_index,
            "normal_bam": normal_bam,
            "normal_bam_index": normal_bam_index,
        })

    # Build final params
    params = {
        "deepsomatic_runs": deepsomatic_runs,
    }

    # Write to output
    output_dir = Path("/cirro/output")
    output_dir.mkdir(parents=True, exist_ok=True)

    params_file = output_dir / "params.json"
    with open(params_file, "w") as f:
        json.dump(params, f, indent=2)

    print(f"✓ Preprocessed {len(deepsomatic_runs)} sample(s)")
    print(f"✓ Parameters written to {params_file}")


if __name__ == "__main__":
    main()
