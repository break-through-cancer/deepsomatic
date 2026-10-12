"""
preprocess.py -- Cirro dataset preprocessing hook for DeepSomatic pipeline.
Runs on the Cirro platform BEFORE Nextflow starts, inspects the dataset's
file manifest, and injects `deepsomatic_runs` into params.

Shared-normal model: every tumor BAM found in the dataset is paired with
the SAME single normal BAM -- exactly one germline/normal sample per dataset.

Classification heuristic: any sample whose name contains "PBMC" is treated
as the normal; everything else is treated as a tumor. Adjust NORMAL_NAME_MARKERS
if normals use different naming.
"""
import json
from cirro.helpers.preprocess_dataset import PreprocessDataset

NORMAL_NAME_MARKERS = ["PBMC"]  # adjust if this cohort's normals aren't PBMC-labeled


# (data extension, index extension), in order of preference: BAM first (the
# original behaviour), then CRAM for datasets that only kept CRAMs.
ALIGNMENT_FORMATS = [(".bam", ".bai"), (".cram", ".crai")]

# Within one format, prefer recalibrated output, then the sorted/mapped
# alignment, then markduplicates output, then anything else.
NAME_PREFERENCE = [".recal.", ".sorted.", ".md."]


def _rank(path):
    name = path.rsplit("/", 1)[-1]
    for i, marker in enumerate(NAME_PREFERENCE):
        if marker in name:
            return i
    return len(NAME_PREFERENCE)


def extract_bams(ds):
    """Pick one alignment + matching index per sample (BAM, else CRAM).

    The keys in the returned dict are still called "bam"/"bai" because the
    workflow params (tumor_bam, tumor_bam_index, ...) use those names; the
    values may be CRAM/CRAI paths, which Parabricks DeepSomatic also accepts.
    """
    df = ds.files.copy()
    df["file"] = df["file"].astype(str)

    bam_map = {}

    for sample, group in df.groupby("sample"):
        files = set(group["file"])

        for data_ext, index_ext in ALIGNMENT_FORMATS:
            # Only accept a data file whose own index sits next to it
            pairs = [
                f for f in files
                if f.endswith(data_ext) and f + index_ext in files
            ]
            if pairs:
                best = sorted(pairs, key=lambda f: (_rank(f), f))[0]
                bam_map[str(sample)] = {"bam": best, "bai": best + index_ext}
                break

    if not bam_map:
        raise ValueError("No BAM/CRAM files with a matching index found in dataset")

    return bam_map


def is_normal(sample_name):
    upper = sample_name.upper()
    return any(marker in upper for marker in NORMAL_NAME_MARKERS)


def main():
    ds = PreprocessDataset.from_running()

    print("=== ds.files preview ===")
    print(ds.files.head(20).to_string(index=False))

    bam_map = extract_bams(ds)

    normal_bam = None
    normal_bai = None
    normal_sample_name = None
    tumor_runs = []

    for sample, files in bam_map.items():
        sample_name = str(sample)

        if is_normal(sample_name):
            if normal_bam is not None:
                raise ValueError(
                    f"Multiple normal-labeled samples found ({normal_sample_name!r} and "
                    f"{sample_name!r}) -- this pipeline expects exactly one shared normal "
                    f"per dataset. Adjust NORMAL_NAME_MARKERS or the dataset contents."
                )
            normal_bam = files["bam"]
            normal_bai = files["bai"]
            normal_sample_name = sample_name
        else:
            tumor_runs.append({
                "sample_name": sample_name,
                "tumor_bam": files["bam"],
                "tumor_bam_index": files["bai"],
            })

    if not tumor_runs:
        raise ValueError("No tumor BAM/CRAM found")

    if normal_bam is None:
        raise ValueError(
            "No normal-labeled sample found -- this pipeline requires a matched normal "
            "shared across all tumor samples (see NORMAL_NAME_MARKERS)"
        )

    deepsomatic_runs = []
    for run in tumor_runs:
        deepsomatic_runs.append({
            "sample_name": run["sample_name"],
            "tumor_bam": run["tumor_bam"],
            "tumor_bam_index": run["tumor_bam_index"],
            "normal_bam": normal_bam,
            "normal_bam_index": normal_bai,
        })

    ds.add_param("deepsomatic_runs", deepsomatic_runs)

    print("\nFinal parameters:")
    print(json.dumps(ds.params, indent=2, default=str))


if __name__ == "__main__":
    main()
