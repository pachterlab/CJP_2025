#!/usr/bin/env python3
from __future__ import annotations

import argparse
import os
import sys
from typing import List

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp


def read_abundance(path: str) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t")
    if "target_id" not in df.columns:
        raise ValueError(f"Missing target_id in {path}")
    return df


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Build an h5ad from kallisto abundance.tsv files."
    )
    parser.add_argument(
        "--manifest",
        default="manifest_samples.tsv",
        help="Manifest TSV with sample_id and fastq_urls.",
    )
    parser.add_argument(
        "--kallisto-dir",
        default="kallisto",
        help="Directory containing per-sample kallisto outputs.",
    )
    parser.add_argument(
        "--out",
        default="gse190001_isoform.h5ad",
        help="Output h5ad filename.",
    )
    parser.add_argument(
        "--layer-tpm",
        action="store_true",
        help="Store TPM in adata.layers['tpm'] (default off).",
    )
    args = parser.parse_args()

    manifest = pd.read_csv(args.manifest, sep="\t")
    if "sample_id" not in manifest.columns:
        raise SystemExit("manifest missing sample_id column")

    sample_ids = manifest["sample_id"].tolist()

    first_path = os.path.join(args.kallisto_dir, sample_ids[0], "abundance.tsv")
    if not os.path.isfile(first_path):
        raise SystemExit(f"Missing abundance.tsv for first sample: {first_path}")

    first = read_abundance(first_path)
    transcripts = first["target_id"].astype(str).tolist()
    n_tx = len(transcripts)
    n_samples = len(sample_ids)

    counts = np.zeros((n_samples, n_tx), dtype=np.float32)
    tpm = np.zeros((n_samples, n_tx), dtype=np.float32) if args.layer_tpm else None

    for i, sid in enumerate(sample_ids):
        path = os.path.join(args.kallisto_dir, sid, "abundance.tsv")
        if not os.path.isfile(path):
            raise SystemExit(f"Missing abundance.tsv for sample {sid}: {path}")
        df = read_abundance(path)

        if df.shape[0] != n_tx or not (df["target_id"].astype(str).tolist() == transcripts):
            raise SystemExit(
                "Transcript lists differ across samples. "
                "Re-run with consistent index or sort/merge in code."
            )

        counts[i, :] = df["est_counts"].astype(np.float32).values
        if args.layer_tpm:
            tpm[i, :] = df["tpm"].astype(np.float32).values

    obs = manifest.set_index("sample_id")
    var = pd.DataFrame(index=pd.Index(transcripts, name="transcript_id"))

    adata = ad.AnnData(X=sp.csr_matrix(counts), obs=obs, var=var)
    if args.layer_tpm:
        adata.layers["tpm"] = sp.csr_matrix(tpm)

    adata.write_h5ad(args.out)
    print(f"Wrote {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
