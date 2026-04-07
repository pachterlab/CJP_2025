#!/usr/bin/env python3
"""Fetch ENA run metadata and build a streaming manifest for a study."""
from __future__ import annotations

import argparse
import io
import sys
import textwrap
import urllib.parse
import urllib.request

import pandas as pd

DEFAULT_FIELDS = [
    "run_accession",
    "fastq_ftp",
    "fastq_md5",
    "fastq_bytes",
    "library_layout",
    "library_strategy",
    "sample_accession",
    "sample_alias",
    "experiment_accession",
    "experiment_title",
    "study_accession",
]


def fetch_filereport(accession: str, fields: list[str]) -> pd.DataFrame:
    base = "https://www.ebi.ac.uk/ena/portal/api/filereport"
    params = {
        "accession": accession,
        "result": "read_run",
        "fields": ",".join(fields),
        "format": "tsv",
    }
    url = f"{base}?{urllib.parse.urlencode(params)}"
    try:
        with urllib.request.urlopen(url) as resp:
            data = resp.read().decode("utf-8")
    except Exception as e:
        raise RuntimeError(
            f"Failed to fetch ENA filereport from {url}. "
            "Check network access or try again later."
        ) from e

    if not data.strip():
        raise RuntimeError("ENA filereport response was empty.")

    df = pd.read_csv(io.StringIO(data), sep="\t")
    if "run_accession" not in df.columns:
        raise RuntimeError("ENA filereport response missing run_accession.")
    return df


def ftp_to_https(ftp_field: str) -> str:
    if pd.isna(ftp_field) or not str(ftp_field).strip():
        return ""
    urls = []
    for entry in str(ftp_field).split(";"):
        entry = entry.strip()
        if not entry:
            continue
        if entry.startswith("ftp://"):
            entry = entry.replace("ftp://", "https://", 1)
        elif entry.startswith("http://") or entry.startswith("https://"):
            pass
        else:
            entry = "https://" + entry
        urls.append(entry)
    return ";".join(urls)


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Fetch ENA filereport and build streaming manifest TSVs."
    )
    parser.add_argument(
        "--study",
        default="PRJNA785113",
        help="Study accession (PRJ/ERP/SRP/DRP). Default: PRJNA785113",
    )
    parser.add_argument(
        "--out-runs",
        default="manifest_runs.tsv",
        help="Output TSV with full ENA run metadata.",
    )
    parser.add_argument(
        "--out-samples",
        default="manifest_samples.tsv",
        help="Output TSV for quantification (sample_id, fastq_urls, layout, etc).",
    )
    args = parser.parse_args()

    df = fetch_filereport(args.study, DEFAULT_FIELDS)

    # Normalize fastq URLs for streaming over HTTPS.
    df["fastq_https"] = df.get("fastq_ftp", "").apply(ftp_to_https)

    # Write the full run report.
    df.to_csv(args.out_runs, sep="\t", index=False)

    # Build a sample-centric manifest using run_accession as sample_id.
    keep_cols = [
        "run_accession",
        "fastq_https",
        "library_layout",
        "library_strategy",
        "sample_accession",
        "sample_alias",
        "experiment_accession",
        "experiment_title",
        "study_accession",
        "fastq_md5",
        "fastq_bytes",
    ]
    for col in keep_cols:
        if col not in df.columns:
            df[col] = ""

    out = df[keep_cols].copy()
    out = out.rename(
        columns={
            "run_accession": "sample_id",
            "fastq_https": "fastq_urls",
        }
    )
    # Provide SRA run accession explicitly for SRA streaming mode.
    out.insert(3, "sra_run", out["sample_id"])
    out.to_csv(args.out_samples, sep="\t", index=False)

    print(
        textwrap.dedent(
            f"""
            Wrote: {args.out_runs}
            Wrote: {args.out_samples}
            Samples: {len(out)}
            """
        ).strip()
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
