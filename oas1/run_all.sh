#!/usr/bin/env bash
set -euo pipefail

STUDY="${STUDY:-PRJNA785113}"
INDEX="${INDEX:-human.idx}"
T2G="${T2G:-human_t2g.txt}"

SCRATCH_DIR="$(mktemp -d -p . scratch.XXXX)"
trap 'rm -rf "${SCRATCH_DIR}"' EXIT

./00_fetch_manifest.py --study "${STUDY}" --out-runs manifest_runs.tsv --out-samples manifest_samples.tsv

SCRATCH_DIR="${SCRATCH_DIR}" INDEX="${INDEX}" T2G="${T2G}" ./01_build_kb_ref.sh

SCRATCH_DIR="${SCRATCH_DIR}" INDEX="${INDEX}" ./02_quant_stream_kallisto.sh manifest_samples.tsv

./03_build_h5ad.py --manifest manifest_samples.tsv --kallisto-dir kallisto --out gse190001_isoform.h5ad --layer-tpm
