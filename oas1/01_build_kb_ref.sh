#!/usr/bin/env bash
set -euo pipefail

INDEX="${INDEX:-human.idx}"
T2G="${T2G:-human_t2g.txt}"
REF_FASTA="${REF_FASTA:-}"
SCRATCH_DIR="${SCRATCH_DIR:-}"

if [[ -z "${SCRATCH_DIR}" ]]; then
  SCRATCH_DIR="$(mktemp -d -p . scratch.XXXX)"
  CLEAN_SCRATCH=1
else
  CLEAN_SCRATCH=0
fi

cleanup() {
  if [[ "${CLEAN_SCRATCH}" == "1" ]]; then
    rm -rf "${SCRATCH_DIR}"
  fi
}
trap cleanup EXIT

if [[ -f "${INDEX}" && -f "${T2G}" ]]; then
  echo "Index and t2g already exist: ${INDEX}, ${T2G}"
  exit 0
fi

if [[ -n "${REF_FASTA}" ]]; then
  if [[ ! -f "${REF_FASTA}" ]]; then
    echo "REF_FASTA not found: ${REF_FASTA}" >&2
    exit 1
  fi
  kb ref -f1 "${REF_FASTA}" -i "${INDEX}" -g "${T2G}"
else
  kb ref -d human -i "${INDEX}" -g "${T2G}"
fi

if [[ ! -s "${INDEX}" || ! -s "${T2G}" ]]; then
  echo "kb ref did not produce expected outputs: ${INDEX}, ${T2G}" >&2
  exit 1
fi

echo "Built: ${INDEX} and ${T2G}"
