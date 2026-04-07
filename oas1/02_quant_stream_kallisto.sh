#!/usr/bin/env bash
set -euo pipefail

MANIFEST="${1:-manifest_samples.tsv}"
INDEX="${INDEX:-gencode.v49.idx}"
OUT_DIR="${OUT_DIR:-kallisto}"
THREADS="${THREADS:-8}"
N_JOBS="${N_JOBS:-1}"
FRAG_LEN="${FRAG_LEN:-200}"
FRAG_SD="${FRAG_SD:-20}"
RETRIES="${RETRIES:-2}"
USE_BIAS="${USE_BIAS:-1}" # 1 to pass --bias, 0 to disable
SCRATCH_DIR="${SCRATCH_DIR:-}"
MODE="${MODE:-ena}" # ena or sra
SRA_TOOL="${SRA_TOOL:-fasterq-dump}"
SRA_TEMP_ROOT="${SRA_TEMP_ROOT:-${SCRATCH_DIR}}"

if [[ ! -f "${MANIFEST}" ]]; then
  echo "Missing manifest: ${MANIFEST}" >&2
  exit 1
fi
if [[ ! -f "${INDEX}" ]]; then
  echo "Missing kallisto index: ${INDEX}" >&2
  exit 1
fi

mkdir -p "${OUT_DIR}"

if [[ -z "${SCRATCH_DIR}" ]]; then
  SCRATCH_DIR="$(mktemp -d -p . scratch.XXXX)"
  CLEAN_SCRATCH=1
else
  CLEAN_SCRATCH=0
fi
if [[ -z "${SRA_TEMP_ROOT}" ]]; then
  SRA_TEMP_ROOT="${SCRATCH_DIR}"
fi

cleanup() {
  if [[ "${CLEAN_SCRATCH}" == "1" ]]; then
    rm -rf "${SCRATCH_DIR}"
  fi
}
trap cleanup EXIT

if command -v pigz >/dev/null 2>&1; then
  DECOMPRESS_CMD="pigz -dc"
else
  DECOMPRESS_CMD="gzip -dc"
fi

run_one() {
  local sample_id="$1"
  local fastq_urls="$2"
  local layout="$3"
  local sra_run="$4"
  local bias_flag=""
  local attempt=1

  if [[ "${MODE}" == "ena" ]]; then
    if [[ -z "${fastq_urls}" ]]; then
      echo "No fastq URLs for ${sample_id}, skipping" >&2
      return 0
    fi
  else
    if [[ -z "${sra_run}" ]]; then
      if [[ "${sample_id}" =~ ^[SED]RR[0-9]+$ ]]; then
        sra_run="${sample_id}"
      else
        echo "No SRA run for ${sample_id}, skipping" >&2
        return 0
      fi
    fi
  fi

  local out_sample="${OUT_DIR}/${sample_id}"
  local log_sample="${out_sample}.log"
  if [[ -f "${out_sample}/abundance.tsv" ]]; then
    echo "Already done: ${sample_id}"
    return 0
  fi

  if [[ "${USE_BIAS}" == "1" ]]; then
    bias_flag="--bias"
  fi

  mkdir -p "${out_sample}"

  while [[ "${attempt}" -le "${RETRIES}" ]]; do
    : > "${log_sample}"
    echo "Running ${sample_id} attempt ${attempt}/${RETRIES}" | tee -a "${log_sample}"
    rm -f "${out_sample}/abundance.tsv" "${out_sample}/abundance.h5" "${out_sample}/run_info.json"

    if [[ "${MODE}" == "ena" ]]; then
      IFS=';' read -r -a url_arr <<< "${fastq_urls}"
      if [[ ${#url_arr[@]} -eq 0 ]]; then
        echo "No fastq URLs parsed for ${sample_id}, skipping" >&2
        return 0
      fi

      local ps_list=()
      local url
      for url in "${url_arr[@]}"; do
        url="${url//[$'\r\n ']/}"
        if [[ -z "${url}" ]]; then
          continue
        fi
        ps_list+=("<(curl -L --fail \"${url}\" | ${DECOMPRESS_CMD})")
      done

      if [[ "${layout}" == "PAIRED" ]]; then
        if [[ ${#ps_list[@]} -ne 2 ]]; then
          echo "Expected 2 FASTQ URLs for paired-end ${sample_id}, got ${#ps_list[@]}" >&2
          return 1
        fi
        eval "kallisto quant -i \"${INDEX}\" -o \"${out_sample}\" -t ${THREADS} ${bias_flag} ${ps_list[0]} ${ps_list[1]}" >> "${log_sample}" 2>&1
      else
        eval "kallisto quant -i \"${INDEX}\" -o \"${out_sample}\" -t ${THREADS} ${bias_flag} --single -l ${FRAG_LEN} -s ${FRAG_SD} ${ps_list[*]}" >> "${log_sample}" 2>&1
      fi
    else
      if ! command -v "${SRA_TOOL}" >/dev/null 2>&1; then
        echo "Missing SRA tool: ${SRA_TOOL}. Install sra-tools or set SRA_TOOL." >&2
        return 1
      fi
      if [[ "${layout}" == "PAIRED" ]]; then
        echo "SRA streaming for paired-end not supported in this script." >&2
        return 1
      fi
      local sra_tmp="${SRA_TEMP_ROOT}/sra_${sample_id}"
      mkdir -p "${sra_tmp}"
      eval "kallisto quant -i \"${INDEX}\" -o \"${out_sample}\" -t ${THREADS} ${bias_flag} --single -l ${FRAG_LEN} -s ${FRAG_SD} <(${SRA_TOOL} --stdout --temp \"${sra_tmp}\" ${sra_run})" >> "${log_sample}" 2>&1
      rm -rf "${sra_tmp}"
    fi

    if [[ -s "${out_sample}/abundance.tsv" ]]; then
      echo "Completed ${sample_id}" | tee -a "${log_sample}"
      return 0
    fi

    if [[ "${attempt}" -eq 1 && "${USE_BIAS}" == "1" ]]; then
      echo "Retrying ${sample_id} with USE_BIAS=0 fallback after failed attempt" | tee -a "${log_sample}"
      bias_flag=""
    fi
    attempt=$((attempt + 1))
    sleep 2
  done

  echo "Failed ${sample_id} after ${RETRIES} attempts" | tee -a "${log_sample}" >&2
  return 1
}

# Export needed for xargs -P
export -f run_one
export INDEX OUT_DIR THREADS FRAG_LEN FRAG_SD RETRIES USE_BIAS DECOMPRESS_CMD MODE SRA_TOOL SRA_TEMP_ROOT

# Read manifest and run (skip header)
# Columns: sample_id fastq_urls library_layout ...

if [[ "${N_JOBS}" -le 1 ]]; then
  tail -n +2 "${MANIFEST}" | while IFS=$'\t' read -r sample_id fastq_urls library_layout sra_run _rest; do
    run_one "${sample_id}" "${fastq_urls}" "${library_layout}" "${sra_run}"
  done
else
  tail -n +2 "${MANIFEST}" | awk -F'\t' '{print $1"\t"$2"\t"$3"\t"$4}' | \
    xargs -P "${N_JOBS}" -I{} bash -c 'IFS=$'"'\t'"' read -r sid urls lay sra <<< "{}"; run_one "$sid" "$urls" "$lay" "$sra"'
fi
