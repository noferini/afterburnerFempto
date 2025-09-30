#!/bin/bash

# Usage:
#   ./scan_RD_R0.sh <input_folder> <n_jobs> <macro_name> <base_output_folder> [job_timeout_minutes]
#
# Example:
#   ./scan_RD_R0.sh input_dir 4 runAll.C test 45
#
# This will run runjobs.sh for:
#   R_D = 1.5, 1.75, ..., 3.5
#   R_0 = 0.80, 0.85, ..., 1.20
# and store results under:
#   <base_output_folder>/Results_R_0_<R0>_R_D_<RD>
# with dots replaced by underscores, e.g.
#   test/Results_R_0_0_8_R_D_1_5
#   test/Results_R_0_1_2_R_D_3_5
#
# It also creates a log file:
#   <base_output_folder>/<basename(base_output_folder)>_log.txt
# containing time, parameters, n_jobs, duration and messages from runjobs.sh.

input_folder="$1"
n_jobs="$2"
macro_name="$3"
base_output="$4"
job_timeout_minutes="${5:-0}"   # optional, in minutes; 0 → no timeout passed

if [[ -z "$input_folder" || -z "$n_jobs" || -z "$macro_name" || -z "$base_output" ]]; then
  echo "Usage: $0 <input_folder> <number_of_jobs> <macro_name> <base_output_folder> [job_timeout_minutes]"
  exit 1
fi

# Simple sanity check on n_jobs
if ! [[ "$n_jobs" =~ ^[0-9]+$ ]]; then
  echo "Error: <number_of_jobs> must be a positive integer."
  exit 1
fi

# Validate timeout minutes
if ! [[ "$job_timeout_minutes" =~ ^[0-9]+$ ]]; then
  echo "Error: [job_timeout_minutes] must be a non-negative integer."
  exit 1
fi

# Normalize base_output (remove trailing slash) and ensure it exists
base_output="${base_output%/}"
mkdir -p "$base_output"

# Prepare log file: <base_output>/<basename(base_output)>_log.txt
base_name="$(basename "$base_output")"
log_file="${base_output}/${base_name}_log.txt"

# Initial log header
{
  echo "======================================================="
  echo "scan_RD_R0.sh started at $(date '+%Y-%m-%d %H:%M:%S')"
  echo "Input folder : $input_folder"
  echo "Base output  : $base_output"
  echo "Macro        : $macro_name"
  echo "Number jobs  : $n_jobs"
  echo "Job timeout  : ${job_timeout_minutes} minutes (0 means no timeout passed to runjobs.sh)"
  echo "R_D values   : 1.625 1.875 2.125 2.375 2.625 2.875 3.125 3.375 3.625"
  echo "R_0 values   : 0.8 0.85 0.9 0.95 1.0 1.05 1.1 1.15 1.2"
  echo "======================================================="
  echo
} >> "$log_file"

# Main scan loops
for R_D in 1.625 1.875 2.125 2.375 2.625 2.875 3.125 3.375 3.625; do
  for R_0 in 0.8 0.85 0.9 0.95 1.0 1.05 1.1 1.15 1.2; do

    # Replace '.' with '_' for folder naming
    R0_tag="${R_0/./_}"
    RD_tag="${R_D/./_}"

    # Folder name: Results_R_0_<R0_tag>_R_D_<RD_tag>
    folder_name="Results_R_0_${R0_tag}_R_D_${RD_tag}"
    outdir="${base_output}/${folder_name}"

    mkdir -p "$outdir"

    echo
    echo "======================================================="
    echo "Running R_D = ${R_D}, R_0 = ${R_0}  →  ${outdir}"
    echo "======================================================="

    start_epoch=$(date +%s)
    start_human=$(date '+%Y-%m-%d %H:%M:%S')

    {
      echo "-------------------------------------------------------"
      echo "Start:      ${start_human}"
      echo "R_D:        ${R_D}"
      echo "R_0:        ${R_0}"
      echo "Output dir: ${outdir}"
      echo "n_jobs:     ${n_jobs}"
      echo "timeout:    ${job_timeout_minutes} minutes"
      echo "-------------------------------------------------------"
    } >> "$log_file"

    # Run runjobs.sh, appending all stdout/stderr into the log file
    if (( job_timeout_minutes > 0 )); then
      ./runjobs.sh "$input_folder" "$n_jobs" "$macro_name" "$outdir" "$R_D" "$R_0" "$job_timeout_minutes" >> "$log_file" 2>&1
    else
      ./runjobs.sh "$input_folder" "$n_jobs" "$macro_name" "$outdir" "$R_D" "$R_0" >> "$log_file" 2>&1
    fi

    status=$?

    end_epoch=$(date +%s)
    end_human=$(date '+%Y-%m-%d %H:%M:%S')
    duration=$(( end_epoch - start_epoch ))

    if (( status == 0 )); then
      result_txt="OK"
    else
      result_txt="FAILED (see messages above; some sub-jobs may have timed out or crashed)"
    fi

    {
      echo "Result:   ${result_txt}"
      echo "End:      ${end_human}"
      echo "Duration: ${duration}s"
      echo
    } >> "$log_file"

    echo "Finished R_D = ${R_D}, R_0 = ${R_0}  (status=${status}, duration=${duration}s)"

  done
done

echo
echo "All (R_D, R_0) points finished."
echo "Full log written to: ${log_file}"
