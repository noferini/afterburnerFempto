#!/bin/bash

# Usage: ./runjobs.sh <input_folder> <number_of_jobs> <macro_name> <output_folder> <R_D> <R_0> [job_timeout_minutes]
#
# Example:
#   ./runjobs.sh input_dir 4 runAll.C output_dir 3.2 1.125 45
#
# If job_timeout_minutes is given and > 0, each coalroot job is wrapped in
#   timeout <job_timeout_minutes*60> ...
# so that individual jobs are killed if they exceed that time, without killing this script.

input_folder="$1"
n_jobs="$2"
macro_name="$3"
output_folder="$4"
R_D="$5"
R_0="$6"
job_timeout_minutes="$7"   # optional, in minutes

########################
#  Argument checking   #
########################

if [[ -z "$input_folder" || -z "$n_jobs" || -z "$macro_name" || -z "$output_folder" || -z "$R_D" || -z "$R_0" ]]; then
  echo "Usage: $0 <input_folder> <number_of_jobs> <macro_name> <output_folder> <R_D> <R_0> [job_timeout_minutes]"
  exit 1
fi

# Validate number of jobs
if ! [[ "$n_jobs" =~ ^[0-9]+$ ]]; then
  echo "Error: <number_of_jobs> must be a positive integer."
  exit 1
fi

# Validate timeout minutes if given
if [[ -n "$job_timeout_minutes" ]]; then
  if ! [[ "$job_timeout_minutes" =~ ^[0-9]+$ ]]; then
    echo "Error: [job_timeout_minutes] must be a non-negative integer."
    exit 1
  fi
fi

# Convert minutes → seconds
if [[ -n "$job_timeout_minutes" && "$job_timeout_minutes" -gt 0 ]]; then
  job_timeout=$(( job_timeout_minutes * 60 ))
else
  job_timeout=0
fi

# Normalize input and output folder paths (remove trailing slashes)
input_folder="${input_folder%/}"
output_folder="${output_folder%/}"

########################
#  CPU / timeout check #
########################

# Determine available CPU cores (macOS or Linux)
if command -v nproc &> /dev/null; then
  available_cores=$(nproc)
elif command -v sysctl &> /dev/null; then
  available_cores=$(sysctl -n hw.ncpu)
else
  echo "Warning: Cannot determine number of CPU cores. Defaulting to 1."
  available_cores=1
fi

if (( n_jobs > available_cores )); then
  echo "Error: Requested number of jobs ($n_jobs) exceeds available CPU cores ($available_cores)."
  echo "Please reduce the number of jobs to $available_cores or fewer."
  exit 1
fi

# Check if `timeout` is available (for per-job timeouts)
have_timeout=0
if command -v timeout &> /dev/null; then
  have_timeout=1
else
  if (( job_timeout > 0 )); then
    echo "Warning: 'timeout' command not found. Job timeout will be ignored."
  fi
fi

########################
#  Prepare directories #
########################

work_dir="$output_folder/job_lists"
mkdir -p "$work_dir"
mkdir -p "$output_folder"

# Save parameters once per output folder (not touched by hadd)
cat > "${output_folder}/parameters.txt" <<EOF
R_D=${R_D}
R_0=${R_0}
EOF

########################
#  Input file listing  #
########################

# Get all .root files
root_files=("$input_folder"/*.root)
total_files=${#root_files[@]}
if (( total_files == 0 )); then
  echo "No .root files found in input folder: ${input_folder}"
  exit 1
fi

files_per_job=$(( (total_files + n_jobs - 1) / n_jobs ))  # Ceiling division

echo "Total files:    $total_files"
echo "Jobs to create: $n_jobs"
echo "Files per job:  $files_per_job"
echo "Macro name:     $macro_name"
echo "Output folder:  $output_folder"
echo "R_D = ${R_D}, R_0 = ${R_0}"

if (( job_timeout > 0 && have_timeout == 1 )); then
  echo "Per-job timeout: ${job_timeout_minutes} minutes (${job_timeout}s)"
else
  echo "Per-job timeout: none"
fi

########################
#  Split into n lists  #
########################

for (( i=0; i<n_jobs; i++ )); do
  list_file="$work_dir/list_$i.txt"
  : > "$list_file"  # truncate if exists
  start=$((i * files_per_job))
  end=$((start + files_per_job))
  
  for (( j=start; j<end && j<total_files; j++ )); do
    echo "${root_files[j]}" >> "$list_file"
  done
done

########################
#  Launch parallel jobs#
########################

pids=()
declare -A job_index_for_pid
declare -A job_start_epoch

echo "Launching jobs..."

for (( i=0; i<n_jobs; i++ )); do
  list_file="$work_dir/list_$i.txt"
  out_folder="$output_folder/"
  output_file="output_$i.root"

  start_human=$(date '+%Y-%m-%d %H:%M:%S')
  start_epoch=$(date +%s)
  job_start_epoch["$i"]=$start_epoch

  echo "JOB $i START at ${start_human}"
  echo "  list    = ${list_file}"
  echo "  output  = ${out_folder}${output_file}"

  if (( job_timeout > 0 && have_timeout == 1 )); then
    echo "  timeout = ${job_timeout_minutes} minutes (${job_timeout}s)"
    timeout "$job_timeout" coalroot -l -b -q \
      "${macro_name}(\"${list_file}\", \"${out_folder}\", \"${output_file}\", ${R_D}, ${R_0})" &
  else
    echo "  timeout = none"
    coalroot -l -b -q \
      "${macro_name}(\"${list_file}\", \"${out_folder}\", \"${output_file}\", ${R_D}, ${R_0})" &
  fi

  pid=$!
  pids+=("$pid")
  job_index_for_pid["$pid"]=$i
done

########################
#  Wait & report       #
########################

overall_exit=0

for pid in "${pids[@]}"; do
  i=${job_index_for_pid["$pid"]}
  if [[ -z "$i" ]]; then
    i="unknown"
  fi

  if wait "$pid"; then
    status="OK"
    exit_code=0
  else
    exit_code=$?
    if (( exit_code == 124 && have_timeout == 1 && job_timeout > 0 )); then
      status="TIMEOUT"
    else
      status="FAILED"
    fi
    overall_exit=1
  fi

  end_human=$(date '+%Y-%m-%d %H:%M:%S')
  end_epoch=$(date +%s)

  start_epoch=${job_start_epoch["$i"]}
  if [[ -n "$start_epoch" ]]; then
    duration=$(( end_epoch - start_epoch ))
  else
    duration="unknown"
  fi

  echo "JOB $i END at ${end_human}"
  echo "  pid      = ${pid}"
  echo "  status   = ${status} (exit code = ${exit_code})"
  echo "  duration = ${duration}s"
done

if (( overall_exit != 0 )); then
  echo "Some jobs failed or timed out; continuing to merging step with whatever outputs were produced."
else
  echo "All jobs finished successfully."
fi

########################
#  Merging step        #
########################

echo "Merging files by stripped base name..."

cd "$output_folder" || exit 1

# Find all .root files and extract unique base names by removing _<number>.root
files=$(ls *.root 2>/dev/null)
if [[ -z "$files" ]]; then
  echo "No .root files produced in ${output_folder}"
  exit 1
fi

base_names=$(echo "$files" | sed -E 's/_[0-9]+\.root$//' | sort | uniq)

for base in $base_names; do
  pattern="${base}_[0-9]*.root"
  final="${base##*/}.root"
  
  echo "Merging files matching: $pattern → $final"
  hadd -f "$final" $pattern
done

echo "All merging complete. Final outputs are in: $output_folder"
