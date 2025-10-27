#!/bin/bash

# Usage: ./runjobs.sh <input_folder> <number_of_jobs> <macro_name> <output_folder> <final_output_file.root>

input_folder="$1"
n_jobs="$2"
macro_name="$3"
output_folder="$4"

# Check arguments
if [[ -z "$input_folder" || -z "$n_jobs" || -z "$macro_name" || -z "$output_folder" ]]; then
  echo "Usage: $0 <input_folder> <number_of_jobs> <macro_name> <output_folder> "
  exit 1
fi

# Normalize input and output folder paths (remove trailing slashes)
input_folder="${input_folder%/}"
output_folder="${output_folder%/}"

# Determine available CPU cores (macOS or Linux)
if command -v nproc &> /dev/null; then
  available_cores=$(nproc)
elif command -v sysctl &> /dev/null; then
  available_cores=$(sysctl -n hw.ncpu)
else
  echo "Warning: Cannot determine number of CPU cores. Defaulting to 1."
  available_cores=1
fi

# Validate number of jobs
if ! [[ "$n_jobs" =~ ^[0-9]+$ ]]; then
  echo "Error: <number_of_jobs> must be a positive integer."
  exit 1
fi

if (( n_jobs > available_cores )); then
  echo "Error: Requested number of jobs ($n_jobs) exceeds available CPU cores ($available_cores)."
  echo "Please reduce the number of jobs to $available_cores or fewer."
  exit 1
fi

# Create working and output directories
work_dir="$output_folder/job_lists"
mkdir -p "$work_dir"
mkdir -p "$output_folder"

# Get all .root files
root_files=("$input_folder"/*.root)
total_files=${#root_files[@]}
files_per_job=$(( (total_files + n_jobs - 1) / n_jobs ))  # Ceiling division

echo "Total files: $total_files"
echo "Jobs to create: $n_jobs"
echo "Files per job: $files_per_job"
echo "Macro name: $macro_name"
echo "Output folder: $output_folder"

# Split files into n lists
for (( i=0; i<n_jobs; i++ )); do
  list_file="$work_dir/list_$i.txt"
  start=$((i * files_per_job))
  end=$((start + files_per_job))
  
  for (( j=start; j<end && j<total_files; j++ )); do
    echo "${root_files[j]}" >> "$list_file"
  done
done

# Run n ROOT jobs in parallel
for (( i=0; i<n_jobs; i++ )); do
  list_file="$work_dir/list_$i.txt"
  out_folder="$output_folder/"
  output_file="output_$i.root"
  echo "Launching job $i with list $list_file → $output_file"
  coalroot -l -b -q "${macro_name}(\"${list_file}\", \"${out_folder}\", \"${output_file}\")" &
done

wait
echo "All jobs finished."
echo "Merging files by stripped base name..."

cd "$output_folder" || exit 1

# Find all .root files and extract unique base names by removing _<number>.root
files=$(ls *.root 2>/dev/null)
base_names=$(echo "$files" | sed -E 's/_[0-9]+\.root$//' | sort | uniq)

for base in $base_names; do
  pattern="${base}_[0-9]*.root"
  final="${base##*/}.root"
  
  echo "Merging files matching: $pattern → $final"
  hadd -f "$final" $pattern
done

echo "All merging complete. Final outputs are in: $output_folder"
