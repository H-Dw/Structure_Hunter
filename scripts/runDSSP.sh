#!/bin/bash

if [ $# -ne 2 ]; then
    echo "Usage: $0 input_directory output_directory"
    exit 1
fi

input_directory="$1"
output_directory="$2"

mkdir -p "$output_directory/dssp"
mkdir -p "$output_directory/dssp_e"
mkdir -p "$output_directory/beta_strain"

script_dir="$(cd "$(dirname "$0")" && pwd)"

pdb_files=$(find "$input_directory" -type f -name "*.pdb")

for pdb_file in $pdb_files; do
    pdb_name="$(basename -- "$pdb_file")"
    # echo "Processing $pdb_name"
    
    dssp -i "$pdb_file" -o "$output_directory/dssp/$pdb_name.dssp"

    # Extract beta-sheet regions

    awk '$5 == "E"' "$output_directory/dssp/$pdb_name.dssp" > "$output_directory/dssp_e/beta_sheet_$pdb_name.dssp"
    awk 'FNR==NR{a[$1]=$0;next} ($6 in a) {print $0}' "$output_directory/dssp_e/beta_sheet_$pdb_name.dssp" "$pdb_file" > "$output_directory/beta_strain/beta_sheet_$pdb_name"
done
