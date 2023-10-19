#!/bin/bash

sra_ids=()

folder_path="$1"

read_sra_ids_from_folder() {
    folder="$1"
    for file in "$folder"/SRR*.txt; do
        while IFS= read -r line || [[ -n "$line" ]]; do
            sra_ids+=("$line")
        done < "$file"
    done
}

read_sra_ids_from_folder "$folder_path"

sra_ids_csv=$(IFS=, ; echo "${sra_ids[*]}")
echo $sra_ids_csv
