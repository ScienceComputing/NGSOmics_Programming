#!/bin/bash

id_input="$1"
id_input=$(echo "$id_input" | tr ',' '\n')
id_arr=()

while read -r id_i; do
  id_arr+=("$id_i")
done <<< "$id_input"


fast_thread=$2

download() {
    sra_id=$1

    prefetch $sra_id -o ../data_b/$sra_id --max-size u
    fasterq-dump $sra_id -O ../data_b -e $fast_thread

    if [ $? -eq 0 ]; then
        echo "Completed downloading: $sra_id" >> ../check_b/fastq_checkpoint.txt
    else
        echo "Error occurred while processing $sra_id. The Preprocess script will resume from the last successful checkpoint."
    fi
}

check_fastq="../check_b/fastq_checkpoint.txt"
if [ -f "$check_fastq" ]
then
    last_success_id=$(tail -n 1 "$check_fastq" | cut -d ' ' -f3)

    idx=0
    for ((i=0; i<${#id_arr[@]}; i++))
    do
        if [ "${id_arr[i]}" = "$last_success_id" ]
        then
            idx=$((i + 1))
            break
        fi
    done

    for ((i=$idx; i<${#id_arr[@]}; i++))
    do
        download "${id_arr[i]}"
    done
else

    for id_i in "${id_arr[@]}"
    do
        download "$id_i"
    done
fi
