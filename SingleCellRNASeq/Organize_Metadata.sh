#!/bin/bash

while getopts i:o:h:t:p:d:l:r:m:f: option
do
    case "${option}" in
        i) inputdir=${OPTARG};;
        o) outputdr=${OPTARG};;
        h) hostid=${OPTARG};;
        t) threads=${OPTARG};;
        p) platform=${OPTARG};;
        d) prime=${OPTARG};;
        f) read_length=${OPTARG};;
        c) covariate=${OPTARG};;
    esac
done
echo "input path: $inputdir"
echo "output path: $outputdr"
echo "host taxonomy id: $hostid"
echo "number of threads: $threads"
echo "single cell sequencing platform: $platform"
echo "barcode locations (on 5' or 3' end of reads): $prime"
echo "average read length: $read_length"
echo "possible confouders: $covariate"

Organize_Metadata.sh -o ~/scrnaseq/output
