#!/bin/bash

read -p "$name, do you want to 1) manually specify or 2) tell me the path of your folder that contains accession list file(s) downloaded from the NCBI SRA run selector (e.g., https://www.ncbi.nlm.nih.gov/Traces/study/?acc=SRP370066&o=acc_s%3Aa)? Please enter a digit, e.g., 1 " sra_id_choice
    if [ $sra_id_choice = 1 ]
    then 
        read -p "$name, please enter the SRA ID(s) of your interest: e.g., SRR18743674,SRR18743675,SRR18743676,SRR18743677,SRR18743678: " sra_id
    elif [ $sra_id_choice = 2 ] 
    then
        read -p "$name, please enter the path of your folder that contains accession list file(s): e.g., /Users/ali4006/Documents/WCM_Project/202307_SingleCell/Kallisto_Bustools/sra_id: " sra_id_path
        sra_id=$(bash organize_sra_ids.sh $sra_id_path)
        echo -e "Here are all your SRA IDs. \n$sra_id"
    fi

py3_version=$(python3 -V 2>&1)
py_version=$(python -V 2>&1)
err="command not found"
if [ $(echo $py_version | grep -c "$err") = 1 ]
then
    python3 QC.py "$sra_id"
elif [ $(echo $py3_version | grep -c "$err") = 1 ]
then
    python QC.py "$sra_id"
else
    python3 QC.py "$sra_id"
fi