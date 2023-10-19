#!/bin/bash

source ../check_b/.bashrc

if ! command -v prefetch &> /dev/null
then
    echo -e "$name, we notice that SRA toolkit is not found in your operating environment. You need to install this program or add the path of your SRA toolkit."
    read -p "$name, do you want to 1) install or 2) specify the path of your SRA toolkit? Please enter a digit, e.g., 1 " sra_choice
        if [ $sra_choice = 1 ]
        then 
            echo "You may visit https://github.com/ncbi/sra-tools/wiki/02.-Installing-SRA-Toolkit to complete the installation, and then rerun this bash script."
            exit
        elif [ $sra_choice = 2 ] 
        then
            read -p "$name, please enter the path of your SRA toolkit: e.g., /Users/ali4006/Documents/WCM_Project/202307_SingleCell/SRATool/sratoolkit.3.0.5-mac64/bin: " sra_path
            export PATH=$PATH:/$sra_path
            echo -e "$name, your SRA toolkit is executable.\n"
        fi
fi

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

echo
read -p "$name, how many threads would you like to use for downloading and compressing the FASTQ files: e.g., 20 " fast_thread
echo

# bash download_check.sh $sra_id $fast_thread

file="..data_b/*.fastq"
if [ -f "$file" ]
then
    echo -e "\n$name, we are helping you compress the FASTQ into FASTQ.gz.\n"
    pigz ../data_b/*.fastq -p $fast_thread
else 
    echo -e "\n$name, we notice that you have compressed FASTQ files at hand.\n"
fi