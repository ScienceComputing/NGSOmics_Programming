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

echo -e "$name, what species are you interested in your study? Please enter a digit to select from the following, e.g., 1"
select species in mouse human
do
    case $species in 
    mouse|human)   
            break
            ;;
    *)
            echo "Invalid choice" 
            ;;
    esac
done


file_1=../data_b/index.idx
file_2=../data_b/t2g.txt
if [ -f $file_1 ] && [ -f $file_2 ]
then
    echo -e "$name, we notice that $file_1 and $file_2 for $species exist.\n"
else 
    echo -e "$name, we are helping you download $species transcriptome index ...\n"
    kb ref -d $species -i "../data_b/index.idx" -g "../data_b/t2g.txt"
fi

echo -e "\n$name, what type of the count matrix file format would you like us to generate using kb python? Please enter a digit to select from the following, e.g., 1 "
select format in "h5ad (recommended)" "raw" "loom"
do
    case $format in 
    "h5ad (recommended)"|"raw"|"loom")   
            break
            ;;
    *)
            echo "Invalid format" 
            ;;
    esac
done
if [ "$format" = "h5ad (recommended)" ]
then
    f_choice="--h5ad"
elif [ "$format" = "raw" ]
then
    f_choice=""
elif [ "$format" = "loom" ]
then
    f_choice="--loom"
fi

echo -e "\nHere we list all the available single cell technologies used in kb python."
output=$(kb --list)
start_line=$(echo "$output" | grep -m 1 "^Custom")
echo "$output" | awk -v start_line="$start_line" 'BEGIN { print_flag = 0 } { if ($0 == start_line) print_flag = 1; if (print_flag) print }'
printf "\n"
read -p "$name, please select one of the above shorthand names of the single cell technology: e.g., 10XV3 " sc_tech
printf "\n"
read -p "$name, what maximum size of memory would you like to use in kb python: e.g., 32 " kb_memory
printf "\n"
read -p "$name, how many threads would you like to use in kb python: e.g., 20 " kb_thread

id=$(echo $sra_id | tr ',' '\n')
id_arr=($id)

for id_i in "${id_arr[@]}";
do  
    file_1=../result_b/"$id_i"/counts_unfiltered/cells_x_genes.mtx
    if [ -f $file_1 ]
    then
        echo -e "$name, we notice that $file_1 exist."
    else 
        echo -e "$name, we are helping you generate the cells x genes count matrix for $id_i ...~~~\n"
        kb count -i ../data_b/index.idx -g ../data_b/t2g.txt -x $sc_tech -m $kb_memory -t $kb_thread \
        ../data_b/"$id_i"_1.fastq.gz ../data_b/"$id_i"_2.fastq.gz $f_choice -o ../result_b/"$id_i" --overwrite
    fi
done