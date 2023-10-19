#!/bin/bash

echo -e "At the start, let us help you examine if SRA toolkit, python, and kb-python are installed in your computer."

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
            echo "export PATH=$PATH:/$sra_path" > ../check_b/.bashrc
            source ../check_b/.bashrc
            echo -e "$name, your SRA toolkit is executable.\n"
        fi
fi

echo
if [ ! command -v python &> /dev/null ] || [! command -v python3 &> /dev/null ]
then
    echo "$name, we notice that Python is not found in your operating environment. You need to install this program or add the path of your Python."
    read -p "$name, do you want to 1) install or 2) specify the path of your Python? Please enter a digit, e.g., 1 " py_choice
        if [ $py_choice = 1 ]
        then 
            echo "You may visit https://www.python.org/downloads/ to complete the installation, and then rerun this bash script."
            exit
        elif [ $py_choice = 2 ] 
        then
            read -p "$name, please enter the path of your python: e.g., /Library/Frameworks/Python.framework/Versions/3.11/bin/python3: " py_path
            export PATH=$PATH:/$py_path
        fi
fi
echo "$name, your Python is executable."
py3_version=$(python3 -V 2>&1)
py_version=$(python -V 2>&1)
err="command not found"
if [ $(echo $py_version | grep -c "$err") = 1 ]
then
    echo -e "Your python version is $py3_version\n"
elif [ $(echo $py3_version | grep -c "$err") = 1 ]
then
    echo -e "You are using $py_version\n"
else
    echo -e "You are using $py3_version\n"
fi

echo
if ! command -v kb &> /dev/null
then
    echo "$name, we notice that kb-python is not found in your operating environment. You need to install this program or add the path of your kb-python."
    read -p "$name, do you want to 1) install or 2) specify the path of your kb-python? Please enter a digit, e.g., 1 " kb_choice
        if [ $kb_choice = 1 ]
        then 
            echo "You may use https://www.kallistobus.tools/kb_usage/kb_usage/ to complete the installation, and then rerun this bash script."
            exit
        elif [ $kb_choice = 2 ] 
        then
            read -p "$name, please enter the path of your kb python: e.g., /Library/Frameworks/Python.framework/Versions/3.11/bin/kb: " kb_path
            export PATH=$PATH:/$kb_path
        fi
fi
echo "$name, your kb-python is executable."
kb_info=$(kb --info 2>&1)
echo -e "You are using $(echo "$kb_info" | grep -m 1 "^kb_python")"

echo
echo -e "In the next step, we are going to create folders for you to store the single cell analysis data and results."
dir_1=../data_b
dir_2=../result_b
dir_3=../result_b/figure
dir_4=../check_b
if [ -d "$dir_1" ] && [ -d "$dir_2" ] && [ -d "$dir_3" ] && [ -d "$dir_4" ]
then
    echo -e "$name, we notice that the directories $dir_1 and $dir_2 and $dir_3 and $dir_4 exist.\n"
else 
    echo -e "$name, we are helping you create $dir_1 and $dir_2 and $dir_3 and $dir_4 ...\n"
    mkdir -p ../data_b ../result_b/figure ../check_b
fi

echo
echo -e "$name, let us help you install python packages necessary for preprocessing single cell data in your environment."

pip install scanpy
pip install python-igraph
pip install louvain
pip install pybiomart