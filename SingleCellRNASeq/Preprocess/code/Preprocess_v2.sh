function display_welcome() {
    clear
    echo
    echo "***************************************************************************"
    echo "      Welcome to ABC Lab's Single Cell RNA-seq Preprocessing Pipeline      "
    echo "***************************************************************************"
    echo
}

function prompt_user() {
    read -p "Please enter your name: " name
    export name=$name
    echo
    echo "Hello, $name! We are happy to provide Single Cell RNA-seq Preprocessing service for you today!"
    echo
}

function service_entry_point() {
    echo
    echo "$name, what would you like us to help you today? Please enter a digit to select from the following, e.g., 1"
    select service in "Check your environment (recommended if you run this pipeline first time)" "Download FASTQ files from the NCBI Sequence Read Archive" "Psedoalign and count" "Perform basic QC" "Not interested in any service"
    do
        case $service in 
        "Check your environment (recommended if you run this pipeline first time)")   
                bash check_environment.sh
                break
                ;;
        "Download FASTQ files from the NCBI Sequence Read Archive")
                bash download_fastq.sh
                break
                ;;
        "Psedoalign and count")
                bash align_count.sh
                break
                ;;
        "Perform basic QC")
                bash perform_qc.sh
                break
                ;;
        "Not interested in any service")
                echo "Bye bye, $name. We wish you a nice day!"
                break
                ;; 
        *)
                echo "Invalid choice" 
                ;;
        esac
    done
    echo
}

function main() {
    display_welcome
    prompt_user
    service_entry_point
}

main


