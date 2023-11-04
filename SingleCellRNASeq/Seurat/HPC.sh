# Log into HPC
ssh aphrodite
ssh curie.pbtech

# Create a job
touch scRNAseq_job.sh
vim scRNAseq_job.sh

# Show the description of all the nodes' CPU cores, memory (in Mb), runtime limits, and partition
sinfo -N -o "%25N %5c %10m %15l %25R"

# Show a description of the nodes in a given partition panda
sinfo -N -o "%25N %5c %10m %15l %25R" -p panda 

# Show how much memory is installed in each compute node affiliated to each partition
sinfo -N -l

# Show how much memory is already in use on each node. Look for the AllocMem entry
scontrol -o show nodes

# Show a simplified output of the nodes with only the first, thirteenth, and fourteenth fields displayed
scontrol -o show nodes | awk '{ print $1, $13, $14}'

# Submit the job
sbatch scRNAseq_job.sh

# Monitor the job status
squeue -u user_ID
squeue -u user_ID -l # Give more advanced output
