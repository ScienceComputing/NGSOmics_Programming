sra_list = ["SRR18743674", "SRR18743675", "SRR18743676", "SRR18743677", "SRR18743678"]
file_name = "embryo"
for sra_id in sra_list:
    prefetch = "prefetch {} -o {}".format(sra_id, work_dir + "/data/" + file_name + "/" + sra_id)
    print ("Currently downloading: " + sra_id)
    subprocess.call(prefetch, shell=True)

    fastq_dump = "fastq-dump {} --outdir {} --gzip --split-files".format(sra_id, file_name)
    print ("Generating fastq for: " + sra_id)
    subprocess.call(fastq_dump, shell=True)
