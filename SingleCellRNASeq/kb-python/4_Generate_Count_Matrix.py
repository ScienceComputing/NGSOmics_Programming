sra_list = ["SRR18743674", "SRR18743675", "SRR18743676", "SRR18743677", "SRR18743678"]
file_name = "embryo"
for sra_id in sra_list:
    kb_quant = "kb count -i mouse.idx -g mouse_t2g.txt -x 10xv3 -t 20 -m 32 \
                {0}/{1}_1.fastq.gz {0}/{1}_2.fastq.gz\
                -o ../result/{0}/raw_{1} --overwrite".format(file_name, sra_id)
    print ("Generating count quantification for: " + sra_id)
    subprocess.call(kb_quant, shell=True)
