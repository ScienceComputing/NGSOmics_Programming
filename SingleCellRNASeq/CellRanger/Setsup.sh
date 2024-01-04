cd opt/
sudo curl -o cellranger-7.2.0.tar.gz "https://cf.10xgenomics.com/releases/cell-exp/cellranger-7.2.0.tar.gz?Expires=1704439627&Key-Pair-Id=APKAI7S6A5RYOXBWRPDA&Signature=EfEgiGLoBDS-008YYSBIMroi4HTM5oloJWeHLnieZqqFSGrngHbP8GvnMNxkQ5wimO0j8to-RWG8lL2NKWujRcSwtb2V3zP8uwj753vZd~TCsreE0WduzFzpIEYwdiDlDotCll31DYCAlVCAdlJCq9tkBHkNJ1l0OSGmEe24H8xDKjVu6rvA5SqYIb-SsInDxK8XmbkBVZXZ5RT6DlZOPxQgr1PaEXv3L-EGd-Kuoox-GAh4B3I0XLwpZpzqRCEB8n70HAnk7snRg808gS7UjGbXvwUuP83MpnKebixxHusLx9c0NIJTRXqWyh-znC~w6HzyUsn266fqWhmBC35V-w__"
sudo curl -O "https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38-and-mm10-2020-A.tar.gz"
sudo tar -xzvf cellranger-7.2.0.tar.gz
export PATH=$PATH:/opt/cellranger-7.2.0
echo $PATH
echo "export PATH=$PATH:/opt/cellranger-7.2.0" > .bashrc
which cellranger
# /opt/cellranger-7.2.0/cellranger
