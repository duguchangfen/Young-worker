# Young-worker
基于质谱数据预测可翻译circRNA的工具
# Usage description
usage: python3 MScirctl1.py [-h] [-ci CIRC_INFO] [-cs CIRC_SEQ_FILE] [-cio CIRC_INFO_OMIT_FILE] [-ocf OUT_CIRC_FILE] [-ds DNA_SEQ] [-dg DNA_GFF] [-fc FIND_CIRC] [-ce CIRCEXPLORER] [-cI CIRI] [-sp SPECIES]

usage: python3 MScirctl2.py [-h] [-fi CIRC_FILE_INPUT] [-pit PATH_INPUT_TSV] [-pip PATH_INPUT_PEPTIDE] [-fip FILE_INPUT_PROTEIN] [-fims FILE_INPUT_PRIDE] [-ec EXCELLENT_COUNT] [-fo FILE_OUT]
# for example:
MScirctl below the examples are shown how to set paramaters:
~$python3.9 MScirctl1.py -cio ./data/circRNA_gff.txt -ds ./data/ath.fa -dg ./data/ath.gff3
~$python3 MScirctl2.py -fi data/MScirctl1_24-03-18.12/MScirctl1_main_result -pip data/peptide.txt -fip data/ath_protein.fa -fims data/MS.txt

# install requirement:
Mscirctl is available at https://github.com, and it can work well without being installed after you set the parameters properly, saving the time for installing in the operating system.

Mscirctl can predicting circRNA full-length sequence, predicting circRNA open reading frame and predicting translatable circRNA.

please transform your important input file format as same as the test input file recognized well by mscirctl. And there are not specific requirements for general input file ，such as genome sequences file (.fa, .fna format) and genome annotation file (.gtf or .gff format), but we recommend to download them from NCBI website.

We introduce IRESfinder (https://github.com/xiaofengsong/IRESfinder) into mstocirc for function IRES_predict, and explain on purpose what we do with their team tool IRESfinder to respect for them.

Thank you for your downloading and running this tool. If you have any questions，contact with us by email (glli@snnu.edu.cn).

Edited by Kai Wang on March 10th, 2024. Thank you....
