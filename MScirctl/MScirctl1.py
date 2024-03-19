#!/usr/bin/python3
###circRNA_mixture.py
import argparse
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.Seq import reverse_complement, transcribe, back_transcribe, translate
import re
import os
import time
from collections import Counter

parser = argparse.ArgumentParser('circRNA removes introns')
parser.add_argument('-ci','--circ_info',help='plantcircRNA数据库注释文件',default='none')
parser.add_argument('-cs','--circ_seq_file',help='plantcircRNA数据库circRNA序列文件',default='none')
parser.add_argument('-cio','--circ_info_omit_file',help='plantcircRNA数据库没有注释文件时，输入的自己生成的注释文件',default='none')
parser.add_argument('-ocf','--out_circ_file',help='输出结果文件的目录',default='./')
parser.add_argument('-ds','--dna_seq',help='对应物种的基因组序列文件')
parser.add_argument('-dg','--dna_gff',help='对应物种的基因组注释文件，gff3格式')
parser.add_argument('-fc','--find_circ',help='find_circ.py预测的结果文件',default='none')
parser.add_argument('-ce','--circexplorer',help='circexplorer2预测的结果文件',default='none')
parser.add_argument('-cI','--CIRI',help='CIRI预测的结果文件',default='none')
parser.add_argument('-sp','--species',help='预测的circRNA物种名称,不能包含标点符号',default='wk')
args = parser.parse_args()

circ_info = args.circ_info
circ_seq_file = args.circ_seq_file
out_circ_file = args.out_circ_file
circ_info_omit_file = args.circ_info_omit_file
dna_seq = args.dna_seq
dna_gff = args.dna_gff
find_circ = args.find_circ
circexplorer = args.circexplorer
CIRI = args.CIRI
species = args.species

def Create_folder():
    now=time.strftime('%y-%m-%d %H:%M:%S',time.localtime())
    ii=1
    out_circ_file = 'MScirctl1_'+now.split()[0]
    while(os.path.exists(out_circ_file+'.%d'%ii)):
        ii=ii+1
    out_circ_file = out_circ_file+'.%d'%ii + '/'
    os.system('mkdir %s'%(out_circ_file))
    return out_circ_file

if out_circ_file == './':
    out_circ_file = Create_folder()
elif out_circ_file != './' and out_circ_file[-1] != '/':
    out_circ_file = out_circ_file + '/'

file_dna_exon_s_gff = open(out_circ_file+'dna_exon_s.gff','w+')
print('结果文件输入',out_circ_file,'文件夹')

dic_gff_info_gene = {}
dic_gff_info_s = {}

def dic_gff_info_gene_s(circ_gff,circ_gene,circ_s):
    if circ_gff not in dic_gff_info_gene.keys():
        dic_gff_info_gene.update({circ_gff:circ_gene})
    elif circ_gff in dic_gff_info_gene.keys():
        if dic_gff_info_gene[circ_gff] == 'NA':
            dic_gff_info_gene.update({circ_gff:circ_gene})
    if circ_gff not in dic_gff_info_s.keys():
        dic_gff_info_s.update({circ_gff:circ_s})
    elif circ_gff in dic_gff_info_s.keys():
        if dic_gff_info_s[circ_gff] == 'NA':
            dic_gff_info_s.update({circ_gff:circ_s})
        elif dic_gff_info_s[circ_gff] != 'NA' and circ_s != 'NA' and dic_gff_info_s[circ_gff] != circ_s:
            dic_gff_info_s.update({circ_gff:'.'})

def find_circ_info_omit():
    file_find = open(find_circ,'r')
    for row_dind in file_find:
        circ_chr = row_dind.strip().split('\t')[0]
        circ_start = str(int(row_dind.strip().split('\t')[1]) + 1)
        circ_end = row_dind.strip().split('\t')[2]
        circ_gene = 'NA'
        circ_s = row_dind.strip().split('\t')[5]
        if circ_s != '-' and circ_s != '+':
            circ_s = 'NA'
        circ_gff = circ_chr + '\t' + circ_start + '\t' + circ_end + '\t' 
        dic_gff_info_gene_s(circ_gff,circ_gene,circ_s)

def circular_circ_info_omit():
    file_circular = open(circexplorer,'r')
    for row_circular in file_circular:
        circ_chr = row_circular.strip().split('\t')[0] 
        circ_start = str(int(row_circular.strip().split('\t')[1]) + 1)
        circ_end = row_circular.strip().split('\t')[2]
        circ_gene = row_circular.strip().split('\t')[14].split(':')[1]
        circ_s = row_circular.strip().split('\t')[5]
        if circ_s != '-' and circ_s != '+':
            circ_s = 'NA'
        circ_gff = circ_chr + '\t' + circ_start + '\t' + circ_end + '\t'
        dic_gff_info_gene_s(circ_gff,circ_gene,circ_s)

def ciri_circ_info_omit():
    file_ciri = open(CIRI,'r')
    for row_ciri in file_ciri :
        circ_chr = row_ciri.strip().split('\t')[1]
        if 'chr' == circ_chr:
            continue
        circ_start = row_ciri.strip().split('\t')[2]
        circ_end = row_ciri.strip().split('\t')[3]
        circ_gene = row_ciri.strip().split('\t')[9].strip(',')
        if 'n' in circ_gene and 'a' in circ_gene:
            circ_gene = 'NA'
        elif ',' in circ_gene:
            circ_gene = row_ciri.strip().split('\t')[9].strip(',').split(',')[0]
        circ_s = row_ciri.strip().split('\t')[10]
        if circ_s != '-' and circ_s != '+':
            circ_s = 'NA'
        circ_gff = circ_chr + '\t' + circ_start + '\t' + circ_end + '\t'
        dic_gff_info_gene_s(circ_gff,circ_gene,circ_s)

def circ_info_omit():
    file_info_omit = open(out_circ_file+'/circ_info_omit.gff','w+')
    num = '1'
    for dic_gff_info_s_key in dic_gff_info_s.keys():
        circ_num = num.rjust(6,'0')
        if dic_gff_info_s[dic_gff_info_s_key] == '.':
            row_info_omit = species + '_circ_' + circ_num + '\t' + dic_gff_info_s_key + dic_gff_info_gene[dic_gff_info_s_key] + '\t' + '+' + '\n' 
            file_info_omit.write(row_info_omit)
            num = str(int(num) + 1)
            circ_num = num.rjust(6,'0')
            row_info_omit = species + '_circ_' + circ_num + '\t' + dic_gff_info_s_key + dic_gff_info_gene[dic_gff_info_s_key] + '\t' + '-' + '\n'
            file_info_omit.write(row_info_omit)
            num = str(int(num) + 1)
        elif dic_gff_info_s[dic_gff_info_s_key] != '.':
            row_info_omit = species + '_circ_' + circ_num + '\t' + dic_gff_info_s_key + dic_gff_info_gene[dic_gff_info_s_key] + '\t' + dic_gff_info_s[dic_gff_info_s_key] + '\n' 
            file_info_omit.write(row_info_omit)
            num = str(int(num) + 1)
    file_info_omit.close()

if find_circ != 'none':
    find_circ_info_omit()
if circexplorer != 'none':
    circular_circ_info_omit()
if CIRI != 'none':
    ciri_circ_info_omit()
if dic_gff_info_gene != {} and dic_gff_info_s != {}:
    circ_info_omit()
    circ_info_omit_file = out_circ_file + '/circ_info_omit.gff'
#软件预测circRNA的结果文件circ_info_omit.gff文件生成完毕

file_dna_chr = os.popen("cat %s |grep '>'|cut -f 1 -d ' '|head -n 20"%(dna_seq)).read()
print(file_dna_chr)

file_dna_gff = open(dna_gff,'r')
file_dna_exon_gff = open(out_circ_file+'dna_exon.gff','w+')
file_dna_exon_s_gff = open(out_circ_file+'dna_exon_s.gff','w+')
file_dna_omit_gene_gff = open(out_circ_file+'dna_gene_id.gff','w+')
list_dna_gff_omit,line_dna_gff_omit_all,line_dna_gff_omit_s_all,line_dna_bed_omit_all = [],'','',''
line_dna_gff_omit_gene_all = ''
#print(file_dna_chr)
for line_dna_gff in file_dna_gff:
    if line_dna_gff[0] != '#':
        seq_type = line_dna_gff.strip().split('\t')[2]
        dna_chr = line_dna_gff.strip().split('\t')[0]
        if 'scaf' not in dna_chr and 'Mt' not in dna_chr and 'Pt' not in dna_chr and 'ChrC' not in dna_chr and 'ChrM' not in dna_chr and 'B73' not in dna_chr and 'iben' not in dna_chr and 'Sy' not in dna_chr and 'Un' not in dna_chr and 'KI' not in dna_chr and 'GL' not in dna_chr and 'MT' != dna_chr and 'X' != dna_chr and 'Y' != dna_chr:
            dna_chr = re.findall(r"\d+\.?\d*",dna_chr)[0]
            if dna_chr[0] == '0' and dna_chr[-1] != '0':
                dna_chr = dna_chr.strip('0')
            if dna_chr == '00':
                dna_chr = '0'
            if 'Chr' not in dna_chr and 'Chr' in file_dna_chr:
                dna_chr = 'Chr' + str(dna_chr)
        elif 'Mt' in dna_chr or 'Pt' in dna_chr or 'Sy' in dna_chr or 'Un' in dna_chr:
            dna_chr = dna_chr.replace('Chr','')
            if 'Chr' not in dna_chr and 'Chr' in file_dna_chr:
                dna_chr = 'Chr' + str(dna_chr)
        #print(dna_chr,'ddddddddddddddddddddddd')
        if seq_type == 'CDS' or seq_type == 'cds' or seq_type == 'exon':
            line_dna_gff_omit = dna_chr + '\t' + line_dna_gff.strip().split('\t')[3] + '\t' + line_dna_gff.strip().split('\t')[4]
            #line_dna_bed_omit = dna_chr + '\t' + str(int(line_dna_gff.strip().split('\t')[3])-1) + '\t' + line_dna_gff.strip().split('\t')[4]
            exon_direction = line_dna_gff.strip().split('\t')[6]
            if line_dna_gff_omit not in list_dna_gff_omit :
                list_dna_gff_omit.append(line_dna_gff_omit)
                file_dna_exon_gff.write(line_dna_gff_omit + '\n')
                line_dna_gff_omit_all = line_dna_gff_omit_all + line_dna_gff_omit + '\n'
                #line_dna_bed_omit_all = line_dna_bed_omit_all + line_dna_bed_omit + '\n'
                file_dna_exon_s_gff.write(line_dna_gff_omit + '\t' + exon_direction + '\n')
                #line_dna_gff_omit_s_all = line_dna_gff_omit_s_all + line_dna_gff_omit + '\t' + exon_direction + '\n'
        elif seq_type == 'gene':
            line_dna_gff_omit_gene = dna_chr + '\t' + line_dna_gff.strip().split('\t')[3] + '\t' + line_dna_gff.strip().split('\t')[4] + '\t' + line_dna_gff.strip().split('\t')[8].split(';')[1].split('=')[1] + '\n'
            file_dna_omit_gene_gff.write(line_dna_gff_omit_gene)
            #line_dna_gff_omit_gene_all = line_dna_gff_omit_gene_all + line_dna_gff_omit_gene

if line_dna_gff_omit_all == '':
    print('dna_gff file is error')

file_dna_exon_gff.close()
file_dna_exon_s_gff.close()
file_dna_omit_gene_gff.close()

#分类gff格式和bed格式的函数
def chr_bed_dff(chr,start,end):
    circ_bed = chr + '\t' + str(int(start) - 1 ) + '\t' + end + '\n'
    circ_gff = chr + '\t' + start + '\t' + end + '\n'
    return [circ_gff,circ_bed]

if circ_info == 'none' and circ_info_omit_file == 'none' and circ_seq_file != 'none':
    file_dna_seq = open(dna_seq,'r')
    dic_id_seq = {}
    dna_id_old = ''
    dna_seq_single = ''
    for line_dna_seq in file_dna_seq:
        if line_dna_seq[0] == '>':
            dna_id_new = line_dna_seq.strip().strip('>').split()[0]
            if dna_id_new != dna_id_old and dna_id_old != '':
                dic_id_seq.update({dna_id_old:dna_seq_single})
                #print(dna_seq_single[0:200])
                dna_seq_single = ''
        elif line_dna_seq[0] != '>':
            dna_seq_single_short = line_dna_seq.strip()
            dna_seq_single = dna_seq_single + dna_seq_single_short
            dna_id_old = dna_id_new

    dic_id_seq.update({dna_id_old:dna_seq_single})
    
    ##################################################
    
    file_seq = open(circ_seq_file,'r')
    line_circ_info_gff_all = ''
    line_circ_info_uniq_gff_all = ''
    line_circ_info_uniq_bed_all = ''
    line_circ_info_id_all = ''
    dic_id_gene = {}
    dic_circ_id_gff = {}
    for line_seq in file_seq :
        if line_seq[0] == '>':
            circ_id = line_seq.strip().strip('>').split(':')[0]
        elif line_seq[0] != '>':
            circ_seq = line_seq.strip()
            if 'N' in circ_seq or 'n' in circ_seq:
                continue
            circ_seq_find = 'no'
            circ_seq_short = circ_seq[0:20]
            for dna_id in dic_id_seq.keys():
                if circ_seq_short in dic_id_seq[dna_id]:
                    circ_seq_location = dic_id_seq[dna_id].find(circ_seq_short)
                    circ_seq_long = len(circ_seq)
                    circ_start_info = str(circ_seq_location+1)
                    circ_end_info = str(circ_seq_location + circ_seq_long)
                    dna_seq_slice = dic_id_seq[dna_id][circ_seq_location:int(circ_end_info)]
                    if dna_seq_slice == circ_seq:
                        line_circ_info_uniq_gff = dna_id + '\t' + circ_start_info + '\t' + circ_end_info
                        line_circ_info_id = line_circ_info_uniq_gff + '\t' + circ_id + '\n'
                        line_circ_info_id_all = line_circ_info_id_all + line_circ_info_id
                        dic_circ_id_gff.update({circ_id:line_circ_info_uniq_gff})
                        line_circ_info_uniq_gff_all = line_circ_info_uniq_gff_all + line_circ_info_uniq_gff + '\n'
                        line_circ_info_uniq_bed = dna_id + '\t' + str(int(circ_start_info) - 1) + '\t' + circ_end_info + '\n'
                        line_circ_info_uniq_bed_all = line_circ_info_uniq_bed_all + line_circ_info_uniq_bed
                        dic_id_gene.update({circ_id:'NA'})
                        circ_seq_find = 'yes'
                    elif dna_seq_slice != circ_seq:
                        if circ_seq in dic_id_seq[dna_id][circ_seq_location:]:
                            line_circ_info_uniq_gff = dna_id + '\t' + circ_start_info + '\t' + circ_end_info 
                            line_circ_info_id = line_circ_info_uniq_gff + '\t' + circ_id + '\n'
                            line_circ_info_id_all = line_circ_info_id_all + line_circ_info_id
                            dic_circ_id_gff.update({circ_id:line_circ_info_uniq_gff})
                            line_circ_info_uniq_gff_all = line_circ_info_uniq_gff_all + line_circ_info_uniq_gff + '\n'
                            line_circ_info_uniq_bed = dna_id + '\t' + str(int(circ_start_info) - 1) + '\t' + circ_end_info + '\n'
                            line_circ_info_uniq_bed_all = line_circ_info_uniq_bed_all + line_circ_info_uniq_bed
                            dic_id_gene.update({circ_id:'NA'})
                            circ_seq_find = 'yes'
                if circ_seq_find == 'yes':
                    break
            if circ_seq_find == 'no':
                print('>',circ_id,'不在基因组文件中')
                
    file_circ_info_id = open(out_circ_file+'info_id.gff','w')
    file_circ_info_id.write(line_circ_info_id_all)
    file_circ_info_id.close()
    os.system('bedtools intersect -a %s%s -b %s%s -wo > %s%s'%(out_circ_file,'info_id.gff',out_circ_file,'dna_gene_id.gff',out_circ_file,'circ_id_gene.gff'))
    file_circ_id_gene = open(out_circ_file+'circ_id_gene.gff','r')
    dic_id_overlap = {}
    for line_circ_id_gene in file_circ_id_gene:
        overlap = line_circ_id_gene.strip().split('\t')[8]
        gene = line_circ_id_gene.strip().split('\t')[7]
        circ_id = line_circ_id_gene.strip().split('\t')[3]
        if circ_id not in dic_id_overlap.keys():
            dic_id_overlap.update({circ_id:overlap})
            dic_id_gene.update({circ_id:gene})
        elif circ_id in dic_id_overlap.keys():
            if int(overlap) > int(dic_id_overlap[circ_id]):
                dic_id_overlap.update({circ_id:overlap})
                dic_id_gene.update({circ_id:gene})
    
    file_134589 = open(out_circ_file+'circ_info_omit.gff','w+')
    line_circ_id_all = ''
    for circ_id in dic_id_gene.keys():
        line_circ_id = circ_id + '\t' + dic_circ_id_gff[circ_id] + '\t' + dic_id_gene[circ_id] + '\tNA\n'
        file_134589.write(line_circ_id)
        #line_circ_id_all = line_circ_id_all + line_circ_id
    
    #file_134589 = open(out_circ_file+'circ_info_omit.gff','w')
    #file_134589.write(line_circ_id_all)
    file_134589.close()
                    
elif circ_info != 'none' and circ_info_omit_file == 'none':
    file_circ_info = open(circ_info,'r')
    line_circ_info_gff_all = ''
    line_circ_info_uniq_gff_all = ''
    line_circ_info_uniq_bed_all = ''
    dic_circ_id_start_seq = {}
    dic_circ_id_end_seq = {}
    dic_id_gene = {}
    dic_circ_bed_circ_id = {}
    for line_circ_info in file_circ_info:
        circ_id = line_circ_info.strip().split('\t')[0]
        if 'circ' not in circ_id:
            continue
        circ_chr_info = line_circ_info.strip().split('\t')[2]
        circ_start_info = line_circ_info.strip().split('\t')[3]
        circ_end_info = line_circ_info.strip().split('\t')[4]
        circ_gene = line_circ_info.strip().split('\t')[7]
        circ_direction = line_circ_info.strip().split('\t')[8]
        circ_start_seq = line_circ_info.strip().split('\t')[11]
        circ_end_seq = line_circ_info.strip().split('\t')[10]
        dic_id_gene.update({circ_id:circ_gene})
        dic_circ_id_start_seq.update({circ_id:circ_start_seq})
        dic_circ_id_end_seq.update({circ_id:circ_end_seq})
        if 'scaf' not in circ_chr_info and 'A' not in circ_chr_info and 'D' not in circ_chr_info and 'iben' not in circ_chr_info and 'Sy' not in circ_chr_info and 'Un' not in circ_chr_info:
            #print(circ_chr_info)
            if 'Mt' in circ_chr_info or 'Pt' in circ_chr_info:
                circ_chr_info = circ_chr_info
            elif 'Mt' not in circ_chr_info and 'Pt' not in circ_chr_info :
                circ_chr_info = re.findall(r"\d+\.?\d*",circ_chr_info)[0]
        if circ_chr_info[0] == '0' and circ_chr_info[-1] != '0':
            circ_chr_info = circ_chr_info.strip('0')
        if circ_chr_info == '00':
            circ_chr_info = '0'
        if 'Chr' not in circ_chr_info and 'Chr' in file_dna_chr:
            circ_chr_info = 'Chr' + str(circ_chr_info)
        elif 'Chr' in circ_chr_info and 'Chr' not in file_dna_chr:
            circ_chr_info = circ_chr_info.replace('Chr','')

        line_circ_info_gff = circ_id + '\t' + circ_chr_info + '\t' + circ_start_info + '\t' + circ_end_info + '\t' + circ_gene + '\t' + circ_direction + '\n'
        line_circ_info_gff_all = line_circ_info_gff_all + line_circ_info_gff
        line_circ_info_uniq_gff = circ_chr_info + '\t' + circ_start_info + '\t' + circ_end_info + '\n'
        line_circ_info_uniq_bed = circ_chr_info + '\t' + str(int(circ_start_info)-1) + '\t' + circ_end_info + '\n'
        line_circ_info_uniq_gff_all = line_circ_info_uniq_gff_all + line_circ_info_uniq_gff
        line_circ_info_uniq_bed_all = line_circ_info_uniq_bed_all + line_circ_info_uniq_bed
        circ_bed_id = circ_chr_info + ':' + str(int(circ_start_info)-1) + '-' + circ_end_info 
        dic_circ_bed_circ_id.update({circ_bed_id:circ_id})
    file_circ_info.close()
    #print(line_circ_info_gff_all[0:200])
    file_134589 = open(out_circ_file+'circ_info_omit.gff','w')
    file_134589.write(line_circ_info_gff_all)
    file_134589.close()
    
elif circ_info == 'none' and circ_info_omit_file != 'none' :
    file_info_omit_self = open(circ_info_omit_file,'r')
    line_circ_info_uniq_gff_all,line_circ_info_uniq_bed_all = '',''
    dic_id_gene = {}
    dic_circ_bed_circ_id = {}
    dic_circ_id_end_seq = {}
    for line_info_omit_self in file_info_omit_self:
        circ_id = line_info_omit_self.strip().split('\t')[0]
        circ_chr = line_info_omit_self.strip().split('\t')[1]
        circ_start = line_info_omit_self.strip().split('\t')[2]
        circ_end = line_info_omit_self.strip().split('\t')[3]
        circ_gene = line_info_omit_self.strip().split('\t')[4]
        dic_id_gene.update({circ_id:circ_gene})
        circ_bed_id = circ_chr + ':' + str(int(circ_start)-1) + '-' + circ_end
        dic_circ_bed_circ_id.update({circ_bed_id:circ_id})
        line_circ_info_uniq_gff_all = line_circ_info_uniq_gff_all + chr_bed_dff(circ_chr,circ_start,circ_end)[0]
        line_circ_info_uniq_bed_all = line_circ_info_uniq_bed_all + chr_bed_dff(circ_chr,circ_start,circ_end)[1]
    file_info_omit_self.close()
        
file_circ_uniq_gff = open(out_circ_file+'circ_uniq.gff','w')
file_circ_uniq_gff.write(line_circ_info_uniq_gff_all)
file_circ_uniq_gff.close()

file_circ_uniq_bed = open(out_circ_file+'circ_uniq.bed','w')
file_circ_uniq_bed.write(line_circ_info_uniq_bed_all)
file_circ_uniq_bed.close()
print('info文件生成')


###########################################################################

if circ_seq_file == 'none':
    os.system('bedtools getfasta -fi %s -bed %s%s -fo %s%s '%(dna_seq,out_circ_file,'circ_uniq.bed',out_circ_file,'circ_seq.fasta'))
    line_id_seq_all = ''
    file_circ_info = open(out_circ_file+'circ_seq.fasta','r')
    for line_circ_info in file_circ_info:
        if line_circ_info[0] == '>':
            circ_bed_new_id = line_circ_info.strip().strip('>')
            circ_new_id = dic_circ_bed_circ_id[circ_bed_new_id]
        elif line_circ_info[0] != '>':
            circ_new_seq = line_circ_info.strip().upper()
            if 'N' in circ_new_seq:
                continue
            if dic_circ_id_end_seq == {} :
                line_id_seq = '>' + circ_new_id + '\n' + circ_new_seq + '\n'
                line_id_seq_all = line_id_seq_all + line_id_seq
                continue
            if circ_new_seq.find(dic_circ_id_start_seq[circ_new_id]) == 0 or circ_new_seq.find(dic_circ_id_end_seq[circ_new_id]) == 0 :
                line_id_seq = '>' + circ_new_id + '\n' + circ_new_seq + '\n'
                line_id_seq_all = line_id_seq_all + line_id_seq
            else :
                print('基因组文件可能与circ_info文件不匹配')
                print(circ_bed_new_id)
                print(circ_new_id)
                print(circ_new_seq)
                print(dic_circ_id_start_seq[circ_new_id])
                print(dic_circ_id_end_seq[circ_new_id])
                
    circ_seq_file = out_circ_file + 'circ_slef_seq.fasta'
    file_circ_seq_new = open(out_circ_file+'circ_slef_seq.fasta','w')
    file_circ_seq_new.write(line_id_seq_all)
    file_circ_seq_new.close()

#file_circ_id_seq = open('/home/li/omics/wk/procedure/try_exon/all_circMeta_DECs','r')
file_circ_id_seq = open(circ_seq_file,'r')
file_id_gene_seq_pr = open(out_circ_file+'plant_circ_pr.txt','w+')
#line_new = ''
list_circ_unknown = []
for line_circ_id_seq in file_circ_id_seq:
    if '>' in line_circ_id_seq:
        #circ_order = '>' + dic_id_order[line_circ_id_seq.strip().split('\t')[0]]
        circ_id = line_circ_id_seq.strip().strip('>').split(':')[0]
    if '>' not in line_circ_id_seq:
        circ_seq = line_circ_id_seq.strip()
        if 'N' not in circ_seq:
            if circ_id in dic_id_gene.keys():
                circ_order_p = '>' + circ_id + ':' + dic_id_gene[circ_id] + '_P'
                circ_order_r = '>' + circ_id + ':' + dic_id_gene[circ_id] + '_R'
                circ_seq_r = reverse_complement(circ_seq)
            elif circ_id not in circ_seq:
                circ_order_p = '>' + circ_id + ':' + 'none' + '_P'
                circ_order_r = '>' + circ_id + ':' + 'none' + '_R'
                circ_seq_r = reverse_complement(circ_seq)
                list_circ_unknown.append(circ_order_p.strip().strip('>').strip('_P'))
            file_id_gene_seq_pr.write(circ_order_p + '\n' + circ_seq + '\n' + circ_order_r + '\n' + circ_seq_r + '\n')
            #line_new = line_new + circ_order_p + '\n' + circ_seq + '\n' + circ_order_r + '\n' + circ_seq_r + '\n'

file_id_gene_seq_pr.close()

####################################################################################################

#需要安装bedtools
os.system('bedtools intersect -a %s%s -b %s%s -wa -wb > %s%s'%(out_circ_file,'circ_uniq.gff',out_circ_file,'dna_exon.gff',out_circ_file,'circ_exon.gff'))


###circ_include_exon_gff.py
#file1 = open('/home/li/public2/wk/work/Oryza_sativa_japonica_plantcircbace/circ_exon.gff','r')
file1 = open(out_circ_file+'circ_exon.gff','r')

x = ''
c = ''
w = ''
for line1 in file1:
    a = line1.strip().split('\t')
    b = a[0] + '_' + a[1] + '_' + a[2]        #circRNA的染色体数和起始终止位置
    if b != x:
        if c != '':
            l = c + '\n'
            w = w + l               #每个circRNA上对应的外显子序列的字符串
        c = b + '\t' + a[0] + '\t' + a[-2] + '-' + a[-1]
        x = b
    elif b == x:
        c = c + ':' + a[-2] + '-' + a[-1]
#file2 = open('/home/li/public2/wk/work/Oryza_sativa_japonica_plantcircbace/circ_exon_exon.gff','w')
file2 = open(out_circ_file+'circ_exon_exon.gff','w')
file2.write(w)
file2.close()

####################################################################################################

###circ_include_exon_gff_uniq.py
#from collections import Counter

def func(lst):
	new_list = []
	for i,j in zip(lst,lst[1:]):
		if j - i > 1:
			new_list.append(lst[:lst.index(j)])
			lst = lst[lst.index(j):]
	new_list.append(lst)
	return new_list
#circRNA匹配到的外显子的正负信息
#file_exon_s = open('/home/li/public2/wk/work/Oryza_sativa_japonica_phytozome/osa_exon_s_sort.gtf','r')
file_exon_s = open(out_circ_file+'dna_exon_s.gff','r')

new_dic = {}
for line_exon_s in file_exon_s:
    new_exon_s = line_exon_s.strip().split('\t')
    new_exon = new_exon_s[0] + '\t' + new_exon_s[1] + '\t' + new_exon_s[2]
    new_s = new_exon_s[3]
    new_dic.update({new_exon:new_s})


if circ_info_omit_file != 'none':
    file_circ_g = open(circ_info_omit_file,'r')
elif circ_info_omit_file == 'none':
    file_circ_g = open(out_circ_file+'circ_info_omit.gff','r')
new_circ_g_dic,new_circ_ID_dic, = {},{}
new_circ_dic = {}
for line_circ_g in file_circ_g:
    new_circ_g = line_circ_g.strip().split('\t')
    new_circ = new_circ_g[1] + '_' + new_circ_g[2] + '_' + new_circ_g[3]
    new_circ_gene = new_circ_g[4]
    new_circ_g_dic.update({new_circ:new_circ_gene})
    new_circ_ID = new_circ_g[0]
    new_circ_ID_dic.update({new_circ:new_circ_ID})
    new_circ_s = new_circ_g[5]
    new_circ_dic.update({new_circ:new_circ_s})

#circ_exon = open('/home/li/public2/wk/work/Oryza_sativa_japonica_plantcircbace/circ_exon_exon.gff','r')
circ_exon = open(out_circ_file+'circ_exon_exon.gff','r')
file_circ_exon_multiple = open(out_circ_file+'circ_exon_uniqer.gff','w+')
file_bed_exon_multiple = open(out_circ_file+'circ_exon.bed','w+')

circ_exon_multiple = ''
ch_exon_multiple,circ_bed_exon_multiple = '',''
for line1 in circ_exon:
    a = line1.strip().split('\t')
    circRNA_s = new_circ_dic[a[0]]
    #circRNA的正负
    ch = a[0].strip().rsplit('_',2)[-3]
    circ_start = a[0].strip().rsplit('_',2)[-2]
    circ_end = a[0].strip().rsplit('_',1)[-1]
    exon_ch = a[1]
    exon_start_end_multiple = a[2].split(':')
    exon_count = len(exon_start_end_multiple)
    set_exon_multiple = set()
    if exon_count == 1:
        ch_exon = exon_ch + '_' + circ_start + '_' +circ_end + '\t' + 'exon:1'
        exon_bed_strat_end = exon_ch + '\t' + a[2].split('-')[0] + '\t' + a[2].split('-')[1]
        ch_bed_ch_strat_end = exon_ch + '\t' + circ_start + '\t' + circ_end 
        ch_bed_exon = exon_ch + '\t' + str(int(circ_start)-1) + '\t' +circ_end + '\t.\t.\t' + circRNA_s + '\t' + new_dic[exon_bed_strat_end] + '\t.\t.'+ '\n'
        #外显子六列数据的bed文件
    elif exon_count != 1:
        ch_exon_multiple,exon_s = '',''
        for exon in exon_start_end_multiple:
            exon_start = exon.split('-')[0]
            exon_end = exon.split('-')[1]
            exon_ch_strat_end = exon_ch + '\t' + exon_start + '\t' + exon_end 
            exon_s = exon_s + new_dic[exon_ch_strat_end]
            #exon正负链的链表
            #发现exon的正负链有问题
            #circRNA的负链是反向的序列
            list_exon = range(int(exon_start),int(int(exon_end)+1))
            set_exon = set(list_exon)
            set_exon_multiple = set.union(set_exon_multiple,set_exon)
        exon_count = Counter(exon_s)
        exon_count_key = exon_count.keys()
        if len(exon_count_key) == 1:
            if '+' in exon_count_key:
                exon_column_4 = exon_count['+']
                exon_column_5 = '0'
                exon_column_6 = '+'
            elif '-' in exon_count_key:
                exon_column_4 = '0'
                exon_column_5 = exon_count['-']
                exon_column_6 = '-'
        elif len(exon_count_key) == 2:
            exon_column_4 = exon_count['+']
            exon_column_5 = exon_count['-']
            if exon_column_4 >= exon_column_5:
                exon_column_6 = '+'
            elif exon_column_4 < exon_column_5:
                exon_column_6 = '-'
        #得到该circRNA的正负链信息
        list_exon_multiple = sorted(list(set_exon_multiple))
        exon_count_uniq = func(list_exon_multiple)
        exon_count_uniq_count = len(exon_count_uniq)
        exon_first,exon_mid,exon_last = '','',''
        exon_bed_first,exon_bed_last,exon_bed_mid = '','',''
        if exon_count_uniq_count == 1:
            ch_exon = exon_ch + '_' + circ_start + '_' +circ_end + '\t' + 'exon:1'
            ch_bed_ch_strat_end = exon_ch + '\t' + str(min(exon_count_uniq[0])) + '\t' + str(max(exon_count_uniq[0]))
            ch_bed_exon = exon_ch + '\t' + str(int(circ_start)-1) + '\t' +circ_end + '\t.\t.\t%s\t%s\t%s\t%s\n'%(circRNA_s,exon_column_6,exon_column_4,exon_column_5)
            #外显子六列数据的bed文件
        elif exon_count_uniq_count != 1:
            for exon_count_2 in range(0,exon_count_uniq_count):
                exon_single = exon_count_uniq[exon_count_2]
                end = int(exon_count_uniq_count) - 1
                if int(exon_count_2) == 0:
                    min1 = circ_start
                    max1 = max(exon_single)
                    exon_first = exon_first  + exon_ch + '_' + str(min1) + '_' + str(max1) + ';'
                    exon_bed_first = exon_bed_first + exon_ch + '\t' + str(int(min1)-1) + '\t' + str(max1) + '\t.\t.\t%s\t%s\t%s\t%s\n'%(circRNA_s,exon_column_6,exon_column_4,exon_column_5)
                    if int(min1) == max1:
                        exon_first,exon_bed_first = '',''
                        print('1',a[0])
                elif int(exon_count_2) == end:
                    min3 = min(exon_single)
                    max3 = circ_end
                    exon_last = exon_last + exon_ch + '_' + str(min3) + '_' + str(max3)
                    exon_bed_last = exon_bed_last + exon_ch + '\t' + str(min3-1) + '\t' + str(max3) + '\t.\t.\t%s\t%s\t%s\t%s\n'%(circRNA_s,exon_column_6,exon_column_4,exon_column_5)
                    if min3 == int(max3):
                        print('3',a[0])
                else:
                    min2 = min(exon_single)
                    max2 = max(exon_single)
                    if int(min1) == max1 and int(exon_count_2)==1:
                        min2 = int(min1)
                        print('11',a[0])
                    exon_mid = exon_mid + exon_ch + '_' + str(min2) + '_' + str(max2) + ';'
                    exon_ch_start_end = exon_ch + '\t' + str(min2-1) + '\t' + str(max2) 
                    exon_bed_mid = exon_bed_mid + exon_ch_start_end + '\t.\t.\t%s\t%s\t%s\t%s\n'%(circRNA_s,exon_column_6,exon_column_4,exon_column_5)
                    if min2 == max2:
                        print('2',a[0])
            ch_exon = exon_first + exon_mid + exon_last + '\t' +'exon:%s'%(exon_count_uniq_count)
            ch_bed_exon = exon_bed_first + exon_bed_mid + exon_bed_last
    if circRNA_s == 'NA':
        circRNA_z_f = '+'
    elif circRNA_s != 'NA':
        circRNA_z_f = circRNA_s
    file_circ_exon_multiple.write(a[0] + '\t' + ch_exon + '\t' + circRNA_z_f + '\n')
    file_bed_exon_multiple.write(ch_bed_exon)
    #circ_exon_multiple = circ_exon_multiple + a[0] + '\t' + ch_exon + '\t' + circRNA_z_f + '\n'
    #circ_bed_exon_multiple = circ_bed_exon_multiple + ch_bed_exon

#file_circ_exon_multiple = open(out_circ_file+'circ_exon_uniqer.gff','w')
#file_circ_exon_multiple.write(circ_exon_multiple)
file_circ_exon_multiple.close()

#file_bed_exon_multiple = open(out_circ_file+'circ_exon.bed','w')
#进行了改变在这一步将包含的exon文件变成了按照circRNA起始位置和终止位置的序列,进行列起始位置减一的操作.
#file_bed_exon_multiple.write(circ_bed_exon_multiple)
file_bed_exon_multiple.close()



#需要安装bedtools
os.system('bedtools getfasta -fi %s -bed %s%s -fo %s%s '%(dna_seq,out_circ_file,'circ_exon.bed',out_circ_file,'dna_exon_bed.fasta'))

####################################################################################################
#Osa_CDS_1_bed.fasta文件需要用bedtools工具生成
###circRNA_exon_connect.py
e = {}
bsj_id_seq_m = ''
#bed文件得出的fasta文件所有需要起始位置再加一和gff文件统一
#for seq_record in SeqIO.parse("/home/li/public2/wk/work/Oryza_sativa_japonica_phytozome/Osa_CDS_1_bed.fasta", "fasta"):
for seq_record in SeqIO.parse(out_circ_file+'dna_exon_bed.fasta','fasta'):
    a = seq_record.id 
    x = re.split('[:-]', a)[0] + "_" + str(int(re.split('[:-]', a)[1]) + 1) + "_" + re.split('[:-]', a)[2].split('(')[0]
    c = seq_record.seq
    b = str(c.upper())               #正链
    d = {x:b}
    e.update(d)                 #正链字典

file1 = open(out_circ_file+'circ_exon_uniqer.gff','r')
file2 = open(out_circ_file+'circ_exon_connect.fasta','w+')
#一条circRNA对应外显子的起始和终止位置和外显子数
circID_exon_multiple_seq = ''
for line1 in file1:
    A = line1.strip().split('\t')
    circID = A[0]
    exon_multiple = A[1].strip().split(';')
    exon_count = A[2].split(':')
    circ_RNA_s = A[-1]
    exon_multiple_p_seq,exon_multiple_n_seq,exon_multiple_pn_seq = '','',''
    for exon_single in exon_multiple:
        #print(exon_single)
        if exon_single in e.keys():
            exon_single_p_seq = str(e[exon_single])
            exon_multiple_pn_seq = exon_multiple_pn_seq + exon_single_p_seq
        if circ_RNA_s == '+':
            #postive
            exon_multiple_p_seq = exon_multiple_pn_seq
            #negative
            exon_multiple_n_seq = reverse_complement(exon_multiple_pn_seq)
        elif circ_RNA_s == '-':
            exon_multiple_p_seq = reverse_complement(exon_multiple_pn_seq)
            exon_multiple_n_seq = exon_multiple_pn_seq
            #这一步是以实际的正链或负链做该序列的正链
        if 'N' not in exon_multiple_p_seq:
            exon_id_multiple_p_seq = '>' + 'exon_' + new_circ_ID_dic[circID] + ':' + new_circ_g_dic[circID] + '_P\n' + exon_multiple_p_seq + '\n'
            exon_id_multiple_n_seq = '>' + 'exon_' + new_circ_ID_dic[circID] + ':' + new_circ_g_dic[circID] + '_R\n' + exon_multiple_n_seq + '\n'
    #circID_exon_multiple_seq是circRNA的正负列集合
    file2.write(exon_id_multiple_p_seq + exon_id_multiple_n_seq)
file2.close()

#############################################################################################################

def bsj_seq(circ_type,circ_ID,circ_seq):
    s_s = len(circ_seq)
    if s_s >= 200:
        bsj_s = circ_seq[-100:] + circ_seq[:99]
    elif s_s < 200:
    #global bsj_123456
        if s_s % 2 == 0:
            bsj_s = circ_seq[-s_s//2:] + circ_seq[:s_s//2-1]
        elif s_s % 2 != 0:
            bsj_s = circ_seq[-s_s//2:] + circ_seq[:s_s//2+1]
    circ_RNA_id_1 = circ_type + circ_ID + '\n'
    bsj_id_seq_1 = '>1_P_' + circ_RNA_id_1 + translate(bsj_s) + '\n'
    bsj_id_seq_2 = '>2_P_' + circ_RNA_id_1 + translate(bsj_s[1:]) + '\n'
    bsj_id_seq_3 = '>3_P_' + circ_RNA_id_1 + translate(bsj_s[2:]) + '\n'
    bsj_id_seq_r1 = '>r1_R_' + circ_RNA_id_1 + translate(reverse_complement(bsj_s)) + '\n'
    bsj_id_seq_r2 = '>r2_R_' + circ_RNA_id_1 + translate(reverse_complement(bsj_s)[1:]) + '\n'
    bsj_id_seq_r3 = '>r3_R_' + circ_RNA_id_1 + translate(reverse_complement(bsj_s)[2:]) + '\n'
    bsj_123456 = bsj_id_seq_1 + bsj_id_seq_2 + bsj_id_seq_3 + bsj_id_seq_r1 + bsj_id_seq_r2 + bsj_id_seq_r3
    return bsj_123456

file_circbase = open(out_circ_file+'plant_circ_pr.txt','r')
file_circ_exon = open(out_circ_file+'circ_exon_connect.fasta','r')
list_base_ID,dic_base = [],{}
for line_base in file_circbase:
    if line_base[0] == '>':
        base_ID = line_base.strip('>').strip()
        base_ID_l = base_ID.strip('P').strip('R').strip('_')
        list_base_ID.append(base_ID_l)
    elif line_base[0] != '>':
        base_seq = line_base.strip()
        dic_base.update({base_ID:base_seq})

list_exon_ID,dic_exon = [],{}
for line_exon in file_circ_exon:
    if line_exon[0] == '>':
        exon_ID = line_exon.strip().replace('>exon_','')
        exon_ID_l = exon_ID.strip('P').strip('R').strip('_')
        list_exon_ID.append(exon_ID_l)
    elif line_exon[0] != '>':
        exon_seq = line_exon.strip()
        dic_exon.update({exon_ID:exon_seq})

set_base_exon = set(list_base_ID)&set(list_exon_ID)
set_only_base = set(list_base_ID) - set(list_exon_ID)
set_only_exon = set(list_exon_ID) - set(list_base_ID)

now_circ_marge = open(out_circ_file+'circ_base_add_exon_connect_PR.fasta','w+')
now_circ_bsj_marge = open(out_circ_file+'circ_base_add_exon_connect_bsj_PR.fasta','w+')
multiple_ID_seq,bsj_all = '',''
for circ_id in set_base_exon:
    circ_id_PP = circ_id + '_P'
    circ_id_RR = circ_id + '_R'
    circ_id_id = circ_id.strip().split(':')[0]
    if dic_base[circ_id_PP] == dic_exon[circ_id_PP] or dic_base[circ_id_PP] == dic_exon[circ_id_RR]:
        seq_base_exon = dic_base[circ_id_PP]
        ID_seq = '>P_exon_' + circ_id + '\n' + seq_base_exon + '\n'
        ID_seq_R = '>R_exon_' + circ_id + '\n' + dic_base[circ_id_RR] + '\n'
        #加入BSJ函数
        circ_type = 'exon_'
        bsj_seq_6 = bsj_seq(circ_type,circ_id_id,seq_base_exon)
        now_circ_bsj_marge.write(bsj_seq_6)
        #bsj_all = bsj_all + bsj_seq_6
    elif dic_base[circ_id_PP] != dic_exon[circ_id_PP] and dic_base[circ_id_PP] != dic_exon[circ_id_RR]:
        seq_base = dic_base[circ_id_PP]
        seq_exon = dic_exon[circ_id_PP]
        ID_seq = '>P_base_' + circ_id + '\n' + seq_base + '\n>P_self_' + circ_id + '\n' + seq_exon + '\n'
        ID_seq_R = '>R_base_' + circ_id + '\n' + dic_base[circ_id_RR] + '\n>R_self_' + circ_id + '\n' + dic_exon[circ_id_RR] + '\n'
        circ_type = 'base_'
        bsj_seq_6_base = bsj_seq(circ_type,circ_id_id,seq_base)
        circ_type = 'self_'
        bsj_seq_6_exon = bsj_seq(circ_type,circ_id_id,seq_exon)
        now_circ_bsj_marge.write(bsj_seq_6_base + bsj_seq_6_exon)
        #bsj_all = bsj_all + bsj_seq_6_base + bsj_seq_6_exon
    now_circ_marge.write(ID_seq + ID_seq_R)
    #multiple_ID_seq = multiple_ID_seq + ID_seq + ID_seq_R


for circ_id in set_only_base:
    circ_id_PP = circ_id + '_P'
    circ_id_RR = circ_id + '_R'
    circ_id_id = circ_id.strip().split(':')[0]
    if circ_id not in list_circ_unknown:
        ID_seq = '>P_Intergenic_' + circ_id + '\n' + dic_base[circ_id_PP] + '\n'
        ID_seq_R = '>R_Intergenic_' + circ_id + '\n' + dic_base[circ_id_RR] + '\n'
        circ_type = 'Intergenic_'
    elif circ_id in list_circ_unknown:
        ID_seq = '>P_unknown_' + circ_id + '\n' + dic_base[circ_id_PP] + '\n'
        ID_seq_R = '>R_unknown_' + circ_id + '\n' + dic_base[circ_id_RR] + '\n'
        circ_type = 'unknown_'
    now_circ_marge.write(ID_seq + ID_seq_R)
    bsj_6 = bsj_seq(circ_type,circ_id_id,dic_base[circ_id_PP])
    now_circ_bsj_marge.write(bsj_6)

for circ_id in set_only_exon:
    print(circ_id)

now_circ_marge.close()
now_circ_bsj_marge.close()


def RNA_protein(RNA_string):

    #start_code = 'ATG'
    #end_code = ['UAA', 'UAG', 'UGA']
    protein_table = {'TTT': 'F', 'CTT': 'L', 'ATT': 'I', 'GTT': 'V', \
    'TTC': 'F', 'CTC': 'L', 'ATC': 'I', 'GTC': 'V', \
    'TTA': 'L', 'CTA': 'L', 'ATA': 'I', 'GTA': 'V', \
    'TTG': 'L', 'CTG': 'L', 'ATG': 'M', 'GTG': 'V', \
    'TCT': 'S', 'CCT': 'P', 'ACT': 'T', 'GCT': 'A', \
    'TCC': 'S', 'CCC': 'P', 'ACC': 'T', 'GCC': 'A', \
    'TCA': 'S', 'CCA': 'P', 'ACA': 'T', 'GCA': 'A', \
    'TCG': 'S', 'CCG': 'P', 'ACG': 'T', 'GCG': 'A', \
    'TAT': 'Y', 'CAT': 'H', 'AAT': 'N', 'GAT': 'D', \
    'TAC': 'Y', 'CAC': 'H', 'AAC': 'N', 'GAC': 'D', \
    'TAA': '*', 'CAA': 'Q', 'AAA': 'K', 'GAA': 'E', \
    'TAG': '*', 'CAG': 'Q', 'AAG': 'K', 'GAG': 'E', \
    'TGT': 'C', 'CGT': 'R', 'AGT': 'S', 'GGT': 'G', \
    'TGC': 'C', 'CGC': 'R', 'AGC': 'S', 'GGC': 'G', \
    'TGA': '*', 'CGA': 'R', 'AGA': 'R', 'GGA': 'G', \
    'TGG': 'W', 'CGG': 'R', 'AGG': 'R', 'GGG': 'G'}

    #start_sit = re.search(start_code, RNA_string)
    protein_list=[]
    index_list=[]

    RNA_string2=RNA_string+RNA_string[:3]
    RNA_string3=RNA_string*3
    pm=0
    while(1):
        pm=(RNA_string2).find('ATG')
        #print(pm)
        rp='y'
        if(pm<0):
            break
        for ii in range(0,len(index_list)):
            if((pm+1-index_list[ii])%3==0 and pm+1-index_list[ii]<len(protein_list[ii])*3):
                RNA_string2=RNA_string2[:pm]+'@'+RNA_string2[pm+1:]
                rp='n'
                break
        if(rp=='n'):
            continue
        RNA_string2=RNA_string2[:pm]+'@'+RNA_string2[pm+1:]
        pm=pm+3
        protein='M'
        for sit in range(pm,((len(RNA_string3)-pm)//3)*3+pm, 3):
            protein = protein + protein_table[RNA_string3[sit:sit+3]]
            #print(protein)
            if (protein_table[RNA_string3[sit:sit+3]]=='*'):
                break
        if(len(protein)>20):
            protein_list.append(protein)
            index_list.append(pm-2)

    return protein_list,index_list
    #RNA_string = open('E:\\bioinfo\data\\rosalind_prot\\rosalind_prot.txt', 'r').read().strip()
    #mRNA_protein(RNA_string)

#预测的开放阅读框
def fagain2(file_circ_seq,file_corf_out):
    #os.system('cp %s %s/circRNA_exon.fasta'%(file_circ_seq,file_corf_out))
    fi=open(file_circ_seq,'r+')
    fo=open(file_corf_out+'/circ_draw.txt','w+')
    while(True):
        si=fi.readline()
        if(si==''):
            break
        if(si[0]=='>'):
            fo.write(si[1:].split(':')[0])
            seqn=fi.readline().strip()
            ll=str(len(seqn))
            fo.write(','+ll+',0-'+ll+','+seqn+'\n')
  
    fi.close()
    fo.close()

#---------------------------------------------
#-------------orf predict---------------------
#---------------------------------------------


def orf_predict(file_circ_seq,file_corf_out):
    fi=open(file_circ_seq,'r+')
    fo=open(file_corf_out+'/circRNA_corf.txt','w+')
    print('predict circular RNA orf ,,,')
    circn=''
    seq_exon=''
    for si in fi.readlines():
        if(si[0]=='>'):
            circn=si.strip().split(':')
            #print(circn[0])
            continue
        seq_exon=si.strip().upper()
        if(len(seq_exon)>10000):
            continue
        if 'N' in seq_exon :           #自己加的
            continue
        orf,index=RNA_protein(seq_exon)
        #rint(seq_exon,orf)
        for ii in range(0,len(orf)):
            fo.write(circn[0]+'-ORF%d(%d,%d):%s\n'%(ii+1,index[ii],index[ii]+len(orf[ii])*3,circn[1])+orf[ii]+'\n')   

    fi.close()
    fo.close()

#fagain2('/home/li/public2/wk/work/Oryza_sativa_japonica_circ_mixture/circ_base_add_exon_connect_PR.fasta','./')
#orf_predict('/home/li/public2/wk/work/Oryza_sativa_japonica_circ_mixture/circ_base_add_exon_connect_PR.fasta','./')

def bsj_Peptide_fragment(bsj_peptide_200):
    if len(bsj_peptide_200)%2 == 0:
        bsj_original_position = int(len(bsj_peptide_200)/2)
    elif len(bsj_peptide_200)%2 != 0:
        bsj_original_position = int((len(bsj_peptide_200)-1)/2)

    #print(bsj_original_position)
    if '*' not in bsj_peptide_200 :
        def_bsj_peptide = bsj_peptide_200
        def_bsj_position = bsj_original_position
    elif '*' in bsj_peptide_200 :
        bsj_middle_peptide = bsj_peptide_200[:bsj_original_position] + '|' + bsj_peptide_200[bsj_original_position:]
        list_bsj_peptide = bsj_middle_peptide.split('*')
        bsj_middle_count = 0
        for bsj_peptide_single in list_bsj_peptide:
            bsj_middle_count += 1
            if '|' in bsj_peptide_single and bsj_middle_count == 1:
                def_bsj_peptide = bsj_peptide_single.replace('|','')
                def_bsj_position = bsj_original_position
            elif '|' in bsj_peptide_single and bsj_middle_count != 1:
                bsj_before_position = bsj_peptide_single.find('|')
                if 'M' in bsj_peptide_single[:bsj_before_position] :
                    M_position = bsj_peptide_single.find('M')
                    def_bsj_peptide = bsj_peptide_single[M_position:].replace('|','')
                    def_bsj_position = bsj_peptide_single[M_position:].find('|')
                    break
                elif 'M' not in bsj_peptide_single[:bsj_before_position]:
                    def_bsj_peptide = ''
                    def_bsj_position = ''
                    break

    return def_bsj_peptide,def_bsj_position
        
def find_all_occurrences(string, sub_string):
    positions = []
    start = 0
    while True:
        position = string.find(sub_string, start)
        if position == -1:
            break
        positions.append(position)
        start = position + 1
    return positions

def orf_bsj(file_circ_bsj,file_corf_out):
  #global bsj_corf_position
  dic_id_bsj_seq = {}
  fi_bsj = open(file_circ_bsj,'r')
  for row_bsj in fi_bsj:
    if row_bsj[0] == '>':
      circ_bsj_id = row_bsj.strip().split('_',1)[1]
    elif row_bsj[0] != '>':
      #print(row_bsj.strip())
      bsj_peptide,bsj_position = bsj_Peptide_fragment(row_bsj.strip())
      #print(bsj_peptide,bsj_position)

      if bsj_peptide != '' and bsj_position != '':
        if circ_bsj_id not in dic_id_bsj_seq.keys():
          dic_id_bsj_seq.update({circ_bsj_id:bsj_peptide + ':' + str(bsj_position)})

        elif circ_bsj_id in dic_id_bsj_seq.keys():
          dic_id_bsj_seq_value_old = dic_id_bsj_seq[circ_bsj_id]
          dic_id_bsj_seq.update({circ_bsj_id:dic_id_bsj_seq[circ_bsj_id] + ';' + bsj_peptide + ':' + str(bsj_position)})
  #print(str(dic_id_bsj_seq)[:500])
  fi_corf = open(file_corf_out+'/circRNA_corf.txt','r')
  fo_bsj_corf = open(file_corf_out+'/circRNA_bsj_corf.fasta','w+')
  for row_corf in fi_corf:
    if row_corf[0] == '>':
      row_corf_id = row_corf.strip()
      row_corf_bsj_id = row_corf.strip().strip('>').split('-')[0]
    elif row_corf[0] != '>':
      row_corf_seq = row_corf.strip()
      if row_corf_bsj_id not in dic_id_bsj_seq.keys():
        continue
      list_bsj_peptide_position = dic_id_bsj_seq[row_corf_bsj_id]
      for bsj_peptide_position_single in list_bsj_peptide_position.split(';'):
        bsj_peptide = bsj_peptide_position_single.split(':')[0]
        bsj_position = bsj_peptide_position_single.split(':')[1]
        if bsj_peptide in row_corf_seq:
          list_positions = find_all_occurrences(row_corf_seq,bsj_peptide)
          bsj_corf_position = ''
          for bsj_peptide_position in list_positions:
            bsj_corf_position = bsj_corf_position + str(int(bsj_position) + bsj_peptide_position) + '|'
          row_corf_id_bsj_seq = row_corf_id + ':' + bsj_corf_position.strip('|') + '\n' + row_corf_seq + '\n'
          fo_bsj_corf.write(row_corf_id_bsj_seq)
  fi_bsj.close()
  fi_corf.close()
  fo_bsj_corf.close()

fagain2(out_circ_file+'circ_base_add_exon_connect_PR.fasta',out_circ_file)
orf_predict(out_circ_file+'circ_base_add_exon_connect_PR.fasta',out_circ_file)
orf_bsj(out_circ_file+'circ_base_add_exon_connect_bsj_PR.fasta',out_circ_file)

if 'circ_info_omit.gff' not in os.listdir(out_circ_file):
    os.system('scp %s %scirc_info_omit.gff'%(circ_info_omit_file,out_circ_file))
os.system('mkdir %sMScirctl1_main_result %sMScirctl1_all_results'%(out_circ_file,out_circ_file))
os.system('scp %scirc_base_add_exon_connect_PR.fasta %scirc_info_omit.gff %scircRNA_bsj_corf.fasta %sMScirctl1_main_result'%(out_circ_file,out_circ_file,out_circ_file,out_circ_file))
os.system('mv %s*.* %sMScirctl1_all_results'%(out_circ_file,out_circ_file))
print('mv %s*.* %sMScirctl1_all_results'%(out_circ_file,out_circ_file))
#os.system("mv !\(%sMScirctl1_all_results|%sMScirctl1_main_result\) %s/* %s/MScirctl1_all_results"%(out_circ_file,out_circ_file,out_circ_file,out_circ_file))
