#!/usr/bin/python3

import time
import os
import numpy as np
from collections import Counter
from copy import deepcopy
import argparse
import matplotlib.pyplot as plt
from matplotlib import colors

parser = argparse.ArgumentParser('circRNA second map')
#parser.add_argument('-cs','--circ_seq',help='circRNA序列,fasta格式')
#parser.add_argument('-cio','--circ_info_omit_file',help='circ_corf.py生成的circRNA注释文件')
#parser.add_argument('-co','--circ_orf',help='circRNA的预测到的开放阅读框')
parser.add_argument('-fi','--circ_file_input',help='MScirctl1的主要结果文件,MScirctl1_main_result')
parser.add_argument('-pit','--path_input_tsv',help='MSGFPlus结果tsv格式的文件夹',default='none')
parser.add_argument('-pip','--path_input_peptide',help='其它质谱软件肽段筛选汇总结果',default='none')
parser.add_argument('-fip','--file_input_protein',help='该物种线性蛋白序列,fasta格式,默认为没有该文件',default='none')
parser.add_argument('-fims','--file_input_PRIDE',help='PRIDE数据库中的质谱信息,自己制作,默认为没有该文件',default='none')
parser.add_argument('-ec','--excellent_count',help='MScirctl2可视化展示最有可能翻译的circRNA,默认为10',default='10')
parser.add_argument('-fo','--file_out',help='输出结果文件夹',default='./')

args = parser.parse_args()
#circ_seq = args.circ_seq
#circ_info_omit_file = args.circ_info_omit_file
#circ_orf = args.circ_orf
circ_file_input = args.circ_file_input
path_input_tsv = args.path_input_tsv
path_input_peptide = args.path_input_peptide
file_input_protein = args.file_input_protein
file_input_PRIDE = args.file_input_PRIDE
excellent_count = args.excellent_count
file_out = args.file_out
excellent_count = int(excellent_count)

if circ_file_input[-1] != '/':
    circ_file_input = circ_file_input + '/'
circ_seq = circ_file_input + 'circ_base_add_exon_connect_PR.fasta'
circ_info_omit_file = circ_file_input + 'circ_info_omit.gff'
circ_orf = circ_file_input + 'circRNA_bsj_corf.fasta'

def Create_folder():
    now=time.strftime('%y-%m-%d %H:%M:%S',time.localtime())
    ii=1
    file_out='MScirctl2_'+now.split()[0]
    while(os.path.exists(file_out+'.%d'%ii)):
        ii=ii+1
    file_out=file_out+'.%d'%ii + '/'
    os.system('mkdir %s'%(file_out))
    return file_out

if file_out == './':
    file_out = Create_folder()
elif file_out != './' and file_out[-1] != '/':
    file_out = file_out + '/'
    
def circ_id_seq():
    dic_circ_id_seq = {}
    file_circ_seq = open(circ_seq,'r')
    for row_circ_seq in file_circ_seq:
        if row_circ_seq[0] == '>':
            circ_name=row_circ_seq[1:].split(':')[0]
        elif row_circ_seq[0] != '>':
            circ_ires_seq = row_circ_seq.strip()[0:300]
            dic_circ_id_seq.update({circ_name:circ_ires_seq})
    return dic_circ_id_seq

def ires_predict(circ_id):
  #list_ires_circ_id = []
  #fj=open(circ_single_seq,'r+')
  fp=open(file_out+'circRNA_ires.fasta','w+')
  fp.write('>' + circ_id + '\n' + dic_circ_id_seq[circ_id][0:300] + '\n')
  fp.close()

  os.system('python3 IRESfinder-master/IRESfinder.py -f '+ file_out + 'circRNA_ires.fasta -o ' + file_out + 'circRNA_ires.res -m 2 -w 50 -s 50 -p '+file_out)
  #print('python3 IRESfinder-master/IRESfinder.py -f '+ file_out + 'circRNA_ires.fasta -o ' + file_out + 'circRNA_ires.res -m 2 -w 50 -s 50 -p '+file_out)
  if os.path.exists(file_out + 'circRNA_ires.res'):
    fk=open(file_out + 'circRNA_ires.res','r+')
    for sk in fk.readlines():
      if ('IRES' not in sk):
        continue
      skt = sk.split('\t')
      circ_id = skt[0].rsplit('_',2)[0]
      circ_ires = skt[1]
      #print(circ_ires)
      if circ_ires == 'IRES':
        #list_ires_circ_id.append(circ_id)
        circ_def_ires = 'ires'
        break
      circ_def_ires = 'non-ires'
    #list_ires_circ_id = list(set(list_ires_circ_id))
    print('-----------IRES analyse finished!!----------------')
    #fj.close()
    fk.close()
    os.system('rm ' +file_out+'circRNA_ires.res')
  else:
    circ_def_ires = 'non-ires'
  os.system('rm ' +file_out+'circRNA_ires.fasta')
  return circ_def_ires
    
def tsv_circ_protaid():
    QValue = ''
    dic_id_peptide = {}
    file_line_circ = open(file_out+"tsv_result.txt", "w+")
    file_line_circ.write('MS_ID\tpeptide\tcORF_ID\tpeptide_score\n')
    for root,dirs,files in os.walk(path_input_tsv):
        if dirs == []:
            for file_tsv in files:
                path_file_tsv = root + '/' + file_tsv
                content_tsv = open(path_file_tsv,'r')
                line_hang_1 = content_tsv.readlines(1)
                for line_content_tsv in content_tsv:
                    file_id = line_content_tsv.strip().split('\t')[0].replace('tsv','mgf')
                    Peptide = line_content_tsv.strip().split('\t')[9].strip().split('.')[1]
                    circ_id_multiple = line_content_tsv.strip().split('\t')[10].strip().split(', ')
                    MSGFScore = float(line_content_tsv.strip().split('\t')[12])
                    EValue = line_content_tsv.strip().split('\t')[14]
                    if len(line_content_tsv.strip().split('\t')) > 16:
                        QValue = line_content_tsv.strip().split('\t')[15]
                        if float(QValue) <= 0.01:
                            for circ_id in circ_id_multiple:
                                line_multiple = file_id + '\t' + Peptide + '\t' + circ_id.strip() + '\t' + str(MSGFScore) + '\n'
                                file_line_circ.write(line_multiple)
                                if file_input_protein == 'none':
                                    dic_id_peptide = dic_circ_mgf_pep_score(circ_id.strip(),file_id,Peptide,str(MSGFScore),dic_id_peptide)
                    elif len(line_content_tsv.strip().split('\t')) <= 16:
                        if float(EValue) <= 0.01:
                            for circ_id in circ_id_multiple:
                                line_multiple = file_id + '\t' + Peptide + '\t' + circ_id.strip() + '\t' + str(MSGFScore) + '\n'
                                file_line_circ.write(line_multiple)
                                if file_input_protein == 'none':
                                    dic_id_peptide = dic_circ_mgf_pep_score(circ_id.strip(),file_id,Peptide,str(MSGFScore),dic_id_peptide)
    file_line_circ.close()
    if QValue != '':
        print('Screening condition is QValue <= 0.01')
    elif QValue == '':
        print('Screening condition is EValue <= 0.01')
    file_line_circ.close()
    print('---------------------tsv_circ_protaid---------------------')
    if file_input_protein == 'none':
        return dic_id_peptide

def dic_circ_mgf_pep_score(circ_id,mgf_id,peptide,MSGF_score,dic_id_peptide):
    if circ_id in dic_id_peptide.keys():
        dic_id_peptide[circ_id] = dic_id_peptide[circ_id] + ';' + mgf_id + '/' + peptide + '/' + MSGF_score 
    if circ_id not in dic_id_peptide.keys():
        dic_id_peptide.update({circ_id:mgf_id + '/' + peptide + '/' + MSGF_score})
    return dic_id_peptide

def tsv_remove_line_protaid():
    file_line_protide = open(file_input_protein,'r')
    line_protide_multiple = ''
    list_line_protide = []
    for line_protide in file_line_protide:
        if line_protide[0] == '>':
            list_line_protide.append(line_protide_multiple)
            line_protide_multiple = ''
        elif line_protide[0] != '>':
            line_protide_single = line_protide.strip()
            line_protide_multiple = line_protide_multiple + line_protide_single
    list_line_protide.append(line_protide_multiple)

    if path_input_tsv != 'none' and path_input_peptide == 'none':
        file_tsv_result = open(file_out+"tsv_result.txt",'r')
    elif path_input_peptide != 'none' and path_input_tsv == 'none':
        file_tsv_result = open(path_input_peptide,'r')

    file_tsv_result_header = file_tsv_result.readline()
    file_remove = open(file_out+'tsv_result_remove_line.txt','w+')
    file_remove.write(file_tsv_result_header)
    list_line_p = []
    list_circ_p = []
    #dic_id_peptide = dic_id_result
    dic_id_peptide = {}
    for line_tsv_result in file_tsv_result:
        #print(line_tsv_result)
        mgf_id = line_tsv_result.strip().split('\t')[0]
        peptide = line_tsv_result.strip().split('\t')[1]
        circ_id = line_tsv_result.strip().split('\t')[2]
        MSGF_score = line_tsv_result.strip().split('\t')[3]
        peptide_type = 'c'
        if peptide in list_line_p:
            peptide_type = 'l'
            continue
        if peptide in list_circ_p:
            peptide_type = 'c'
            file_remove.write(line_tsv_result)
            dic_id_peptide = dic_circ_mgf_pep_score(circ_id,mgf_id,peptide,MSGF_score,dic_id_peptide)
            continue
        for ath_line_protide in list_line_protide:
            if peptide in ath_line_protide :
                list_line_p.append(peptide)
                peptide_type = 'l'
                break
        if peptide_type == 'c':
            list_circ_p.append(peptide)
            file_remove.write(line_tsv_result)
            dic_id_peptide = dic_circ_mgf_pep_score(circ_id,mgf_id,peptide,MSGF_score,dic_id_peptide)
    print(len(list_line_p),'线性')
    print(len(list_circ_p),'环状')
    file_remove.close()
    print('------------------tsv_remove_line_protaid------------------')
    return dic_id_peptide

def func(lst):       #将不连序的列表变成几个连续的列表
	new_list = []
	for i,j in zip(lst,lst[1:]):
		if j - i > 1:
			new_list.append(lst[:lst.index(j)])
			lst = lst[lst.index(j):]
	new_list.append(lst)
	return new_list

def PRIDE_experiment():
    file_organization = open(file_input_PRIDE,'r')
    dic_raw_organization = {}
    dic_experiment = {}
    list_control_experiment_raw_id = []
    dic_PXD_control_raw,dic_PXD_experiment_raw = {},{}
    for line_organization in file_organization:
        PXD_id = line_organization.strip().split('\t')[0]
        raw_id = line_organization.strip().split('\t')[1].replace('.raw','.mgf').replace('.RAW','.mgf')
        if '20181214_LFQ' in raw_id:
            raw_id = raw_id.replace('+','%2B')
        if 'Q35D185' in raw_id:
            raw_id = raw_id.replace('%28','(').replace('%29',')')
        organization = line_organization.strip().split('\t')[5]
        plant_experiment = line_organization.strip().split('\t')[6]
        experiment_class = line_organization.strip().split('\t')[7]
        dic_raw_organization.update({raw_id:organization})
        if plant_experiment == '':
            plant_experiment = 'none'
        dic_experiment.update({PXD_id:plant_experiment})
        if experiment_class == 'control' or experiment_class == 'Control': 
            if PXD_id not in dic_PXD_control_raw.keys():
                dic_PXD_control_raw.update({PXD_id:raw_id})
            elif PXD_id in dic_PXD_control_raw.keys():
                dic_PXD_control_raw_value = dic_PXD_control_raw[PXD_id] + ';' + raw_id
                dic_PXD_control_raw.update({PXD_id:dic_PXD_control_raw_value})
            list_control_experiment_raw_id.append(raw_id)
        elif experiment_class == 'experiment' or experiment_class == 'Experiment':
            if PXD_id not in dic_PXD_experiment_raw.keys():
                dic_PXD_experiment_raw.update({PXD_id:raw_id})
            elif PXD_id in dic_PXD_experiment_raw.keys():
                dic_PXD_experiment_raw_value = dic_PXD_experiment_raw[PXD_id] + ';' + raw_id
                dic_PXD_experiment_raw.update({PXD_id:dic_PXD_experiment_raw_value})
            list_control_experiment_raw_id.append(raw_id)
    print('---------------PRIDE_experiment----------------')
    return dic_experiment,dic_raw_organization,list_control_experiment_raw_id,dic_PXD_control_raw,dic_PXD_experiment_raw

def mgf_amount_peptide():
    file_remove_line = open(file_out+'tsv_result_remove_line.txt','r')
    dic_mgf_circ_amount = {}
    for row_remove_line in file_remove_line:
        mgf_id = row_remove_line.strip().split('\t')[0]
        circ_id = row_remove_line.strip().split('\t')[2].strip().split('(')[0]
        if mgf_id in list_control_experiment_raw_id:
            if mgf_id not in dic_mgf_circ_amount.keys():
                dic_small = {circ_id:'1'}
                dic_mgf_circ_amount.update({mgf_id:dic_small})
            elif mgf_id in dic_mgf_circ_amount.keys():
                if circ_id not in dic_small.keys():
                    dic_small.update({circ_id:'1'})
                elif circ_id not in dic_small.keys():
                    dic_small_value = str(int(dic_small[circ_id]) + 1)
                    dic_small.update({circ_id:dic_small_value})
                dic_mgf_circ_amount.update({mgf_id:dic_small})
    print('-------------mgf_amount_peptide---------------')    
    return dic_mgf_circ_amount

def merge_dictionary(list_control_experiment_circ,dic_small,dic_big):
    for circ_id in dic_small.keys():
        if circ_id not in list_control_experiment_circ:
            list_control_experiment_circ.append(circ_id)
        amount = dic_small[circ_id]
        if circ_id not in dic_big:
            dic_big.update({circ_id:amount})
        elif circ_id in dic_big:
            amount_old = dic_big[circ_id]
            dic_big.update({circ_id:str(int(amount_old)+int(amount))})
    return list_control_experiment_circ,dic_big

def PRIDE_experiment_class():
    path_experiment = file_out+'PRIDE_file_experiment_class/'
    path_png = file_out+'PRIDE_file_experiment_class_png/'
    path_edge = file_out+'PRIDE_file_experiment_class_edge/'
    os.system('mkdir %s %s %s'%(path_experiment,path_png,path_edge))

    for PXD_id in dic_PXD_control_raw.keys() :
        if PXD_id not in dic_PXD_experiment_raw.keys():
            continue
        file_experiment = open(path_experiment+PXD_id,'w+')
        head = 'circ_ORF\tcontrol\texperiment\n'
        file_experiment.write(head)
        list_control_mgf_id = dic_PXD_control_raw[PXD_id].split(';')
        list_experiment_mgf_id = dic_PXD_experiment_raw[PXD_id].split(';')
        dic_control_PXD_circ_amount = {}
        dic_experiment_PXD_circ_amount = {}
        list_control_experiment_circ = []
        for mgf_control_id in list_control_mgf_id:
            if mgf_control_id in dic_mgf_circ_amount.keys():
                dic_control_circ_amount = dic_mgf_circ_amount[mgf_control_id]
                list_control_experiment_circ,dic_control_PXD_circ_amount = merge_dictionary(list_control_experiment_circ,dic_control_circ_amount,dic_control_PXD_circ_amount)
            elif mgf_control_id not in dic_mgf_circ_amount.keys():
                print(mgf_control_id,'不在对照组中')
                continue
                #dic_control_circ_amount = dic_mgf_circ_amount[mgf_control_id]
                #list_control_experiment_circ,dic_control_PXD_circ_amount = merge_dictionary(list_control_experiment_circ,dic_control_circ_amount,dic_control_PXD_circ_amount)
                #for circ_id in dic_control_circ_amount.keys():
                #    amount = dic_control_circ_amount[circ_id]
                #    if circ_id not in dic_control_PXD_circ_amount:
                #        dic_control_PXD_circ_amount.update({circ_id:amount})
                #    elif circ_id in dic_control_PXD_circ_amount:
                #        amount_old = dic_control_PXD_circ_amount[circ_id]
                #        dic_control_PXD_circ_amount.update({circ_id:str(int(amount_old)+int(amount))})
        for mgf_experiment_id in list_experiment_mgf_id:
            if mgf_experiment_id in dic_mgf_circ_amount.keys():
                dic_experiment_circ_amount = dic_mgf_circ_amount[mgf_experiment_id]
                list_control_experiment_circ,dic_experiment_PXD_circ_amount = merge_dictionary(list_control_experiment_circ,dic_experiment_circ_amount,dic_experiment_PXD_circ_amount)
            elif mgf_experiment_id not in dic_mgf_circ_amount.keys():
                print(mgf_experiment_id,'不在实验组中')
                continue
        for circ_id in list_control_experiment_circ :
            if circ_id in dic_control_PXD_circ_amount.keys():
                control_amount = dic_control_PXD_circ_amount[circ_id]
            elif circ_id not in dic_control_PXD_circ_amount.keys():
                control_amount = '0'
            if circ_id in dic_experiment_PXD_circ_amount.keys():
                experiment_amount = dic_experiment_PXD_circ_amount[circ_id]
            elif circ_id not in dic_experiment_PXD_circ_amount.keys():
                experiment_amount = '0'
            row_experiment_class = circ_id + '\t' + control_amount + '\t' + experiment_amount + '\n'
            file_experiment.write(row_experiment_class)
        file_experiment.close()
    #print(str(dic_control_PXD_circ_amount)[0:200])
    #print(str(dic_experiment_PXD_circ_amount)[0:200])
    #print(str(dic_experiment_PXD_circ_amount)[0:200])
    for file_experiment in os.listdir(path_experiment):
        path_file_experiment = path_experiment + file_experiment
        path_file_png = path_png + file_experiment + '.png'
        path_file_edge = path_edge + file_experiment
        #print('Rscript degeR_single_for.R %s %s %s '%(path_file_experiment,path_file_png,path_file_edge))
        #os.system('Rscript /home/li/public2/wk/procedure/R_msgf/degeR_single_for.R %s %s %s '%(path_file_experiment,path_file_png,path_file_edge))
        try :
            os.system('Rscript degeR_single_for.R %s %s %s '%(path_file_experiment,path_file_png,path_file_edge))
        except :
            print('It may be that there are too few samples of differential expression')
    print('-----------------------PRIDE_experiment_class---------------------')

def edge_predict():
    dic_circ_PXD_edge = {}
    path_edge = file_out+'PRIDE_file_experiment_class_edge/'
    for file_id_edge in os.listdir(path_edge):
        file_edge = open(path_edge+file_id_edge,'r')
        for line_edge in file_edge:
            if 'TRUE' in line_edge:
                if float(line_edge.strip().split('\t')[1]) > 0 :
                    diff = 'up'
                elif float(line_edge.strip().split('\t')[1]) < 0 :
                    diff = 'down'
                #list_file_id.append(file_id_edge)
                corf_id = line_edge.strip().split('\t')[0]
                FC = line_edge.strip().split('\t')[1]
                experiment = dic_experiment[file_id_edge]     ####
                if corf_id not in dic_circ_PXD_edge.keys():
                    value = file_id_edge + ':' + experiment + '\t' + FC + '\t' + diff
                    dic_circ_PXD_edge.update({corf_id:value})
                elif corf_id in dic_circ_PXD_edge.keys():
                    value_old = dic_circ_PXD_edge[corf_id]
                    value = value_old.strip().split('\t')[0] + ';' + file_id_edge + ':' + experiment + '\t' + value_old.strip().split('\t')[1] + ';' + FC + '\t' + value_old.strip().split('\t')[2] + ';' + diff
                    dic_circ_PXD_edge.update({corf_id:value})
    return dic_circ_PXD_edge

def circ_crf_linepro_map():
    os.system('makeblastdb -in %s -dbtype prot -out %sprotein_db.faa -parse_seqids'%(file_input_protein,file_out))
    os.system('blastp -query %s -db %sprotein_db.faa -out %scorf_8_result.txt -outfmt 6 -num_threads 8 -word_size 4 -evalue 0.05'%(circ_orf,file_out,file_out))
    file_corf_result = open(file_out+'corf_8_result.txt','r')
    dic_circ_len_line_lacation,list_corc_id= {},[]
    for line_corf_result in file_corf_result:
        circ_corf_id = line_corf_result.strip().split('\t')[0]
        circ_corf_line_tatio = line_corf_result.strip().split('\t')[2]
        circ_corf_start = line_corf_result.strip().split('\t')[6]
        circ_corf_end = line_corf_result.strip().split('\t')[7]
        circ_corf_line_len = line_corf_result.strip().split('\t')[3]
        if float(circ_corf_line_tatio) >= 95 and int(circ_corf_line_len) >= 20:
            if circ_corf_id not in list_corc_id:
                list_corc_id.append(circ_corf_id)
                dic_circ_len_line_lacation_value = circ_corf_start + ',' + circ_corf_end + ':'
            elif circ_corf_id in list_corc_id:
                list_dic_value = dic_circ_len_line_lacation[circ_corf_id].split(':')
                dic_value_xin = circ_corf_start + ',' + circ_corf_end
                if dic_value_xin not in list_dic_value:
                    dic_circ_len_line_lacation_value = dic_circ_len_line_lacation[circ_corf_id] + dic_value_xin + ':'
                elif dic_value_xin in list_dic_value:
                    continue
            dic_circ_len_line_lacation.update({circ_corf_id:dic_circ_len_line_lacation_value})
    return dic_circ_len_line_lacation
    print('-----------------circ_crf_linepro_map----------------')

def circ_forward_reverse():
    dic_circ_strand = {}
    file_circ_info_omit = open(circ_info_omit_file,'r')
    for row_circ_info_omit in file_circ_info_omit:
        circ_id = row_circ_info_omit.strip().split('\t')[0]
        circ_strand = row_circ_info_omit.strip().split('\t')[5]
        dic_circ_strand.update({circ_id:circ_strand})
    return dic_circ_strand

def circ_orf_bsj_location():
    dic_circ_orf_bsj_location = {}
    flie_circ_orf = open(circ_orf,'r')
    for row_circ_orf in flie_circ_orf:
        if row_circ_orf[0] == '>':
            circ_orf_id = row_circ_orf.strip().strip('>').split('(')[0]
            bsj_location = row_circ_orf.strip().split(':')[2]
            if circ_orf_id not in dic_circ_orf_bsj_location.keys():
                dic_circ_orf_bsj_location.update({circ_orf_id:bsj_location})
            elif circ_orf_id in dic_circ_orf_bsj_location.keys():
                dic_value = dic_circ_orf_bsj_location[circ_orf_id]
                dic_circ_orf_bsj_location.update({circ_orf_id:dic_value + '|' + bsj_location})
    return dic_circ_orf_bsj_location

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

def circ_orf_tsv_peptide_map(circ_orf_id,corf):
    #print(str(dic_id_peptide)[:200])
    if circ_orf_id not in dic_id_peptide.keys():
        #print(circ_orf_id,'没有肽段匹配')
        return 'none','none','none','none','none','none','none'
    mgf_pep_score_assemble = dic_id_peptide[circ_orf_id]
    dic_pro_sore = {}
    dic_mgf_score = {}
    list_mgf_id = []
    list_peptide_start_end_multiple = []
    list_organization = []
    mgf_score_all = ''
    corf_s = len(corf.strip('*'))
    bsj_location = dic_circ_orf_bsj_location[circ_orf_id.split('(')[0]]
    for mgf_pep_score_single in mgf_pep_score_assemble.strip(';').split(';'):
        mgf_id = mgf_pep_score_single.strip().split('/')[0]
        if '20181214_LFQ' in mgf_id:
            mgf_id = mgf_id.replace('+','%2B')
        if 'Q35D185' in mgf_id:
            mgf_id = mgf_id.replace('%28','(').replace('%29',')')
        list_mgf_id.append(mgf_id)    #msf名字的列表
        pro = mgf_pep_score_single.strip().split('/')[1]
        score = float(mgf_pep_score_single.strip().split('/')[2].replace('-',''))
        """
        if mgf_id in dic_raw_organization.keys():
            organization = dic_raw_organization[mgf_id]
        elif mgf_id not in dic_raw_organization.keys():
            organization = 'none'
        list_organization.append(organization)
        dic_organization_number = Counter(list_organization)
        str_dic_organization_number = str(dic_organization_number).strip('Counter({').strip('})').replace("'",'').replace(', ',',').replace(': ',':')
        """
        #dic_mgf_score是mgf文件总的肽段打分字典
        if mgf_id not in dic_mgf_score.keys():
            dic_mgf_score.update({mgf_id:int(score)})
        elif mgf_id in dic_mgf_score.keys():
            dic_mgf_score_value = dic_mgf_score[mgf_id] + int(score)
            dic_mgf_score.update({mgf_id:dic_mgf_score_value})
        
        #pro_number_score短肽出现的数目和它的打分平均值
        value = []
        if pro not in dic_pro_sore.keys():
            value.append(score)
            dic_pro_sore.update({pro:value})
        elif pro in dic_pro_sore.keys():
            value = dic_pro_sore[pro]
            value.append(score)
            dic_pro_sore.update({pro:value})
            
        #新加的将corf断点位置的肽段和不是断点位置的肽段整合进corf中
        if pro in corf:
            positions = find_all_occurrences(corf, pro)
            for start in positions:
                #start = corf.find(pro)
                end = start + len(pro) + int(1)
                list_peptide_start_end = range(start,end)
                list_peptide_start_end_multiple.extend(list_peptide_start_end)
        pro_number_score = ''
        number_all = 0
        for pro in dic_pro_sore.keys():
            #print(dic_pro_sore[pro].strip().split(','))
            score = np.mean(dic_pro_sore[pro])
            number = len(dic_pro_sore[pro])
            number_all = number_all + number
            pro_number_score = pro_number_score + pro + '|' + str(number) + '|' + str('%.3f'%score) + ';'
    #mgf_score_all每个质谱给出该肽段的打分值
    for mgf_id in dic_mgf_score.keys():
        mgf_score_all = mgf_score_all + mgf_id + ':' + str(dic_mgf_score[mgf_id]) + ';'

    raw_sum = len(set(list_mgf_id))
    list_set_peptide_start_end_multiple = list(set(list_peptide_start_end_multiple))
    corf_peptide_vatio = (len(list_set_peptide_start_end_multiple)-1)/corf_s
    list_set_peptide_start_end_multiple.sort()
    new_list = func(list_set_peptide_start_end_multiple)
    if new_list == [] or new_list[0] == []:
        print(circ_id,'\t',corf,'\t',pro,'没有完全匹配的肽端')
    elif new_list != []:
        merge_end = 0
        marge_peptide_multiple = ''
        merge_x_start = ''
        for list_p in new_list:
            merge_start = min(list_p)
            if merge_start > 0 :
                merge_x_start = 'X'*(merge_start - merge_end)
                #用来补充circRNA中没有被比对到的序列,用X代替
                marge_peptide_multiple = marge_peptide_multiple + merge_x_start
            merge_end = max(list_p)
            marge_peptide = corf[merge_start:merge_end]
            marge_peptide_multiple = marge_peptide_multiple + marge_peptide
            #marge_peptide_multiple 是有X和被比对到的序列的一条序列
        if merge_end < corf_s:
            merge_x_end = 'X'*(corf_s - merge_end)
            marge_peptide_multiple = marge_peptide_multiple + merge_x_end
            #增加了尾部是否被比对到的序列
    if '*' in corf :
        marge_peptide_multiple = marge_peptide_multiple + '*'
        
    list_bsj_location_single = []
    if '|' in bsj_location:
        for bsj_location_single in bsj_location.split('|'):
            list_bsj_location_single.append(int(bsj_location_single))
        list_bsj_location_single.sort()
        if bsj_location.count('|') == 1:
            #print(marge_peptide_multiple)
            #print(list_bsj_location_single)
            #print(circ_orf_id)
            #print('++++++++++++++++++++++')
            marge_peptide_bsj_multiple = marge_peptide_multiple[:list_bsj_location_single[0]] + '|' + marge_peptide_multiple[list_bsj_location_single[0]:list_bsj_location_single[1]] + '|' + marge_peptide_multiple[list_bsj_location_single[1]:]
        elif bsj_location.count('|') == 2:
            marge_peptide_bsj_multiple = marge_peptide_multiple[:list_bsj_location_single[0]] + '|' + marge_peptide_multiple[list_bsj_location_single[0]:list_bsj_location_single[1]] + '|' + marge_peptide_multiple[list_bsj_location_single[1]:list_bsj_location_single[2]] + '|' + marge_peptide_multiple[list_bsj_location_single[2]:]
            print(circ_orf_id,'bsj位点有两个')
        elif bsj_location.count('|') >= 2:
            print(circ_orf_id,'bsj位点超过两个,需要检查')
            return 'none','none','none','none','none','none','none'
    elif '|' not in bsj_location:
        marge_peptide_bsj_multiple = marge_peptide_multiple[:int(bsj_location)] + '|' + marge_peptide_multiple[int(bsj_location):]
    #print(marge_peptide_bsj_multiple)
    #print(corf_peptide_vatio)
    #print(bsj_location)
    #print(raw_sum)
    #print(number_all)
    #print(mgf_score_all)

    return pro_number_score,marge_peptide_bsj_multiple,corf_peptide_vatio,bsj_location,raw_sum,number_all,mgf_score_all
    

def circ_peptide_result():
    list_circ_orf_id_function = []
    list_pro_long_score = []
    file_first_predict = ''
    file_circ_orf = open(circ_orf,'r')
    file_circ_predict = open(file_out + 'circ_predict.txt','w+')
    file_circ_predict.write('circRNA_ID\tgene\tstrand\tcirc_ORF_ID\tires\tcORF_seq\tcirc_line_position\tcORF_map\tpeptide_coverage\tcORF_bsj_position\tMS_count\tpeptide_count\tpeptide_count_score\tMS_score\tcORF_score\n')
    for row_circ_orf in file_circ_orf:
        if row_circ_orf[0] == '>':
            circ_orf_id = row_circ_orf.strip().strip('>')
            circ_gene = row_circ_orf.strip().strip('>').split(':')[1]
        elif row_circ_orf[0] != '>':
            #控制函数circ_orf_tsv_peptide_map每个circRNA开放阅读框只运行一次
            if circ_orf_id.split('(')[0] in list_circ_orf_id_function:
                continue
            list_circ_orf_id_function.append(circ_orf_id.split('(')[0])
            #circ_orf_seq是circRNA的开放阅读框
            circ_orf_seq = row_circ_orf.strip()
            #circ_line_location:circRNA与线性蛋白序列一致的蛋白序列位置
            if circ_orf_id in dic_circ_len_line_lacation.keys():
                circ_line_location = dic_circ_len_line_lacation[circ_orf_id]
            elif circ_orf_id not in dic_circ_len_line_lacation.keys():
                circ_line_location = 'none'
            #总结果的写出
            #print(circ_orf_id,circ_orf_seq)
            pro_number_score,marge_peptide_bsj_multiple,corf_peptide_vatio,bsj_location,raw_sum,number_all,mgf_score_all = circ_orf_tsv_peptide_map(circ_orf_id,circ_orf_seq)
            
            #新加字符串circRNA翻译产物打分值
            pro_long_score = 0
            for pro_number_score_single in pro_number_score.strip(';').split(';'):
                if pro_number_score_single == 'none':
                    continue
                #print(pro_number_score_single)
                pro_count = pro_number_score_single.strip().split('|')[1]
                pro_score = pro_number_score_single.strip().split('|')[2]
                pro_long_score = pro_long_score + float(pro_count)*float(pro_score)
            list_pro_long_score.append(pro_long_score)

            if marge_peptide_bsj_multiple == 'none' and corf_peptide_vatio == 'none' :
                continue
            if 'X|' in marge_peptide_bsj_multiple or '|X' in marge_peptide_bsj_multiple:
                marge_peptide_bsj_multiple_replace = marge_peptide_bsj_multiple.replace('|X','X').replace('X|','X')
                if '|' not in marge_peptide_bsj_multiple_replace:
                    print(circ_orf_id,'匹配肽段未穿过circRNA断点位置')
                    continue
            #预测该circRNA是否包含ires元件
            circ_ires = ires_predict(circ_orf_id.strip().split('-')[0])
            #circ_ires = 'none_predict'
            if circ_orf_id in dic_circ_len_line_lacation.keys():
                circ_orf_line_location = dic_circ_len_line_lacation[circ_orf_id]
            elif circ_orf_id not in dic_circ_len_line_lacation.keys():
                circ_orf_line_location = 'none'
                
            circ_info_id = circ_orf_id.strip().split('-')[0].split('_',2)[2]
            if circ_info_id in dic_circ_strand:
                circ_strand = dic_circ_strand[circ_info_id]
            elif circ_info_id not in dic_circ_strand:
                print(circ_orf_id,'不在info_omit文件中')
            
            if circ_strand == '+' and circ_orf_id[0] == 'P' :
                ll = 1
            elif circ_strand == '-' and circ_orf_id[0] == 'R' and circ_orf_id.split('_')[1] != 'self' :
                ll = 1
            elif circ_strand == '-' and circ_orf_id[0] == 'P' and circ_orf_id.split('_')[1] == 'self' :
                ll = 1
            elif circ_strand == 'NA':
                ll = 1
            else :
                continue
            row_circ_predict = circ_info_id + '\t' + circ_gene + '\t' + circ_strand + '\t' + circ_orf_id.strip().split('(')[0] + '\t' + circ_ires + '\t' + circ_orf_seq + '\t' + circ_orf_line_location + '\t' + marge_peptide_bsj_multiple + '\t' + str(corf_peptide_vatio) + '\t' + bsj_location + '\t' + str(raw_sum) + '\t' + str(number_all) + '\t' + pro_number_score + '\t' + mgf_score_all + '\t' + str(pro_long_score) + '\n'
            file_circ_predict.write(row_circ_predict)
            file_first_predict = file_first_predict + row_circ_predict
    file_circ_predict.close()
    list_pro_long_score.sort(reverse=True)
    return list_pro_long_score
    
def circ_predict_PRIDE():
    file_circ_predict = open(file_out + 'circ_predict.txt','r')
    file_circ_predict_header = file_circ_predict.readline()
    file_circ_predict_PRIDE = open(file_out + 'circ_predict_PRIDE.txt','w+')
    file_circ_predict_PRIDE.write(file_circ_predict_header.strip('\n') + '\tpeptide_tissue\tdifference_MS_ID\tlogFC\texpression_level\n')
    for row_circ_predict_PRIDE in file_circ_predict:
        list_organization = []
        #print(row_circ_predict_PRIDE)
        circ_orf_id = row_circ_predict_PRIDE.strip().split('\t')[3]
        list_mgf_score = row_circ_predict_PRIDE.strip().split('\t')[13].strip(';').split(';')
        for mgf_score in list_mgf_score:
            mgf_id = mgf_score.strip().split(':')[0]
            if mgf_id in dic_raw_organization.keys():
                organization = dic_raw_organization[mgf_id]
            elif mgf_id not in dic_raw_organization.keys():
                organization = 'none'
            list_organization.append(organization)
            
        dic_organization_number = Counter(list_organization)
        str_dic_organization_number = str(dic_organization_number).strip('Counter({').strip('})').replace("'",'').replace(', ',',').replace(': ',':')
        
        if circ_orf_id not in dic_circ_PXD_edge.keys():
            PRIDE_class = 'none\tnone\tnone'
        elif circ_orf_id in dic_circ_PXD_edge.keys():
            PRIDE_class = dic_circ_PXD_edge[circ_orf_id]
        row_circ_predict_PRIDE = row_circ_predict_PRIDE.strip() + '\t' + str_dic_organization_number + '\t' + PRIDE_class + '\n'
        file_circ_predict_PRIDE.write(row_circ_predict_PRIDE)
    file_circ_predict.close()
    file_circ_predict_PRIDE.close()

def draw_msgf_second(file_input,excellent_count):
    file_input_content = open(file_input,'r')
    file_input_content_header = file_input_content.readline()
    list_start_end = ['(',')','|']
    if not os.path.exists(file_out+'Peptide_coverage_map'):
        os.system("mkdir %sPeptide_coverage_map/"%(file_out))
    for line_input in file_input_content:
        list_circ_line_multiple = []
        list_protaid_multiple = []
        dic_protaid_score = {}
        corf_id = line_input.strip().split('\t')[0]
        circ_line_position = line_input.strip().split('\t')[6]
        if circ_line_position != 'none':
            for list_line_section in circ_line_position.strip(':').split(':'):
                list_circ_line_position_single = list(range(int(list_line_section.split(',')[0])-1,int(list_line_section.split(',')[1])-1))
                if list_circ_line_position_single not in list_circ_line_multiple:
                    list_circ_line_multiple.append(list_circ_line_position_single)
        circ_gene = line_input.strip().split('\t')[1]
        corf_seq = line_input.strip().split('\t')[5]
        corf_seq_len = len(corf_seq)
        circ_bsj_position = line_input.strip().split('\t')[9].split('|')
        circ_protaid = line_input.strip().split('\t')[12]
        pro_long_score = line_input.strip().split('\t')[14]
        float_circ_seq_coverage = float(line_input.strip().split('\t')[8])*100
        #if float(line_input.strip().split('\t')[8]) < float(pro_cover_limit):
        #    continue
        if float(pro_long_score) not in list_pro_long_score[:excellent_count]:
            continue
        #float_4_circ_seq_coverage = '%.2f'%float_circ_seq_coverage
        #print(float_circ_seq_coverage)
        #print(float_4_circ_seq_coverage)
        circ_seq_coverage = str('%.2f'%float_circ_seq_coverage) + '%'
        #print(circ_seq_coverage,'+++++++++++++++++++')
        for circ_protaid_single in circ_protaid.strip(';').split(';'):
            #print(circ_protaid_single)
            protaid_single = circ_protaid_single.split('|')[0]
            protaid_num = circ_protaid_single.split('|')[1]
            protaid_score = circ_protaid_single.split('|')[2]
            protaid_start = corf_seq.find(protaid_single)
            protaid_end = protaid_start + len(protaid_single)
            protaid_num_score = int(protaid_num)*int(float(protaid_score))
            if protaid_single in corf_seq:
                list_protaid_single = list(range(protaid_start,protaid_end+1))
                dic_protaid_score.update({str(protaid_start) + '_' + str(protaid_end):protaid_num_score})
                list_protaid_multiple.append(list_protaid_single)
                #print(list_protaid_single)

        #print(dic_protaid_score)
        list_protaid_score_multiple = [0]*corf_seq_len
        for AA in range(0,corf_seq_len):
            for list_protaid_single in list_protaid_multiple:
                #print(list_protaid_single)
                protaid_num_score = dic_protaid_score[str(list_protaid_single[0]) + '_' + str(list_protaid_single[-1])]
                for AAA in list_protaid_single:
                    if AA == AAA:
                        #print('int(list_protaid_score_multiple[AA]',list_protaid_score_multiple[AA])
                        #print('list_protaid_single',list_protaid_single)
                        #print('dic_protaid_score',dic_protaid_score)
                        #print('int(dic_protaid_score[list_protaid_single]',int(dic_protaid_score[list_protaid_single]))
                        list_protaid_score_multiple[AA] = int(list_protaid_score_multiple[AA]) + int(dic_protaid_score[str(list_protaid_single[0]) + '_' + str(list_protaid_single[-1])])
        #print(list_protaid_score_multiple,'_____________________-')
        changecolor = colors.Normalize(vmin=1, vmax=100.0)

        #画一个方框
        y_list = []
        fig = plt.figure()
        if corf_seq_len <= 150:
            column_num = 20
        if corf_seq_len <= 400:
            column_num = 30
        if corf_seq_len > 400:
            column_num = 40
        corf_list_int = len(corf_seq)//column_num
        corf_list_remainder = len(corf_seq) - corf_list_int*column_num
        #print(corf_list_int)
        #print(corf_list_remainder)
        x = list(range(0,column_num))*corf_list_int + list(range(0,corf_list_remainder))
        #print(x)
        yyy_list = list(range(1,corf_list_int+1))
        yyy_list.sort(reverse=True)
        for yy in yyy_list:
            yy_list = [yy]*column_num
            #print(yy_list,'-----------')
            y_list.extend(yy_list)
        y = y_list + [0]*corf_list_remainder
        #print(y)
        list_circ_line_single_all = []
        for aa in range(0,len(corf_seq)):
            if list_circ_line_multiple != []:
                for list_circ_line_single in list_circ_line_multiple:
                    x_start = x[int(list_circ_line_single[0])] - 0.5
                    y_start = y[int(list_circ_line_single[0])]
                    plt.text(x_start,y_start,list_start_end[0],fontsize=300/column_num,color='purple')
                    x_end = x[int(list_circ_line_single[-1])] + 0.6
                    y_end = y[int(list_circ_line_single[-1])]
                    plt.text(x_end,y_end,list_start_end[1],fontsize=300/column_num,color='purple')
                    if aa in list_circ_line_single:
                        plt.text(x[aa],y[aa],corf_seq[aa],fontsize=300/column_num,color='purple')
                        list_circ_line_single_all.append(aa)
            if aa not in list_circ_line_single_all:
                plt.text(x[aa],y[aa],corf_seq[aa],fontsize=300/column_num,color='red')
            #if aa not in list_circ_line_single_all: 
            #plt.scatter(x[aa]+0.5, y[aa]+0.1, s = 90, marker = 's', c = list_protaid_score_multiple ,norm = changecolor,cmap= 'viridis')

            #plt.fill(x[aa],y[aa], fc='yellow')
            if str(aa) in circ_bsj_position:
                plt.text(x[aa]-0.5,y[aa],list_start_end[2],fontsize=600/column_num,color='black')
        plt.scatter([item + 0.5 for item in x], [item + 0.1 for item in y], s = 90, marker = 's', c = list_protaid_score_multiple,norm = changecolor,cmap= 'Blues')# ,norm = changecolor,cmap= 'Blues')
        plt.xlim(xmax=column_num+column_num*0.1,xmin=0)   #x轴
        #plt.tick_params(axis='x', width=0)   #隐藏x轴坐标
        #plt.xticks(alpha=0)                  #隐藏X轴数字
        y_max = corf_list_int*0.1
        if corf_list_int <= 20:
            y_max = 2
        plt.ylim(ymax=corf_list_int+y_max,ymin=-y_max*2)
        #plt.tick_params(axis='y', width=0)
        #plt.yticks(alpha=0)
        plt.text(1, -y_max/2, 'cORF:'+corf_id)
        plt.text(1, -y_max/1, 'Gene_name:'+circ_gene)
        plt.text(1, -y_max/0.66, 'Sequence cover:'+ circ_seq_coverage)
        #print(corf_id,'---------',circ_seq_coverage)
        plt.colorbar()                       #色度条
        plt.scatter(x,y,marker='o',alpha=0.0000006,cmap='viridis')
        plt.axis('off')
        #plt.savefig('try_ath_ms_draw/'+corf_id+'.pdf', bbox_inches='tight')
        plt.savefig(file_out+'Peptide_coverage_map/'+corf_id+'.pdf', bbox_inches='tight')
        plt.show()
        print(corf_id.strip().split(':')[0],':绘图完成')

def circ_size(circ_len):
    siz = 0.5
    if circ_len < 200:
        siz = 10
    elif circ_len < 400:
        siz = 4.4
    elif circ_len < 600:
        siz = 2.2
    elif circ_len < 800:
        siz = 1.5
    elif circ_len < 1000:
        siz = 1.2
    elif circ_len <= 1200:
        siz = 0.4
    return siz

def draw_circ_corf():
    dic_circ_corf_txt = {}
    file_circ_predict = open(file_out + 'circ_predict.txt','r')
    file_circ_predict_header = file_circ_predict.readline()
    #print(list_pro_long_score[:excellent_count])
    for row_circ_predict in file_circ_predict:
        circ_gene = row_circ_predict.strip().split('\t')[1]
        circ_corf_id = row_circ_predict.strip().split('\t')[3].split('-')[0]
        circ_corf_seq = row_circ_predict.strip().split('\t')[5].split('-')[0]
        circ_corf_seq_x = row_circ_predict.strip().split('\t')[7]
        pro_long_score = row_circ_predict.strip().split('\t')[14]

        if float(pro_long_score) not in list_pro_long_score[:excellent_count]:
            continue
        if circ_corf_id not in dic_circ_corf_txt.keys():
            dic_circ_corf_txt.update({circ_corf_id:circ_gene+'\t'+circ_corf_seq+'\t'+circ_corf_seq_x})
        elif circ_corf_id in dic_circ_corf_txt.keys():
            dic_value = dic_circ_corf_txt[circ_corf_id] + ';' + circ_gene+'\t'+circ_corf_seq+'\t'+circ_corf_seq_x
            dic_circ_corf_txt.update({circ_corf_id:dic_value})
   
    #print(str(dic_circ_corf_txt)[0:2000])
    dic_draw_circ_id_seq = {}
    file_circ_seq = open(circ_seq,'r')
    for row_circ_seq in file_circ_seq:
        if row_circ_seq[0] == '>':
            circ_id = row_circ_seq[1:].split(':')[0]
        elif row_circ_seq[0] != '>':
            circ_sequence = row_circ_seq.strip()
            if circ_id in dic_circ_corf_txt.keys():
                dic_draw_circ_id_seq.update({circ_id:circ_sequence})
                
    if not os.path.exists(file_out+'circRNA_Peptide_coverage_map'):
        os.system("mkdir %scircRNA_Peptide_coverage_map/"%(file_out))
    a,b = (0.0,0.0)
    #draw_sum = 10
    #draw_actual = 0
    for circ_corf_id in dic_circ_corf_txt.keys():
        #if draw_actual > draw_sum:
        #    break
        #draw_actual = draw_actual + 1
        r = 5
        circ_sequence = dic_draw_circ_id_seq[circ_corf_id]
        circ_seq_len = len(circ_sequence)
        circ_surplus = 0
        siz = circ_size(circ_seq_len)
        if circ_seq_len > 1200:
            circ_surplus = circ_seq_len - 1200
            circ_sequence = circ_sequence[0:599] + '--' + circ_sequence[-600:]
            circ_seq_len = len(circ_sequence)
        #
        fig = plt.figure()
        plt.text(0.0,5.0,'|',size=siz+1,color='red')
        circ_list_d = np.linspace(5/2*np.pi,1/2*np.pi,circ_seq_len+1)
        xc = a + r * np.cos(circ_list_d)
        yc = b + r * np.sin(circ_list_d)
        for i in range(0,circ_seq_len):
            plt.text(xc[i],yc[i],circ_sequence[i],size=siz/1.2,rotation=-i/360,color='black')
        if circ_surplus != 0 :
            plt.text(-0.5,-5.5,'surplus:'+str(circ_surplus)+'nucleotide')
        
        for row_circ_predict in dic_circ_corf_txt[circ_corf_id].split(';'):
            circ_gene = row_circ_predict.strip().split('\t')[0]
            circ_corf_seq = row_circ_predict.strip().split('\t')[1]
            circ_corf_seq_len = len(circ_corf_seq)
            circ_corf_seq_x = row_circ_predict.strip().split('\t')[2]
            aa = -1
            list_circ_corf_map = []
            for AA in circ_corf_seq_x:
                if AA != '|':
                    aa = aa + 1
                    if AA != 'X' and AA != '*':
                        list_circ_corf_map.append(aa)
            circ_corf_len = len(circ_corf_seq)
            circ_bsj_position_first = circ_corf_seq_x.find('|')
            #
            corf_list_d = np.linspace((0.25+((circ_bsj_position_first+1)*3)/(circ_seq_len+1))*2*np.pi,(0.25-(circ_corf_len-circ_bsj_position_first)*3/(circ_seq_len+1))*2*np.pi,circ_corf_len+1)
            r = r - 0.5
            if circ_corf_seq_len*3 >= circ_seq_len:
                r_poor = (circ_corf_seq_len*3)/circ_seq_len
                r_spiral = np.linspace(r,r-r_poor*0.3,circ_corf_seq_len+1)
                xco = []
                yco = []
                for i in range(0,circ_corf_seq_len+1):
                    xco_single = a + r_spiral[i]*np.cos(corf_list_d[i])
                    xco.append(xco_single)
                    yco_single = b + r_spiral[i]*np.sin(corf_list_d[i])
                    yco.append(yco_single)
            elif circ_corf_seq_len*3 < circ_seq_len:
                xco = a + r * np.cos(corf_list_d)
                yco = b + r * np.sin(corf_list_d)
            ang = np.linspace(-((circ_bsj_position_first*3+1)/(circ_seq_len+1))*360,(circ_corf_len-circ_bsj_position_first)*3/(circ_seq_len+1)*360,circ_corf_seq_len+1)
            
            for i in range(0,circ_corf_seq_len):
                if i not in list_circ_corf_map:
                    plt.text(xco[i],yco[i],circ_corf_seq[i],fontsize=siz+1.2,rotation=-ang[i],color='black')
                if i in list_circ_corf_map:
                    plt.scatter(xco[i],yco[i],marker='o',alpha=0.000000000006,color='blue',s=siz*10+15)
                    plt.text(xco[i],yco[i],circ_corf_seq[i],fontsize=siz+1.3,rotation=-ang[i],color='blue')
                    #plt.text(xco[i],yco[i],'__',fontsize=siz+3.5,rotation=-ang[i],color='blue')
        plt.text(-1.0,0.0,circ_corf_id,color='black',fontsize=siz+2.2)
        plt.text(-1.0,-0.5,'Gene:'+circ_gene,color='black',fontsize=siz+2.2)
        plt.axis("equal")
        plt.axis('off')
        plt.savefig(file_out + 'circRNA_Peptide_coverage_map/'+circ_corf_id+'.pdf',bbox_inches='tight')
        print(circ_corf_id+':绘图完成')


#circ_RNA名字和序列对应的字典
dic_circ_id_seq = circ_id_seq()
#整合tsv格式文件中的数据
if file_input_protein == 'none' and path_input_peptide != 'none':
    dic_id_peptide = tsv_circ_protaid()
#去除tsv文件中与线性蛋白序列一致的肽段
elif file_input_protein != 'none':
    if path_input_tsv != 'none' and path_input_peptide == 'none':
        tsv_circ_protaid()
    dic_id_peptide = tsv_remove_line_protaid()

#整合circRNA开放阅读框中多个bsj位点
dic_circ_orf_bsj_location = circ_orf_bsj_location()
#circRNA开放阅读框中与线性RNA翻译蛋白一致的蛋白序列
if file_input_protein != 'none':
    dic_circ_len_line_lacation = circ_crf_linepro_map()
elif file_input_protein == 'none':
    dic_circ_len_line_lacation = {}
#circRNA序列的正负链信息
dic_circ_strand = circ_forward_reverse()
#circRNA的可翻译预测结果
list_pro_long_score = circ_peptide_result()

#有质谱数据的分类文件
if file_input_PRIDE != 'none':
    dic_experiment,dic_raw_organization,list_control_experiment_raw_id,dic_PXD_control_raw,dic_PXD_experiment_raw = PRIDE_experiment()
    dic_mgf_circ_amount = mgf_amount_peptide()
    #print(dic_mgf_circ_amount)
    PRIDE_experiment_class()
    dic_circ_PXD_edge = edge_predict()
    circ_predict_PRIDE()

#肽覆盖率图
draw_msgf_second(file_out + 'circ_predict.txt',excellent_count)

#circRNA和开放阅读框
draw_circ_corf()


