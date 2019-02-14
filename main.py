#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
version:3.0
by:lss
"""

import numpy as np
import pandas as pd
import sys
reload(sys)
sys.setdefaultencoding("utf8")
import os
import re
import time
import collections
import smtplib
from email.mime.text import MIMEText
from email.mime.multipart import MIMEMultipart
from email.mime.base import MIMEBase
from email import encoders
import shutil
import subprocess

from lib.prepare import *

m_path = os.path.split(os.path.abspath(sys.argv[0]))[0]
print m_path
if len(sys.argv) > 1:
    cmd = sys.argv[1]

if not os.path.exists('config'):
    #shutil.copy(m_path+'/src/config', '.')
    print """
    Examples: 
    prepare config
    format: (delimiter \t)
    -sn excelfile 
    -strategy   IDT-PANEL 
    -filtmode   hard 
    -dir    absolute_path_of_rawdata 
    -md5    absolute_path_of_rawdata_of_md5_file 
    -R1 suffix_of_raw_R1 
    -R2 suffix_of_raw_R2
    Note:
    make sure excel colnames contain are 批次   到样日期    公司编码    原始编码    项目    性别    年龄    科室    浓度 (ng/ul)    体积 (μL)   总量(μg)    OD260/OD280 OD260/OD230 质检结果    上机日期    下机日期    Reads   Total_base(Mbases)  Q30(%)
    Reads, Total_base(Mbases), Q30(%) and they are not null; the rawdata and md5 fileare in the same single dir 
    """
    sys.exit()
    
log = open('log', 'a')

kwargs={'-sn':'', '-strategy':'', '-filtmode':'', '-dir':'', '-md5':'', '-R1':'', '-R2':''}

with open('config') as f:
    for l in f:
        l = l.strip('\n').split('\t')
        if len(l) != 0:
            for j in kwargs.keys():
    		    if l[0] == j:
    	    		kwargs.update({j:l[1]})
        
sn=kwargs.get('-sn')
strategy=kwargs.get('-strategy')
filtmode=kwargs.get('-filtmode')
dir = kwargs.get('-dir')
md5 = kwargs.get('-md5')
R1 = kwargs.get('-R1')
R2 = kwargs.get('-R2')
dirname = os.getcwd().split('/')[-2:]
if u'分析任务单' in sn:
    sep = u'分析任务单'
else:
    sep = u'产品信息表'

excel = pd.read_excel(sn)
excel = append_sample_txt(excel)
merge_columns = excel.columns.tolist()
d1 = append_gender(excel)
d1 = append_relation(d1)
d1 = append_pedigree(d1)
d1 = append_phenotype(d1)
d1 = append_father_mother(d1)
print d1
pedigree_dict = generate_pedigree(d1) # dictionary of pedigree samples, keys:names of patient samples, values:names of the whole family members including the patient sample
trio = pedigree_dict.keys() 
single = generate_single(d1)

#####Prepare ################################################################################
if 'prepare' in cmd:
    log.write('Start\t'+time.asctime( time.localtime(time.time()) )+'\n')
    generate_dbevn(strategy)
    print_panel()
    print '-filtmode: '+filtmode
    
    report_dirname = os.getcwd().split('/')[-1]
    for i in ['ped', 'peddy', 'gvcf', 'annotation', 'raw', 'CNV', report_dirname, 'depth']:
        if not os.path.exists(i):
            os.mkdir(i)
    if not os.path.exists(report_dirname+'/vcf'):
        os.makedirs(report_dirname+'/vcf')
    if d1.shape[1] > 1:
        if not os.path.exists(report_dirname+'/CNV'):
            os.makedirs(report_dirname+'/CNV')

    if u'公司编码' in d1.columns and u'原始编码' in d1.columns:
        d1_tmp = d1[[u'公司编码', u'原始编码']]
    d1_tmp.to_csv('rename', index=None, header=None, sep="\t")
    
    for i in d1['sample']:
    	if not os.path.exists(str(i)):
    		os.mkdir(str(i))
    d1['sample'].to_csv('list', header=None, index=None)
    d1['sample'].to_csv('CNV/id', header=None, index=None)
    d1['CNV'] = '/BP12_share/sslyu/bam/'+d1['sample']+'.sort.mkdup.bam'
    d1['CNV'].to_csv('CNV/list', header=None, index=None)
    for i in ['run2_sentieon.sh', 'cmd.sh', 'cp_rawdata.sh']:
        shutil.copy(m_path+'/src/'+i,os.getcwd())
    
    run2_sentieon = open(m_path+'/src/run2_sentieon.sh').readlines()
    for i in range(21,31):
        print run2_sentieon[i]
    
    if os.path.exists('annotation/annotation.sh'):
        os.remove('annotation/annotation.sh')
    if len(single) > 0:
        list_single = open('list_single', 'w')
        if strategy == "WES":
            cmd_single = "clean mapping precalling gvcf index_single gtgvcf vqsr left-normalize depth qcsum qcvcf annotation"
        else:
            cmd_single = "clean mapping precalling gvcf index_single gtgvcf hardfilt left-normalize depth qcsum qcvcf annotation"
        for i in single:
            sample = d1.loc[i]['sample'] # index:姓名
            print 'single '+sample
            ID = d1.loc[i][u'公司编码']
            list_single.write(sample+'\n')
            if os.path.exists(sample+'/cmd.sh'):
                os.remove(sample+'/cmd.sh')
            generate_cp(sample, ID, dir, md5, R1, R2, strategy)
            generate_cmd_part1(sample)
            generate_cmd_part2(sample, dirname)
            generate_single_ped(sample,d1)
            generate_run2_sentieon(sample,cmd_single)
        list_single.close()
    if len(trio) > 0:
        shutil.copy(m_path+'/src/trio_cmd.sh', os.getcwd())
        if not os.path.exists('trio'):
            os.mkdir('trio')
        cmd_trio = "clean mapping precalling gvcf index_trio depth qcsum"
        list_trio = open('list_trio', 'w')
        ofile = open('trio/list', 'w')
        for k in pedigree_dict:
            generate_ped(k,d1)
            sample = d1.loc[k]['sample']
            ofile.write(sample+'\n')
            generate_trio_cmd(k, filtmode, pedigree_dict, d1, dirname)
            for i in pedigree_dict[k]:
                sample = d1.loc[i]['sample']
                print 'trio '+sample
                ID = d1.loc[i][u'公司编码']
                list_trio.write(sample+'\n')
                generate_cp(sample, ID, dir, md5, R1, R2, strategy)
                generate_cmd_part1(sample)
                generate_run2_sentieon(sample,cmd_trio)
        list_trio.close()
        ofile.close()
    os.system('sort annotation/annotation.sh |uniq > annotation/tmp; mv annotation/tmp annotation/annotation.sh')        
    if os.path.exists('list_single'):
        os.system('for i in `cat list_single`; do cat ped/$i.mendel.ped >> peddy_tmp; done')
    if os.path.exists('trio/list'):
        os.system('for i in `cat trio/list`; do cat ped/$i.mendel.ped >> peddy_tmp; done')
    os.system('sort peddy_tmp |uniq |awk \'{OFS="\t"; print \'0\',$2,$3,$4,$5,$6}\' > ped/peddy.ped; rm peddy_tmp')
    with open('ped/peddy.ped') as f:
        for l in f:
            print l.strip('\n')
    log.write('Prepare\t'+time.asctime( time.localtime(time.time()) )+'\n')
###############################################################################################################

#######################################################################################################
if 'cp' in cmd:
    if u'公司编码' in d1.columns and u'原始编码' in d1.columns:
        d1_tmp = d1[[u'公司编码', u'原始编码']]
    d1_tmp.to_csv('rename', index=None, header=None, sep="\t")

    if len(single) > 0:
        for i in single:
            sample = d1.loc[i]['sample']
            ID = d1.loc[i][u'公司编码']
            generate_cp(sample, ID, dir, md5, R1, R2, strategy)
    if len(trio) > 0:
        for k in pedigree_dict:
            for i in pedigree_dict[k]:
                sample = d1.loc[i]['sample']
                ID = d1.loc[i][u'公司编码']
                generate_cp(sample, ID, dir, md5, R1, R2, strategy)
    log.write('cp\t'+time.asctime( time.localtime(time.time()) )+'\n')
#########################################################################################################

#########################################################################################################
if 'ped' in cmd:
    if len(single) > 0:
        list_single = open('list_single', 'w')
        for i in single:
            sample = d1.loc[i]['sample']
            generate_single_ped(sample,d1)
            list_single.write(sample+'\n')
        list_single.close()
    if len(trio) > 0:
        list_trio = open('list_trio', 'w')
        for k in pedigree_dict:
            generate_ped(k,d1)
            for i in pedigree_dict[k]:
                sample = d1.loc[i]['sample']
                list_trio.write(sample+'\n')
        list_trio.close()
    if os.path.exists('list_single'):
        os.system('for i in `cat list_single`; do cat ped/$i.mendel.ped >> peddy_tmp; done')
    if os.path.exists('trio/list'):
        os.system('for i in `cat trio/list`; do cat ped/$i.mendel.ped >> peddy_tmp; done')
    os.system('sort peddy_tmp |uniq |awk \'{OFS="\t"; print \'0\',$2,$3,$4,$5,$6}\' > ped/peddy.ped; rm peddy_tmp')
    with open('ped/peddy.ped') as f:
        for l in f:
            print l.strip('\n')
    log.write('ped\t'+time.asctime( time.localtime(time.time()) )+'\n')
####################################################################################################

#####################################################################################################
if 'run' in cmd:
#######Run
##### 1 copy rawdata ###########
    jobID_d = {}
    for sample in d1['sample']:
        sample = str(sample)
        jobID = submit(sample, 'cp_rawdata.sh')
        jobID_d[sample] = jobID
    print jobID_d
    
    ##### 2 mapping and gvcf ######
    check_copy = copy_check(jobID_d)
    while check_copy:
        log.write('Copy\t'+time.asctime( time.localtime(time.time()) )+'\n')
        break
    while check_copy:
        single_jobID = {}
        trio_jobID = {}
        for i in single:
            sample = d1.loc[i]['sample'] # index:姓名
            sample = str(sample)
            jobID = submit(sample, 'cmd.sh')
            single_jobID[sample] = jobID
        trio_jobID = cmd_trio_submit(d1)
        break
    print single_jobID
    print trio_jobID
    if os.path.exists('list_single') and os.path.exists('list_trio'):
        jobID_trio = cmd_trio_check(trio_jobID, d1)
        log.write('Sentieon_trio\t'+time.asctime( time.localtime(time.time()) )+'\n')
        check_cmd_single = cmd_single_check(single_jobID)
        if len(jobID_trio) > 0:
            check_cmd_trio = True
        while check_cmd_trio and check_cmd_single:
            log.write('Sentieon_single\t'+time.asctime( time.localtime(time.time()) )+'\n')
            break
    elif os.path.exists('list_single'):
        check_cmd_single = cmd_single_check(single_jobID)
        while check_cmd_single:
            log.write('Sentieon\t'+time.asctime( time.localtime(time.time()) )+'\n')
            break
    else:
        jobID_trio = cmd_trio_check(trio_jobID, d1)
        if len(jobID_trio) > 0:
            check_cmd_trio = True
        while check_cmd_trio:
            log.write('Sentieon\t'+time.asctime( time.localtime(time.time()) )+'\n')
            break
    
    if os.path.exists('list_trio'):
        check_trio = trio_check(jobID_trio)
        while check_trio:
            log.write('Trio\t'+time.asctime( time.localtime(time.time()) )+'\n')
            break

##########################################################################################################

#### 4 pedigree and sex check ############
if 'check' in cmd:
    m_path2 = '\/'.join(m_path.split('/'))
    print m_path2
    cmd = 'sed -i "14s/basedir.*$/basedir='+m_path2+'\/src/" '+m_path+'/src/trio_cmd.sh'
    print cmd
    os.system('cd peddy; sed -i "14s/basedir.*$/basedir='+m_path2+'\/src/" '+m_path+'/src/trio_cmd.sh; sbatch '+m_path+'/src/trio_cmd.sh; cd ..')
    while True:
        time.sleep(60)
        if os.path.exists('peddy/pedtrio/results/'):
            break
    peddy_check()
    log.write('Check\t'+time.asctime( time.localtime(time.time()) )+'\n')

##### 5 jiaofuxinxi and database ############
if 'excel' in cmd:
    qcsum = generate_qcsum(d1)
    qcvcf = generate_qcvcf(d1)
    jiaofuxinxi_full = generate_jiaofuxinxi(d1,merge_columns, qcsum,strategy,sep,sn)
    navicat = generate_navicat(d1, jiaofuxinxi_full, sn, sep, strategy, qcsum, qcvcf)
    if not os.path.exists('unmapped_duplicate'):
        subprocess.call(["sh",m_path+"/src/unmapped_duplicate.sh","list"])
    table_map_all = pd.read_table('unmapped_duplicate', header=None, index_col=0)    
    table_map_all.columns = ['Unmapped reads', 'UNPAIRED_READ_DUPLICATES', 'READ_PAIR_DUPLICATES', 'Percent duplication']
    generate_jiesuan(navicat,table_map_all,strategy)
    log.write('Excel\t'+time.asctime( time.localtime(time.time()) )+'\n')

##### 6 qcsum check ############
if 'qcsum' in cmd:
    qcsum_check(strategy)
    log.write('qcsum\t'+time.asctime( time.localtime(time.time()) )+'\n')

###### 7 calculate CNV ###########
if 'cnv' in cmd:
    if not os.path.exists('CNV/id') and not os.path.exists('CNV/list'):
        d1['sample'].to_csv('CNV/id', header=None, index=None)
        d1['CNV'] = '/BP12_share/sslyu/bam/'+d1['sample']+'.sort.mkdup.bam'
        d1['CNV'].to_csv('CNV/list', header=None, index=None)

    ofile = open('step5_backup.sh', 'w')
    backup(ofile, dirname, sn, sep)
    ofile.close()
    os.system('sh step5_backup.sh')
    log.write('CNV\t'+time.asctime( time.localtime(time.time()) )+'\n')

###### 8 Send Email ###############
if 'email' in cmd:
    n = d1.shape[0]
    create_email(n)
    if 'email1' in cmd:
        send_email(1, n)
    if 'email2' in cmd:
        send_email(2, n)
    log.write('Email\t'+time.asctime( time.localtime(time.time()) )+'\n')
######################################################################################

###### 9 calculating exon coverage #############################################
if 'exon' in cmd:
    for i in ['tsv', 'pdf']:
        if not os.path.exists(i):
            os.mkdir(i)
    for i in d1['sample']:
        i = str(i)
        exon_plot(i)
    log.write('Exon\t'+time.asctime( time.localtime(time.time()) )+'\n')
##########################################################################################

    log.write('End\t'+time.asctime( time.localtime(time.time()) )+'\n')

log.close()
