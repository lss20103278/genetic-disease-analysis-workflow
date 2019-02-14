#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
version:2.6
by:lss
"""

import numpy as np
import pandas as pd
import sys
reload(sys)
sys.setdefaultencoding("utf8")
import os
import re
import collections
import subprocess
import time

m_path = os.path.split(os.path.realpath(sys.argv[0]))[0]

def append_sample_excel(df):
    ID = []
    if u'公司编码' in df.columns and u'原始编码' in df.columns:
        for i in range(len(df[u'原始编码'])):
            if df[u'原始编码'][i] == 'nan':
                ID.append(str(df.iloc[i][u'公司编码']))
            else:
                ID.append(str(df.iloc[i][u'原始编码']))
    df['sample'] = ID
    return df

def append_sample_txt(df):
    ID = []
    if u'公司编码' in df.columns and u'原始编码' in df.columns:
        for i in range(len(df[u'原始编码'].isnull())):
            if df[u'原始编码'].isnull()[i]:                
                ID.append(df.iloc[i][u'公司编码'])
            else:
                ID.append(df.iloc[i][u'原始编码'])
    df['sample'] = ID
    return df

def append_gender(df):
    df['gender'] = df[u'性别'].apply(lambda x:'1' if u'男' in x else '2' if u'女' in x else 'unknown')
    return df

def generate_pedigree(df):
    pedigree = {}
    df.index = df[u'项目']
    for i in df[u'项目']:
        pedigree[i] = []
    for k in pedigree:
        for i in df[u'项目']:
            if k in i: # the name of the child should be totally included in other family member's names
                pedigree[k].append(i)
    list_keys = pedigree.keys()
    for k in list_keys:
        if len(pedigree[k]) == 1:
            pedigree.pop(k)
    df.index = df[u'项目']
    return pedigree

def append_relation(df):
    df.index = df[u'项目']
    if 'relationship' not in df.columns:
        df['relationship'] = None
    pedigree_dict = generate_pedigree(df)
    trio = pedigree_dict.keys()
    for i in df.index:
        for j in trio:
            if j in i:
                if j == i:
                    df.loc[i,'relationship'] = u'子'
                elif u'父' in i:
                    df.loc[i,'relationship'] = u'父'  #there is no u'父' in i
                elif u'母' in  i:
                    df.loc[i,'relationship'] = u'母'  #there is no u'母' in i
                else:
                    df.loc[i,'relationship'] = 'other'
    df[u'家系关系'] = None
    df[u'家系关系'] = df['relationship']
    return df

def append_pedigree(df):
    pedigree_dict = generate_pedigree(df)
    df.index = df[u'项目']
    df['familyname'] = df['sample']
    for i in df.index:
        for k in pedigree_dict:
            if k in i:  # the name of the child should be totally included in other family member's names
                df.loc[i,'familyname'] = df.loc[k]['sample']
    df[u'样本间关系'] = None
    for i in df.index:
        for k in pedigree_dict:
            if k in i:
                df.loc[i,u'样本间关系'] = 'fam'+df.loc[k]['sample']
    return df                    

def flatten_list(nested):
    if isinstance(nested, list):
        for sublist in nested:
            for item in flatten_list(sublist):
                yield item
    else:
        yield nested

def append_phenotype(df):
    pedigree_dict = generate_pedigree(df)
    trio = pedigree_dict.keys()
    try:
        single = [i for i in df[u'项目'] if i not in flatten_list(pedigree_dict.values())] # dict.values()
    except:
        pass
    df['phenotype1'] = '0'
    df['phenotype2'] = None
    if len(trio) > 0 and len(single) > 0:
        for i in df.index: 
            if i in trio:
                df.loc[i,'phenotype2'] = '2'
            elif i in single:
                df.loc[i,'phenotype2'] = '2'
            else:
                df.loc[i,'phenotype2'] = '1'
    elif len(trio) > 0:
        for i in df.index:
            if i in trio:
                df.loc[i,'phenotype2'] = '2'
            else:
                df.loc[i,'phenotype2'] = '1'
    else:
        for i in df.index:
            df.loc[i,'phenotype2'] = '2'
    return df

def generate_single(df):
    pedigree_dict = generate_pedigree(df)
    trio = pedigree_dict.keys()
    trio_m = []
    for i in df[u'项目']:
        for k in pedigree_dict:
            trio_m.extend(pedigree_dict[k])
    single = [i for i in df[u'项目'] if i not in trio_m]
    return single

def append_father_mother(df):
    df['father'] = '0'
    df['mother'] = '0'
    pedigree_dict = generate_pedigree(df)
    trio = pedigree_dict.keys()
    try:
        single = [i for i in df[u'项目'] if i not in pedigree_dict.values()[0]]
    except:
        single = df[u'项目'].tolist()
    df.index = df[u'项目']
    for i in trio:
        relation = ''
        for j in pedigree_dict[i]:
            relation = relation+df.loc[j]['relationship']
        for j in df.index:
            if j in pedigree_dict[i]:
                if u'父' in relation and u'母' in relation:
                    if j == i:
                        for k in pedigree_dict[i]:
                            if u'父' in k:
                                father_name = k
                            if u'母' in k:
                                mother_name = k
                        df.loc[j,'father'] = df.loc[father_name]['sample']
                        df.loc[j,'mother'] = df.loc[mother_name]['sample']
                    else:
                        if len(pedigree_dict[i]) == 3:
                            df.loc[j,'father'] = '0'
                            df.loc[j,'mother'] = '0'
                        else:
                            if u'姐' in j or u'妹' in j or u'哥' in j or u'弟' in j:
                                for k in pedigree_dict[i]:
                                    if u'父' in k:
                                        father_name = k
                                    if u'母' in k:
                                        mother_name = k
                                df.loc[j,'father'] = df.loc[father_name]['sample']
                                df.loc[j,'mother'] = df.loc[mother_name]['sample']
                            elif u'爷' not in relation and u'奶' not in relation and u'外婆' not in relation and u'外公' not in relation:
                                if u'父' in j:
                                    df.loc[j,'father'] = '0'
                                    df.loc[j,'mother'] = '0'
                                if u'母' in j:
                                    df.loc[j,'father'] = '0'
                                    df.loc[j,'mother'] = '0'
                            else:
                                father = raw_input('father ID of '+j+'(if not exist, please enter 0): ')
                                df.loc[j,'father'] = father
                                mother = raw_input('mother ID of '+j+'(if not exist, please enter 0):')
                                df.loc[j,'mother'] = mother
                else:
                    if j == i:
                        if u'父' in relation:
                            for k in pedigree_dict[i]:
                                if u'父' in k:
                                    father_name = k
                            df.loc[j,'father'] = df.loc[father_name]['sample']
                        if u'母' in relation:
                            for k in pedigree_dict[i]:
                                if u'母' in k:
                                    mother_name = k
                            df.loc[j,'mother'] = df.loc[mother_name]['sample']
                    elif u'父' in j:
                        if u'爷' not in relation and u'奶' not in relation and u'外婆' not in relation and u'外公' not in relation:
                            df.loc[j,'father'] = '0'
                            df.loc[j,'mother'] = '0'
                    elif u'母' in j:
                        if u'爷' not in relation and u'奶' not in relation and u'外婆' not in relation and u'外公' not in relation:
                            df.loc[j,'father'] = '0'
                            df.loc[j,'mother'] = '0'
                    else:
                        father = raw_input('father ID of '+j+'(if not exist, please enter 0):')
                        df.loc[j,'father'] = father
                        mother = raw_input('mother ID of '+j+'(if not exist, please enter 0):')
                        df.loc[j,'mother'] = mother
    return df                        

def generate_single_ped(k,d):
    d.index = d['sample']
    ped_mendel = d.loc[k][['sample', 'father', 'mother', 'gender', 'phenotype2']]
    fname = d.loc[k]['sample']+'/'+d.loc[k]['sample']+'.mendel.ped'
    fname = 'ped/'+d.loc[k]['sample']+'.mendel.ped'
    ped = open(fname, 'w')
    ped.write(str(d.loc[k]['sample'])+'\t'+'\t'.join(ped_mendel)+'\n')
    ped.close()
    d.index = d[u'项目']

def generate_ped(k,d):
    pedigree_dict = generate_pedigree(d)
    ped_f = d.loc[pedigree_dict[k]][['sample', 'father', 'mother', 'gender', 'phenotype1']]
    ped_f.index = ped_f['sample']
    ped_mendel = d.loc[pedigree_dict[k]][['sample', 'father', 'mother', 'gender', 'phenotype2']]
    ped_mendel.index = ped_mendel['sample']
    f1name = 'trio/ped/'+d.loc[k]['sample']+'.ped'
    f2name = 'trio/ped/'+d.loc[k]['sample']+'.mendel.ped'
    f1name = 'ped/'+d.loc[k]['sample']+'.ped'
    f2name = 'ped/'+d.loc[k]['sample']+'.mendel.ped'
    ped1 = open(f1name, 'w')
    ped2 = open(f2name, 'w')
    relation = []
    for i in pedigree_dict[k]:
        relation.append(d.loc[i]['relationship'])
    a = ' '.join(relation)
    key1 = d.loc[k]['sample']
    ped1.write(key1+'\t'+'\t'.join(ped_f.loc[key1])+'\n')
    ped2.write(key1+'\t'+'\t'.join(ped_mendel.loc[key1])+'\n')
    if u'父' in a and u'母' in a:
        for j in [u'父', u'母']:
            for i in pedigree_dict[k]:
                if i != k and j in i:
                    key2 = d.loc[i]['sample']
                    ped1.write(key1+'\t'+'\t'.join(ped_f.loc[key2])+'\n')
                    ped2.write(key1+'\t'+'\t'.join(ped_mendel.loc[key2])+'\n')
    	n_c_f_m = [i for i in pedigree_dict[k] if u'父' not in i and u'母' not in i and i != k]
        for i in n_c_f_m:
            key2 = d.loc[i]['sample']
            ped1.write(key1+'\t'+'\t'.join(ped_f.loc[key2])+'\n')
            ped2.write(key1+'\t'+'\t'.join(ped_mendel.loc[key2])+'\n')
    elif u'父' in a or u'母' in a:
        parent = ''
        if u'父' in a:
            parent = u'父'
        else:
            parent = u'母'
        for i in pedigree_dict[k]:
            if parent in i and i != k:
                key2 = d.loc[i]['sample']
                ped1.write(key1+'\t'+'\t'.join(ped_f.loc[key2])+'\n')
                ped2.write(key1+'\t'+'\t'.join(ped_mendel.loc[key2])+'\n')
        n_c_f_m = [i for i in pedigree_dict[k] if parent not in i and i != k]
        for i in n_c_f_m:
            key2 = d.loc[i]['sample']
            ped1.write(key1+'\t'+'\t'.join(ped_f.loc[key2])+'\n')
            ped2.write(key1+'\t'+'\t'.join(ped_mendel.loc[key2])+'\n')
    else:
    	n_c = [i for i in pedigree_dict[k] if i != k]
    	for i in n_c:
    	    key2 = d.loc[i]['sample']
    	    ped1.write(key1+'\t'+'\t'.join(ped_f.loc[key2])+'\n')
    	    ped2.write(key1+'\t'+'\t'.join(ped_mendel.loc[key2])+'\n')
    ped1.close()
    ped2.close()

def generate_dbevn(strategy):
    dbevn = open('dbevn.sh', 'w')
    if strategy != 'none':
        with open(m_path+'/src/dbevn.sh') as f:
            for l in f:
                if l.startswith('#panel'):
                    if strategy in l:
                        dbevn.write(l[1:])
                    else:
                        dbevn.write(l)
                else:
                	dbevn.write(l)
    else:
        with open(m_path+'/src/dbevn.sh') as f:
            for l in f:
                dbevn.write(l)
        panel = raw_input('please enter the absolute path of the bed:')
        dbevn.write('panel=\"'+panel+'\"\n')
    dbevn.close()

def generate_run2_sentieon(path,str):
    global m_path
    run2_sentieon = open(path+'/run2_sentieon.sh', 'w')
    n = 0
    with open(m_path+'/src/run2_sentieon.sh') as f:
        for l in f:
            n = n+1
            if n>0 and n<= 16:
                run2_sentieon.write(l)
            elif n>17 and n<24:
                run2_sentieon.write(l)
            elif n>= 25:
                run2_sentieon.write(l)
            elif n==24:
                run2_sentieon.write('cmd="'+str+'"\n')
            elif n==17:
                run2_sentieon.write('basedir='+m_path+'/src\n')
            else:
                continue
    run2_sentieon.close()

def generate_cp(sample,ID,dir,md5,R1,R2,strategy):
    cp = open(sample+'/cp_rawdata.sh', 'w')
    cp.write('#!/bin/sh\n\n#SBATCH -J '+sample+'\n#SBATCH -p BP10\n#SBATCH -N 1\n#SBATCH -n 4\n\n')
    cp.write('sample='+sample+'\nID='+ID+'\n\ndir='+dir+'\nmd5='+md5+'\nR1='+R1+'\nR2='+R2+'\nstrategy='+strategy+'\n\n')
    #cp.write('mkdir -p /anaData/anaData004/children_hos_genetic/rawdata/2018/$strategy/$ID\n')
    #cp.write('awk \'{if(match($2,"\'$ID\'")){print $0}}\' $md5 >> /anaData/anaData004/children_hos_genetic/rawdata/2018/$strategy/$ID/$ID.md5\n')
    #cp.write('cp $dir/*$ID*$R1 /anaData/anaData004/children_hos_genetic/rawdata/2018/$strategy/$ID\n')
    #cp.write('cp $dir/*$ID*$R2 /anaData/anaData004/children_hos_genetic/rawdata/2018/$strategy/$ID\n')
    cp.write('mkdir -p /anaData/anaData004/children_hos_genetic/rawdata/2019/$strategy/$ID\n')
    cp.write('awk \'{if(match($2,"\'$ID\'")){print $0}}\' $md5 >> /anaData/anaData004/children_hos_genetic/rawdata/2019/$strategy/$ID/$ID.md5\n')
    cp.write('cp $dir/*$ID*$R1 /anaData/anaData004/children_hos_genetic/rawdata/2019/$strategy/$ID\n')
    cp.write('cp $dir/*$ID*$R2 /anaData/anaData004/children_hos_genetic/rawdata/2019/$strategy/$ID\n')
    cp.write('cp $dir/*$ID*$R1 ../raw/$sample\_R1.fastq.gz\n')
    cp.write('cp $dir/*$ID*$R2 ../raw/$sample\_R2.fastq.gz\n')

    cp.write('for j in $R1 $R2; do md5sum $dir/*$ID*$j >> ../original_md5; done\n')
    cp.write('for j in R1 R2; do md5sum ../raw/$sample\_$j*.fastq.gz >> ../cp_md5; done\n')
    #cp.write('for j in $R1 $R2; do md5sum /anaData/anaData004/children_hos_genetic/rawdata/2018/$strategy/$ID/*$ID*$j >> ../backup_md5; done\n')
    cp.write('for j in $R1 $R2; do md5sum /anaData/anaData004/children_hos_genetic/rawdata/2019/$strategy/$ID/*$ID*$j >> ../backup_md5; done\n')
    cp.write('for j in $R1 $R2; do a=`grep $ID ../original_md5 |grep $j |cut -d" " -f 1`; b=`grep $ID $md5 |grep $j |cut -d" " -f 1`; if [ "$a" != "$b" ]; then echo -e "original "$ID"_"$j" is problematic" >> ../check_md5.log; echo -e "original "$ID"_"$j" is problematic"; fi; done\n')
    cp.write('pattern_original=($R1 $R2); pattern_copied=(R1 R2); for j in 0 1; do original=`grep $ID ../original_md5 |grep ${pattern_original[$j]} |cut -d" " -f 1`; copied=`grep $sample ../cp_md5 |grep ${pattern_copied[$j]} |cut -d" " -f 1`; if [ "$original" != "$copied" ]; then echo -e $sample"_"$ID"_"${pattern_original[$j]}" is not copied correctly" >> ../check_md5.log; echo -e $sample"_"$ID"_"${pattern_original[$j]}" is not copied correctly"; fi; done\n')
    cp.write('for j in $R1 $R2; do original=`grep $ID ../original_md5 |grep $j |cut -d" " -f 1`; backup=`grep $ID ../backup_md5 |grep $j |cut -d" " -f 1`; if [ "$original" != "$backup" ]; then echo -e $ID"_"$j" is not backup correctly" >> ../check_md5.log; echo -e $ID"_"$j" is not backup correctly"; fi; done\n')
    cp.close()
    
def generate_cmd_part1(sample):
    single_cmd = open(sample+'/cmd.sh', 'w')
    single_cmd.write('#!/bin/sh\n\n#SBATCH -J '+sample+'\n#SBATCH -p BP10\n#SBATCH -N 1\n#SBATCH -n 4\n\n')
    single_cmd.write('sample='+sample+'\n\n')
    single_cmd.write('sh run2_sentieon.sh $sample $SLURM_NPROCS\n')
    single_cmd.write('samtools depth -b /DATA/sslyu/Project/Genetics_children_hospital/all_ch_freq_gene.bed /BP12_share/sslyu/bam/'+sample+'.sort.mkdup.bam > ../depth/'+sample+'.depth\n')
    single_cmd.write('scp /BP12_share/sslyu/bam/'+sample+'.sort.mkdup.bam* klyang@122.112.248.194:/media/bpdata/data/exon_temp2/\n')
    single_cmd.close()

def generate_cmd_part2(sample, dirname):
    single_cmd = open(sample+'/cmd.sh', 'a')
    single_cmd.write('mv 2_mapping/$sample.depth.sample_gene_summary ../annotation\n')
    single_cmd.write('mv 3_variants/$sample.{ann.hg19_multianno.txt,CADD,maf,link} ../annotation\n')
    single_cmd.write('cd ../annotation\n')
    single_cmd.write('echo $sample >> list\n')
    single_cmd.write('python '+m_path+'/src/score_re.py -sn $sample\n')
    single_cmd.write('python '+m_path+'/src/annotation_filt.py -sn $sample -c_f_m c\n')
    single_cmd.write('mv '+sample+'_children_hospital*xlsx ../'+dirname[-1]+'\n')
    single_cmd.write('mv ../'+sample+'/3_variants/'+sample+'.vcf ../'+dirname[-1]+'/vcf\n')
    single_cmd.write('echo "end `date`" >> ../$sample/finished\n')
    single_cmd.close()
    ofile = open('annotation/annotation.sh', 'a')
    ofile.write('python '+m_path+'/src/annotation_filt.py -sn '+sample+' -c_f_m c\n')
    ofile.close()

def generate_trio_cmd(k, filtmode, pedigree_dict, d1, dirname):              
    sample = d1.loc[k]['sample']
    path = 'trio/'+sample
    if not os.path.exists(path):
        os.makedirs(path)
    trio_cmd = open(path+'/trio_cmd.sh', 'w')
    trio_cmd.write('#!/bin/sh\n\n#SBATCH -J '+sample+'\n#SBATCH -p BP10\n#SBATCH -N 1\n#SBATCH -n 4\n\nsample='+sample+'\n\n')
    trio_cmd.write('filtmode=\''+filtmode+'\'\npath=`pwd`\nvar=`echo $path | awk -F \'/\' \'{print $NF}\'`\n')
    trio_cmd.write('if [ $var = \'peddy\' ]\nthen\nsh '+m_path+'/src/trio.sh $sample $SLURM_NPROCS\n')
    trio_cmd.write('else\nsh '+m_path+'/src/trio.sh $sample $SLURM_NPROCS $filtmode\n')
    trio_cmd.write('[ -e ../../$sample/2_mapping/$sample.depth.sample_gene_summary ] && mv ../../$sample/2_mapping/$sample.depth.sample_gene_summary ../../annotation\n')
    trio_cmd.write('mv segtrio/$sample.{ann.hg19_multianno.txt,CADD,maf,link} ../../annotation\n')
    trio_cmd.write('cd ../../annotation\n')
    trio_cmd.write('echo $sample >> list\n')
    trio_cmd.write('python '+m_path+'/src/score_re.py -sn $sample\n')
    relation = []
    for j in pedigree_dict[k]:
        relation.append(d1.loc[j]['relationship'])
    a = ' '.join(relation)
    print a
    annotationfile = open('annotation/annotation.sh', 'a')
    if u'父' in a and u'母' in a:
        trio_cmd.write('python '+m_path+'/src/annotation_filt.py -sn $sample -c_f_m c_f_m\n')
        annotationfile.write('python '+m_path+'/src/annotation_filt.py -sn '+sample+' -c_f_m c_f_m\n')
    else:
        if u'父' in a:
            trio_cmd.write('python '+m_path+'/src/annotation_filt.py -sn $sample -c_f_m c_f\n')
            annotationfile.write('python '+m_path+'/src/annotation_filt.py -sn '+sample+' -c_f_m c_f\n')
        else:
            trio_cmd.write('python '+m_path+'/src/annotation_filt.py -sn $sample -c_f_m c_m\n')
            annotationfile.write('python '+m_path+'/src/annotation_filt.py -sn '+sample+' -c_f_m c_m\n')
    trio_cmd.write('mv '+sample+'_children_hospital*xlsx ../'+dirname[-1]+'\n')
    trio_cmd.write('mv ../trio/'+sample+'/sep/*.vcf ../'+dirname[-1]+'/vcf\n')
    trio_cmd.write('echo "end `date`" >> ../$sample/finished\n')
    trio_cmd.write('fi\n')
    trio_cmd.close()
    annotationfile.close()

def print_panel():
    with open('dbevn.sh') as f:
        for l in f:
            if l.startswith('panel'):
                print l.strip('\n')

def append_fenxiliucheng(d1):
    trio = [item for item, count in collections.Counter(d1['familyname'].tolist()).items() if count > 1]
    single = [item for item, count in collections.Counter(d1['familyname'].tolist()).items() if count == 1]
    ped_v = []
    for j in d1['sample']:
        if j in single:
            ped_v.append(0)
        else:
            ped_v.append(1)
    d1['pedigree'] = ped_v
    d1[u'分析流程'] = d1['pedigree'].map(lambda x:'single' if x==0 or x=='nan' else 'trio')
    return d1

def generate_qcsum(d1):
    qcsum = pd.DataFrame(columns=['sample', '1Xcov', '20Xcov', 'AvgDepth', 'duplication'])
    for j in d1['sample']:
        fn = j+'/2_mapping/'+j+'.qcsum'
        df = pd.read_table(fn, dtype=str)
        qcsum = pd.concat([qcsum,df], ignore_index=True)
    qcsum.index = qcsum['sample']
    qcsum.to_csv('all.qcsum', sep="\t", index=None)
    return qcsum

def generate_qcvcf(d1):
    trio = [item for item, count in collections.Counter(d1['familyname'].tolist()).items() if count > 1]
    single = [item for item, count in collections.Counter(d1['familyname'].tolist()).items() if count == 1]
    qcvcf = pd.DataFrame(columns=['sample', 'indel_genomic', u'snp_genomic', u'ts/tv_genomic',u'snp_exonic', u'ts/tv_exonic'])
    qcvcf_list = []
    d1.index = d1['sample']
    for j in d1['sample'].tolist():
        if j in single:
            qcvcf_list.append(j)
        elif j in trio:
            qcvcf_list.append(j)
        else:
            continue
    for j in qcvcf_list:
        if j in single:
            fn = j+'/3_variants/'+j+'.qcvcf'
        else:
            fn = 'trio/'+j+'/triotmp/'+j+'.qcvcf'
        df = pd.read_table(fn, dtype=str)
        qcvcf = pd.concat([qcvcf,df], ignore_index=True)
    qcvcf.index = qcvcf['sample']
    qcvcf.to_csv('all.qcvcf', sep="\t", index=None)
    return qcvcf

def generate_jiaofuxinxi(d1,merge_columns, qcsum,strategy,sep,sn):
    d1 = append_fenxiliucheng(d1)
    d1_merge = d1[merge_columns]
    qcsum_merge = qcsum.drop(u'duplication', axis=1)
    jiaofuxinxi_table = pd.merge(d1_merge, qcsum_merge, on="sample") # notice the type of the value
    jiaofuxinxi_full = pd.merge(d1, qcsum_merge, on="sample")
    l2 = [u'策略', 'vcf', u'质检结果']
    l3 = [strategy, 'available', u'合格']
    for j in range(len(l2)):
        if l2[j] not in jiaofuxinxi_table.columns:
            jiaofuxinxi_table.insert(len(jiaofuxinxi_table.columns.tolist()), l2[j], l3[j])
            jiaofuxinxi_full.insert(len(jiaofuxinxi_full.columns.tolist()), l2[j], l3[j])
    jiaofuxinxi_table = jiaofuxinxi_table.drop('sample', axis=1)
    ofile = sn.split(sep)[0]+'交付信息表'+u'）'.join(sn.split(sep)[1].split(')'))
    wb=pd.ExcelWriter(ofile,engine='openpyxl')
    jiaofuxinxi_table.to_excel(wb,index=False)
    wb.save()
    return jiaofuxinxi_full

def generate_navicat(d1, jiaofuxinxi_full, sn, sep, strategy, qcsum, qcvcf):
    header = pd.read_excel(''+m_path+'/doc/navicat_header.xlsx')
    #total_table = pd.read_excel(''+m_path+'/bk/navicat_bak.xlsx')
    total_table = pd.read_excel('/DATA/sslyu/Project/Genetics_children_hospital/navicat_bak.xlsx')
    colname3 = header.columns
    sample_n = d1.shape[0]
    start = int(total_table['order'].tolist()[-1])+1
    header[u'order'] = range(int(start), int(start)+sample_n)
    jiaofuxinxi_full.rename(columns={u'原始编码':u'样本ID', u'科室':u'申请科室'}, inplace=True)
    #jiaofuxinxi_full_columns_ix = [0,2,3,-2,4,5,6,-12,7,16,17,18,8,9,10,11,12,13,13,-6,-1]
    jiaofuxinxi_full_columns_ix = [0,2,3,-2,4,5,6,-12,7,17,18,19,8,9,10,11,12,13,13,-6,-1]
    jiaofuxinxi_full_columns = []
    for i in jiaofuxinxi_full_columns_ix:
        jiaofuxinxi_full_columns.append(jiaofuxinxi_full.columns.tolist()[i])
    print '\t'.join(jiaofuxinxi_full_columns)
    header[[u'批次', u'样本编号', u'样本ID', u'策略', u'姓名', u'性别', u'年龄',u'样本间关系', u'申请科室',u'exp_reads',u'exp_bases(Mb)', u'exp_Q30', u'exp_浓度(ng/ul)', u'exp_体积(ul)',u'exp_质量(ug)', u'exp_OD260/OD280', u'exp_OD260/OD230', u'exp_QC结论',u'ana_qc_comment',u'ana_pepline', u'ana_vcf']] = jiaofuxinxi_full[jiaofuxinxi_full_columns[:21]]
    #header[[u'批次', u'样本编号', u'样本ID', u'策略', u'姓名', u'性别', u'年龄',u'样本间关系', u'申请科室',u'exp_reads',u'exp_bases(Mb)', u'exp_Q30', u'exp_浓度(ng/ul)', u'exp_体积(ul)',u'exp_质量(ug)', u'exp_OD260/OD280', u'exp_OD260/OD230', u'exp_QC结>论',u'ana_qc_comment',u'ana_pepline', u'ana_vcf']] = jiaofuxinxi_full[jiaofuxinxi_full_columns] ValueError: Wrong number of items passed 2, placement implies 1
    header.index = d1.index # very important
    header[u'家系关系'] = None
    header[u'家系关系'] = d1[u'家系关系']
    header[[u'process_送样日期', u'process_上机日期', u'process_下机日期']] = None
    header[[u'process_送样日期', u'process_上机日期', u'process_下机日期']] = d1[[u'到样日期', u'上机日期',u'下机日期']]
    header[[u'ana_coverage1X(%)',u'ana_coverage20X(%)', u'ana_avg_depth(X)', u'ana_duplication(%)']] = qcsum[[u'1Xcov', u'20Xcov', u'AvgDepth', u'duplication']]
    header[[u'ana_indel(genomic)', u'ana_snp(genomic)', u'ana_ts/tv(genomic)',u'ana_snp(exonic)', u'ana_ts/tv(exonic)']] = qcvcf[[u'indel_genomic', u'snp_genomic', u'ts/tv_genomic',u'snp_exonic', u'ts/tv_exonic']]
    annotation_ver = "3.0"
    header[u'ana_annotation'] = annotation_ver
    header[u'ana_reporter'] = u'吕珊珊'
    #bk_rawdata = u'/anaData/anaData004/children_hos_genetic/rawdata/2018/'+strategy+'/'
    bk_rawdata = u'/anaData/anaData004/children_hos_genetic/rawdata/2019/'+strategy+'/'
    for i in d1[u'公司编码']:
        os.system('ls '+bk_rawdata+i)
    print bk_rawdata
    header[u'bk_rawdata'] = bk_rawdata
    #bk_reports = u'/DATA/BPshare/opm/遗传病-交大附属儿医/2018/'+'/'.join(os.getcwd().split('/')[-2:])
    bk_reports = u'/DATA/BPshare/opm/遗传病-交大附属儿医/2019/'+'/'.join(os.getcwd().split('/')[-2:])
    header[u'bk_reports'] = bk_reports
    #result = total_table.append(header)
    total_table.index = total_table[u'姓名']
    header.index = header[u'姓名']
    for i in header.index:
        if i in total_table.index:
            total_table = total_table.drop(i)
    result = pd.concat([total_table,header])

    ofile = sn.split(sep)[0]+'navicat'+u'）'.join(sn.split(sep)[1].split(')'))
    wb=pd.ExcelWriter(ofile, engine='openpyxl')
    header.to_excel(wb,index=False)
    wb.save()
    ofile1 = m_path+'/doc/navicat.xlsx'
    wb1=pd.ExcelWriter(ofile1, engine='openpyxl')
    result.to_excel(wb1,index=False)
    wb1.save()
    #os.system('cp '+m_path+'/doc/navicat.xlsx '+m_path+'/bk/navicat_bak.xlsx')
    os.system('cp '+m_path+'/doc/navicat.xlsx /DATA/sslyu/Project/Genetics_children_hospital/navicat_bak.xlsx')
    result.to_csv('/DATA/sslyu/Project/Genetics_children_hospital/navicat_bak.csv', index=None, sep="\t", encoding="utf-8")
    return header

def backup(ofile, dirname, sn, sep):
    #ofile.write('sh '+m_path+'/src/unmapped_duplicate.sh list\n')
    ofile.write('cd CNV/\n')
    #ofile.write('cat '+m_path+'/src/CNV_Read.R\n')
    print """
    .libPaths("/home/ana005/R/x86_64-pc-linux-gnu-library/3.4")
    library(ExomeDepth)
    """
    if 'IDT-PANEL' in dirname[-1]:
        print """
        read.table("/home/ana005/data/annovar/bed/IDT/IDT_gene.bed")->bed
        """
    elif 'WES' in dirname[-1]:
        print """
        read.table("/home/ana005/data/annovar/bed/wes/wes_gene.bed")->bed
        """
    else:
        print """
        read.table("/home/ana005/data/annovar/bed/Agilent/S06588914/S06588914_Regions_cnv.bed")->bed
        """
    print """
    colnames(bed)=c("chromosome","start","end","name")
    read.table("list")->bamfile
    unlist(bamfile)->bamfile
    bamfile=as.character(bamfile)
    bamfile <- as.vector(bamfile)
    bam.counts <- getBamCounts(bed.frame=bed, bam.files=bamfile, referenceFasta="/SSD750/PB1/db1/Homo/refseq/hg19.fa")
    save.image("bam.counts.Rdata")
    """
    ofile.write('/DATA/ypliu/opt/R-3.4.3/bin/R # note: choose the right bed file\n\n')
    ofile.write('/DATA/ypliu/opt/R-3.4.3/bin/R --slave --vanilla < '+m_path+'/src/CNV_Run.R > log 2>&1\n')
    ofile.write('cd ..\n')
    snfile = sn.split(sep)[0]+'交付信息表'+u'）'.join(sn.split(sep)[1].split(')'))
    ofile.write('cp *交付信息表* '+dirname[-1]+'\n')
    ofile.write('mv CNV/*csv '+dirname[-1]+'/CNV\n')
    ofile.write('ls -R '+dirname[-1]+'\n\n')
    ofile.write('zip -r '+dirname[-1]+'.zip '+dirname[-1]+'\n')
    #ofile.write('if [ ! -e /DATA/BPshare/opm/遗传病-交大附属儿医/2018/'+dirname[-2]+' ]\n')
    ofile.write('if [ ! -e /DATA/BPshare/opm/遗传病-交大附属儿医/2019/'+dirname[-2]+' ]\n')
    ofile.write('then\n')
    #ofile.write('mkdir -p /DATA/BPshare/opm/遗传病-交大附属儿医/2018/'+dirname[-2]+'\n')
    ofile.write('mkdir -p /DATA/BPshare/opm/遗传病-交大附属儿医/2019/'+dirname[-2]+'\n')
    ofile.write('fi\n')
    #ofile.write('cp -r '+dirname[-1]+' /DATA/BPshare/opm/遗传病-交大附属儿医/2018/'+dirname[-2]+'\n')
    ofile.write('cp -r '+dirname[-1]+' /DATA/BPshare/opm/遗传病-交大附属儿医/2019/'+dirname[-2]+'\n')
    ofile.write('scp list klyang@122.112.248.194:/media/bpdata/data/exon_temp2/list_'+dirname[-1]+'\n')
    ofile.write('ssh kly\n\n')
    return ofile

def create_email(n):
    n = str(n)
    if os.getcwd().startswith('/DATA'):
        dirname = "/".join(os.getcwd().split('/')[-2:])
        dirname = "/DATA/sslyu/Project/Genetics_children_hospital/2018/"+dirname
    else:
        dirname = os.getcwd()
    
    sn = os.getcwd().split('/')[-1].split('_')[1]
    batch = os.getcwd().split('/')[-1].split('_')[0]
    strategy = os.getcwd().split('/')[-1].split('_')[2]

    if sn != "SR":
        if strategy == 'IDT-PANEL':
            sn = sn+'儿医-'+n+'样-IDT4500基因检测分析（'+batch+'批）'
        if strategy == 'WES':
            sn = sn+'儿医-'+n+'样-全外显子检测分析（'+batch+'批）'
        if strategy == 'Agilent-wes':
            sn = sn+'儿医-'+n+'样-安捷伦亚全外显子检测分析（'+batch+'批）'
    else:
        sn = batch+'儿医-'+n+'样-医学全外显子检测分析'
        
    #for f in os.listdir('.'):
    #    if 'xlsx' in f:
    #        if strategy == 'IDT-PANEL':
    #            sn = sn+'儿医-'+n+'样-IDT4500基因检测分析（'+f.split('_')[1].split('.')[0]+'）'
    #            break
    #        if strategy == 'WES':
    #            sn = sn+'儿医-'+n+'样-全外显子检测分析（'+f.split('_')[1].split('.')[0]+'）'
    #            break
    #    if 'eml' in f:
    #        sn = f.split(']')[1].split('.eml')[0]
    
    ofile = open('email.txt','w')
    ofile.write('以下为儿童遗传病分析结果：\n批次：'+sn+'结果\n\n位置\n'+dirname+'\n\n')
    if os.path.exists('peddy_check.log'):
        ofile.write('备注：\n')
        with open('peddy_check.log') as f:
            for l in f:
                ofile.write(l)
        ofile.write('\n')
    if os.path.exists('qcsum_check.log'):
        ofile.write('备注：\n')
        with open('qcsum_check.log') as f:
            for l in f:
                ofile.write(l)
        ofile.write('\n')                
    ofile.write('祝好，\n吕珊珊\n\n')
    ofile.close()

def send_email(choose, n):
    n = str(n)
    import smtplib
    from email.mime.text import MIMEText
    from email.mime.multipart import MIMEMultipart
    from email.mime.base import MIMEBase
    from email import encoders
    dirname = '/'.join(os.getcwd().split('/')[-2:])
    sn = os.getcwd().split('/')[-1].split('_')[1]
    batch = os.getcwd().split('/')[-1].split('_')[0]
    strategy = os.getcwd().split('/')[-1].split('_')[2]

    if sn != "SR":
        if strategy == 'IDT-PANEL':
            sn = sn+'儿医-'+n+'样-IDT4500基因检测分析（'+batch+'批）'
        if strategy == 'WES':
            sn = sn+'儿医-'+n+'样-全外显子检测分析（'+batch+'批）'
        if strategy == 'Agilent-wes':
            sn = sn+'儿医-'+n+'样-安捷伦亚全外显子检测分析（'+batch+'批）'
    else:
        sn = batch+'儿医-'+n+'样-医学全外显子检测分析'

    #for f in os.listdir('.'):
    #    if 'xlsx' in f:
    #        if strategy == 'IDT-PANEL':
    #            sn = sn+'儿医-'+n+'样-IDT4500基因检测分析（'+f.split('_')[1].split('.')[0]+'）'
    #            break
    #        if strategy == 'WES':
    #            sn = sn+'儿医-'+n+'样-全外显子检测分析（'+f.split('_')[1].split('.')[0]+'）'
    #            break
    #    if 'eml' in f:
    #        sn = f.split(']')[1].split('.eml')[0]
    for f in os.listdir('.'):
        if 'zip' in f:
            filename = f
    email_user='sslv@basepair.cn'
    email_send1='lss@sibs.ac.cn'
    email_send2=['yyzhou@basepair.cn', 'jjia@basepair.cn', 'gli@basepair.cn', 'yjsun@basepair.cn'] # in smtplib module, when multi recipients, need a list
    subject = sn+'结果'
    msg = MIMEMultipart()
    msg['From'] = email_user
    if choose == 1:
        msg['To'] = ", ".join(email_send1) # in email module, need a string
    else:
        msg['To'] = ", ".join(email_send2) # in email module, need a string
    msg['Subject'] = subject
    body = open('email.txt').readlines()
    body = ''.join(body)
    msg.attach(MIMEText(body, 'plain', 'utf-8'))
    attachment = open(filename, 'rb')
    part = MIMEBase('application', 'octet-stream')
    part.set_payload((attachment).read())
    encoders.encode_base64(part)
    part.add_header('Content-Disposition', "attachment; filename="+filename)
    msg.attach(part)
    text = msg.as_string()
    mail = smtplib.SMTP('smtp.basepair.cn', 25)
    mail.ehlo()
    mail.starttls()
    mail.login(email_user, 'azby06078961!')
    if choose == 1:
        mail.sendmail(email_user, email_send1, text)
    else:
        mail.sendmail(email_user, email_send2, text)
    mail.quit()

def peddy_check():
    ofile = open('peddy_check.log', 'w')
    sex = pd.read_csv('peddy/pedtrio/results/peddy.sex_check.csv', dtype=str)
    for i in sex.index:
        if sex.loc[i]['error'] == 'True':
            ofile.write(sex.loc[i]['sample_id']+u'分析单上性别为'+sex.loc[i]['ped_sex']+u'，预测性别为'+sex.loc[i]['predicted_sex']+'\n')
            #ofile.write('the sex of '+sex.loc[i]['sample_id']+' is not the same with the predicted result\n')   
    ped = pd.read_csv('peddy/pedtrio/results/peddy.ped_check.csv', dtype=str)
    for i in ped.index:
        if ped.loc[i]['parent_error'] == 'True':
            ofile.write(ped.loc[i]['sample_a']+' and '+ped.loc[i]['sample_b']+' are not child-parent\n')
    ofile.close()
    if os.path.getsize('peddy_check.log') == 0:
        os.remove('peddy_check.log')

def qcsum_check(strategy):
    ofile = open('qcsum_check.log', 'w')
    d2 = pd.read_table('all.qcsum', dtype=str)
    d2.index = d2['sample']
    if strategy == 'PANEL' or strategy == 'Agilent-PANEL' or strategy == 'IDT-PANEL':
        for i in d2.index:
            if float(d2.loc[i]['20Xcov']) < 95:
                ofile.write('the 20X coverage of '+i+' is less than 95%\n')
    else:
        for i in d2.index:
            if float(d2.loc[i]['20Xcov']) < 90:
                ofile.write('the 20X coverage of '+i+' is less than 90%\n')
    ofile.close()
    if os.path.getsize('qcsum_check.log') == 0:
        os.remove('qcsum_check.log')

def generate_jiesuan(navicat,table_map_all,strategy):
    table1 = navicat[[u'姓名', 'exp_reads', 'exp_bases(Mb)', 'exp_Q30', 'ana_coverage1X(%)', 'ana_coverage20X(%)', 'ana_avg_depth(X)', 'ana_indel(genomic)', 'ana_snp(genomic)', 'ana_ts/tv(genomic)', 'ana_snp(exonic)', 'ana_ts/tv(exonic)', u'process_送样日期', u'申请科室']]
    table1.columns = ['Sample', 'Total Reads', 'Total Bases(Mb)', 'Q30', '1Xcov', '20Xcov', 'Ave Depth', 'Indel_genomic', 'Snp_genomic', 'ts/tv_genomic', 'snp_exonic', 'ts/tv_exonic',u'送检日期', u'送检科室']
    index = []
    for i in range(navicat.shape[0]):
        if navicat[u'样本ID'].isnull()[i]:
            index.append(navicat[u'样本编号'][i])
        else:
            index.append(navicat[u'样本ID'][i])
    table1.index = index
    table2 = pd.concat([table1,table_map_all], axis=1)
    table2.insert(table2.columns.tolist().index('Unmapped reads'), 'Valid reads', table2['Total Reads'])
    print table2
    mapping_ratio = table2.apply(lambda x:round((int(x['Valid reads'])-int(x['Unmapped reads']))/float(x['Valid reads'])*100,2), axis=1)
    table2.insert(table2.columns.tolist().index('Unmapped reads')+1, 'Mapping Ratio', [str(i)+'%' for i in mapping_ratio])
    tmp_list = ['Valid reads','Unmapped reads','Mapping Ratio','UNPAIRED_READ_DUPLICATES','READ_PAIR_DUPLICATES','Percent duplication']
    tmp = table2[tmp_list]
    table2 = table2.drop(tmp_list, axis=1)
    for i in range(len(tmp_list)):
        table2.insert(table2.columns.tolist().index('Q30')+1+i, tmp_list[i], tmp[tmp_list[i]])
    table2.insert(table2.columns.tolist().index(u'送检日期'), u'序号', range(1,table2.shape[0]+1))
    if strategy == 'WES':
        table2[u'检测项目'] = u'全外显子测序'
        table2[u'金额/元'] = u'2970'
    else:
        table2[u'检测项目'] = u'IDT-PANEL'
        table2[u'金额/元'] = u'0'
    table2.to_csv('jiesuan_table.csv', sep="\t")
    ofile = 'jiesuan_table.xlsx'
    wb = pd.ExcelWriter(ofile,engine='openpyxl')
    table2.to_excel(wb,index=False)
    wb.save()

def exon_plot(sample):
    bed = pd.read_table('/DATA/sslyu/Project/Genetics_children_hospital/all_ch_freq_gene.bed', header=None)
    bed.index = bed[3]
    gene = set([i.split('_exon')[0] for i in bed[3]])
    cdslen = bed.apply(lambda x:x[2]-x[1], axis=1)
    loci = pd.read_table('/DATA/sslyu/Project/Genetics_children_hospital/all_ch_freq_gene.loci', header=None)
    depth = pd.read_table('depth/'+sample+'.depth', header=None)
    depth[1] = depth[1].apply(str)
    depth['loci'] = depth[0]+'-'+depth[1]
    loci[1] = loci[1].apply(str)
    loci['loci'] = loci[0]+'-'+loci[1]
    d = pd.merge(depth, loci, how='left', on='loci')
    count = d.groupby('2_y')['2_x'].sum()
    mean_depth = count.divide(cdslen)
    mean_depth = mean_depth.apply(lambda x:round(x,2))
    df = pd.DataFrame({'exons':mean_depth.index,sample:mean_depth}, columns=['exons',sample])
    df.to_csv('tsv/'+sample+'_exon.depth.tsv', index=None, sep="\t")

def submit(sample,script):
    if script != "trio_cmd.sh":
        os.chdir(sample)
        p1 = subprocess.Popen(["sbatch", script], stdout=subprocess.PIPE)
        os.chdir('..')
    else:
        os.chdir('trio/'+sample)
        p1 = subprocess.Popen(["sbatch", script], stdout=subprocess.PIPE)
        os.chdir('../..')
    jobID = p1.communicate()[0].strip('\n').split()[-1]
    return jobID

def check_sacct(jobID_l):
    n = 0
    for i in jobID_l:
        p1 = subprocess.Popen(["sacct", "-u", "sslyu"], stdout=subprocess.PIPE)
        p2 = subprocess.Popen(["grep",i], stdin=p1.stdout, stdout=subprocess.PIPE)
        p3 = subprocess.Popen(["grep","RUNNING"], stdin=p2.stdout, stdout=subprocess.PIPE)
        sacct = p3.communicate()[0]
        #print sacct+' is'
        if sacct == '':
            n = n+1
    return n

def check_slurm(jobID_d, script):
    n = 0
    if script != "trio_cmd.sh":
        for i in jobID_d:
            p1 = subprocess.Popen(["grep", "-ci", "error", i+"/slurm-"+jobID_d[i]+".out"], stdout=subprocess.PIPE)
            sacct = p1.communicate()[0]
            if sacct.split('\n')[0] == '0':
                n = n+1
            else:
                print i+" "+script+" has a problem"
                sys.exit()
    else:
        for i in jobID_d:
            p1 = subprocess.Popen(["grep", "-ci", "error", "trio/"+i+"/slurm-"+jobID_d[i]+".out"], stdout=subprocess.PIPE)
            sacct = p1.communicate()[0]
            if sacct.split('\n')[0] == '0':
                n = n+1
            else:
                print i+" "+script+" has a problem"
                sys.exit()
    return n            

def cmd_trio_submit(d1):
    pedigree = generate_pedigree(d1)
    cmd_trio_jobID = {}
    for k in pedigree:
        cmd_trio_jobID[k] = {}
        for v in pedigree[k]:
            print v
            sample_v = d1.loc[v]['sample']
            jobID = submit(sample_v, 'cmd.sh')
            cmd_trio_jobID[k][sample_v] = jobID
    return cmd_trio_jobID

def cmd_trio_check(cmd_trio_jobID, d1):
    pedigree = generate_pedigree(d1)
    triocmd_jobID = {}
    while True:
        time.sleep(60)
        for k in cmd_trio_jobID.keys():
            n_sacct_cmd = check_sacct(cmd_trio_jobID[k].values())
            if n_sacct_cmd == len(cmd_trio_jobID[k]):
                #list_trio = []
                #for v in pedigree[k]:
                #    print v
                #    list_trio.append(d1.loc[v]['sample'])
                n_slurm_cmd = check_slurm(cmd_trio_jobID[k], 'cmd.sh')
                if n_slurm_cmd == len(cmd_trio_jobID[k]):
                    sample = d1.loc[k]['sample']
                    jobID = submit(sample, 'trio_cmd.sh')
                    triocmd_jobID[sample] = jobID
                    cmd_trio_jobID.pop(k)
        if len(cmd_trio_jobID) == 0:
            break
    return triocmd_jobID            

def copy_check(jobID_d):
    copy_sacct = False
    copy_slurm = False
    while True:
        time.sleep(60)
        #print jobID_l
        n_sacct = check_sacct(jobID_d.values())
        #print "copy_check n_sacct %s" % n_sacct
        if n_sacct == len(jobID_d):
            copy_sacct = True
            break
    while copy_sacct:
        n_slurm = check_slurm(jobID_d, 'cp_rawdata.sh')
        if n_slurm == len(jobID_d):
            copy_slurm = True
            break
    return copy_sacct and copy_slurm        

#def cmd_trio_check(trio_jobID):
#    trio_list = pd.read_table('list_trio', header=None)[0].tolist()
#    cmd_sacct_trio = False
#    cmd_slurm_trio = False
#    while True:
#        time.sleep(60)
#        n_sacct_cmd = check_sacct(trio_jobID)
#        if n_sacct_cmd == len(trio_jobID):
#            cmd_sacct_trio = True
#            break
#    while cmd_sacct_trio:
#        n_slurm_cmd = check_slurm(trio_list, trio_jobID, 'cmd.sh')
#        if n_slurm_cmd == len(trio_jobID):
#            cmd_slurm_trio = True
#            break
#    return cmd_sacct_trio and cmd_slurm_trio

def trio_check(jobID_trio):
    #list_trio = pd.read_table('trio/list', header=None)[0].tolist()
    trio_sacct = False
    trio_slurm = False
    while True:
        time.sleep(60)
        jobID_trio_l = jobID_trio.values()
        n_sacct_trio = check_sacct(jobID_trio_l)
        if n_sacct_trio == len(jobID_trio_l):
            trio_sacct = True
            break
    while trio_sacct:
        n_slurm_trio = check_slurm(jobID_trio, 'trio_cmd.sh')
        if n_slurm_trio == len(jobID_trio):
            trio_slurm = True    
            break
    return trio_sacct and trio_slurm

def cmd_single_check(single_jobID):
    cmd_sacct_single = False
    cmd_slurm_single = False
    #list_single = pd.read_table('list_single', header=None)[0].tolist()
    while True:
        time.sleep(60)
        n_sacct_single = check_sacct(single_jobID.values())
        if n_sacct_single == len(single_jobID):
            cmd_sacct_single = True
            break
    while cmd_sacct_single:
        n_slurm_single = check_slurm(single_jobID, 'cmd.sh')
        if n_slurm_single == len(single_jobID):
            cmd_slurm_single = True
            break
    return cmd_sacct_single and cmd_slurm_single

