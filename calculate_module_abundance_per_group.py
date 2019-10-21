#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
__author__ = 'Yaxin Xue'
__license__ = "GPL"
__version__ = "1.0.0"
__email__ = "yaxin.xue@uib.no"

"""
import sys, os, re
import argparse
import pandas as pd
import numpy as np

# define functions
def parse_mo_info(km_file):
    km_file = open(km_file)
    mo_info = dict()
    mo_map_ko = dict()
    ko_map_mo = dict()
    ko_info = dict()
    current_mo = ''
    for inp in km_file:
        inp = inp.strip()
        # parse module info
        if inp.startswith('D'):
            data = re.split('\s\s+', inp)
            mo_id = data[1]
            mo_name = data[2]
            if mo_id not in mo_info:
                mo_info[mo_id] = mo_name
                current_mo = mo_id
                mo_map_ko[mo_id] = []
        # parse kegg info
        if inp.startswith('E'):
            data = re.split('\s\s+', inp)
            ko_id = data[1]
            ko_name = data[2]
            mo_map_ko[current_mo].append(ko_id)
            ko_info[ko_id] = ko_name
            if ko_id not in ko_map_mo:
                ko_map_mo[ko_id] = []
                ko_map_mo[ko_id].append(current_mo)
            else:
                ko_map_mo[ko_id].append(current_mo)
    km_file.close()
    return(mo_info, mo_map_ko, ko_info, ko_map_mo)

def convert_gene_to_contig_map(gene_map_ko):
    gene_map_ko = pd.read_csv(gene_map_ko, sep = '\t', header=None)
    gene_map_ko.columns = ['gene_id', 'kegg_id']
    # convert gene id to contig id
    gene_ids = gene_map_ko.gene_id.str.split('_', expand=True)
    ctg_map_ko = gene_map_ko
    ctg_map_ko['Contig_ID'] = gene_ids.iloc[:,0].astype(str)+'_'+gene_ids.iloc[:,1]
    ctg_map_ko = ctg_map_ko.dropna()
    return(ctg_map_ko)

def detect_module_by_ko_overlap(ctg_map_ko, mo_map_ko):
    ctg_has_ko = dict()
    for idx,row in ctg_map_ko.iterrows():
        ctg_id = row['Contig_ID']
        kegg_id = row['kegg_id']
        if ctg_id not in ctg_has_ko:
            ctg_has_ko[ctg_id] = []
            ctg_has_ko[ctg_id].append(kegg_id)
        else:
            ctg_has_ko[ctg_id].append(kegg_id)
    list_ctg_has_mo = []
    for key,value in ctg_has_ko.items():
        ctg_id = key
        ko_sets = set(value)
        for key,value in mo_map_ko.items():
            mo_id = key
            mo_set = set(value)
            if list(ko_sets & mo_set):
                data =[ctg_id, mo_id]
                list_ctg_has_mo.append(data)
    # transfer to data frame
    ctg_has_mo = pd.DataFrame(list_ctg_has_mo)
    ctg_has_mo.columns = ['Contig_ID', 'module_id']
    return(ctg_has_mo)


def cal_bin_mo_ab_divnum(group_ctg_wg_map_ko, mo_map_ko):
    group_bin_id = group_ctg_wg_map_ko.Bin_ID.unique()
    list_bin_mo_cov = []
    for bin_id in group_bin_id:
        bin_group_ctg_wg_map_ko = group_ctg_wg_map_ko.loc[group_ctg_wg_map_ko.Bin_ID==bin_id]
        bin_ko_cov = bin_group_ctg_wg_map_ko.loc[:,['kegg_id', 'cov_sum']]
        # calculate average ko coverage
        bin_ko_cov_avg = bin_ko_cov.groupby('kegg_id', as_index=False)['cov_sum'].mean()
        # check modules with each bin: each bin shares kos
        bin_ko_set = set(bin_ko_cov.kegg_id)
        for key,value in mo_map_ko.items():
            mo_id = key
            mo_set = set(value)
            mo_ko_num = len(mo_set)
            #if mo_set.issubset(bin_ko_set):
            if list(mo_set & bin_ko_set):
                ko_data = list(mo_set & bin_ko_set)
                mo_cov_sum = float(0)
                for ko_id in ko_data:
                    mo_cov_sum = mo_cov_sum+float(bin_ko_cov_avg.loc[bin_ko_cov_avg.kegg_id==ko_id].cov_sum)
                mo_cov_avg = mo_cov_sum/mo_ko_num
                data = [bin_id, mo_id, mo_cov_avg]
                list_bin_mo_cov.append(data)
    group_bin_mo_ab = pd.DataFrame(list_bin_mo_cov)
    group_bin_mo_ab.columns= ['Bin_ID', 'module_id', 'abundance']
    return(group_bin_mo_ab)

def init():
    example_text = '''Example:

    python calculate_module_abundance_per_group.py -gk gene_ko_anno_ghostkoala.txt -ko ko00002.keg -cg NORM_group.csv

    '''
    parser = argparse.ArgumentParser(description='Calculate module abundance per group with a MAG-centric view\n',epilog=example_text, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser._optionals.title = 'Mandatory Arguments'
    parser.add_argument('-gk', required=True, help='Gene-KO annotation table')
    parser.add_argument('-ko', required=True, help='KEGG module database')
    parser.add_argument('-cg', required=True, help='Coverage table with grouping informtation')
    return(parser)

def main():
    # parse arguments
    parser = init()
    args = parser.parse_args()
    kegg_module_db = args.ko
    gene_map_ko = args.gk
    ctg_wg = pd.read_csv(args.cg)
    ##### Start analysis #####

    #parse KEGG module database
    (mo_info, mo_map_ko, ko_info, ko_map_mo) = parse_mo_info(kegg_module_db)
    mo_info_pd = pd.DataFrame.from_dict(mo_info, orient='index')
    mo_info_pd = mo_info_pd.reset_index()
    mo_info_pd.columns = ['module_id','module_anno']
    #convert gene-map-ko to contig-map-ko
    ctg_map_ko = convert_gene_to_contig_map(gene_map_ko)
    #detect modules based on KO overlap
    ctg_has_mo = detect_module_by_ko_overlap(ctg_map_ko, mo_map_ko)
    # build df for all detected modules
    all_mo_ab = pd.DataFrame(ctg_has_mo['module_id'].unique())
    all_mo_ab.columns = ['module_id']
    # get all group ids
    groups = np.sort(ctg_wg.group.unique())
    # calculate module abundane for each group
    for grp in groups:
        # select group
        group_ctg_wg = ctg_wg[ctg_wg.group==grp]
        # merge KO info
        group_ctg_wg_map_ko = group_ctg_wg.merge(ctg_map_ko, how='inner', on ='Contig_ID')
        # sum all samples' cov
        group_ctg_wg_map_ko['cov_sum'] = group_ctg_wg_map_ko.drop(['Contig_ID','Bin_ID','group'], axis = 1).sum(axis=1)
        # calculate module abundance in each bin
        group_bin_mo_ab = cal_bin_mo_ab_divnum(group_ctg_wg_map_ko, mo_map_ko)
        # calculate module abundance in eech group
        group_bin_mo_ab.groupby('module_id', as_index=False)
        group_mo_ab = group_bin_mo_ab.groupby('module_id', as_index=False)['abundance'].sum()
        group_mo_ab.columns = ['module_id',grp]
        # save group abundance
        all_mo_ab = all_mo_ab.merge(group_mo_ab, how='left', on='module_id')

    #add module annotation
    all_mo_ab = all_mo_ab.merge(mo_info_pd, how='left', on = 'module_id')
    #output final module abundance results
    all_mo_ab.to_csv('module_abundance_per_group.csv',index=False)

if __name__ == "__main__":
   main()
