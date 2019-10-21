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
pd.set_option('mode.chained_assignment', None)
class MyParser(argparse.ArgumentParser):
   def error(self, message):
      sys.stderr.write('error: Check your parameters!\n')
      sys.stderr.write('error: %s\n' % message)
      self.print_help()
      sys.exit(2)

def init():
    example_text = '''Example:

    python merge_and_norm_cov.py -i ./mapping_stats/ -s samples_mapnum.txt -m scaf2bin.txt

    '''
    parser = argparse.ArgumentParser(description='Merge and normalize coverage with TPM \n',epilog=example_text, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser._optionals.title = 'Mandatory Arguments'
    parser.add_argument('-i', required=True, help='Path of your mapping statistics folder')
    parser.add_argument('-s', required=True, help='Sample mapping read number table')
    parser.add_argument('-m', required=True, help='Mapping information of Contig-MAG table')
    return(parser)


def main():
    parser = init()
    args = parser.parse_args()
    # get parameters
    sample_readstats = args.s
    scaf2bin = args.m
    mapping_stats_dir = args.i +'/'
    # parse sample mapped reads
    sample_ids = []
    sample_reads = dict()
    for inp in open(sample_readstats):
        data = re.split(r'\s+', inp.strip())
        sample_ids.append(data[0])
        sample_reads[data[0]] = float(data[1])
    # parse scaf2bin
    scaf2bin_df = pd.read_csv(scaf2bin,sep='\t')
    scaf2bin_df.columns = ['Contig_ID','Bin_ID']
    # save info
    all_cov_df_wb = scaf2bin_df
    all_cov_norm_df_wb = scaf2bin_df
    # parse each sample
    all_cov_df = pd.DataFrame()
    for sample_id in sample_ids:
        sample_id = sample_id.strip()
        # get coverage file
        cov_df_fn = ''.join([mapping_stats_dir, sample_id, '.covstats.txt'])
        cov_df = pd.read_csv(cov_df_fn, sep='\t')
        cov_df['Contig_ID'] = cov_df.iloc[:,0].str.split(' ', expand=True).iloc[:,0]
        # save contig id and average fold, and renamed with sample id
        sel_cols = ['Contig_ID', 'Avg_fold']
        sel_cov_df = cov_df[sel_cols]
        sel_cov_df_cp = sel_cov_df
        sel_cov_df.columns = ['Contig_ID', ''.join([sample_id,'.raw_fold'])]
        # normlized by 1 million scale factor
        millions = 1000000
        scale_factor = sample_reads[sample_id]/millions
        sel_cols = ['Contig_ID']
        sel_cov_norm_df = cov_df[sel_cols]
        #cov_df2= cov_df.copy()
        sel_cov_norm_df['Norm_fold'] = cov_df['Avg_fold']/scale_factor
        #sel_cov_norm_df = sel_cov_norm_df.drop('Avg_fold')
        sel_cov_norm_df.columns = ['Contig_ID', ''.join([sample_id,'.norm_fold'])]
        # merge and save
        all_cov_df_wb = all_cov_df_wb.merge(sel_cov_df, on = 'Contig_ID',how = 'inner')
        all_cov_norm_df_wb = all_cov_norm_df_wb.merge(sel_cov_norm_df, on = 'Contig_ID',how = 'inner')
    # save coverage table with bin ids
    all_cov_df_wb.to_csv('ctg_bin_cov.csv', index=False, sep='\t')
    all_cov_norm_df_wb.to_csv('ctg_bin_cov.norm.csv', index=False, sep='\t')


if __name__ == "__main__":
   main()
