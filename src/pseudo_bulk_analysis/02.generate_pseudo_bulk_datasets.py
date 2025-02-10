#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import pandas
import os

if __name__ == "__main__":
    
    data_dir = '../../results/pseudo_bulk_analysis/pseudo_bulk_RNAseq/'
    output_dir = '../../results/pseudo_bulk_analysis/pseudo_bulk_RNAseq_datasets/'
    
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        
    batches = os.listdir(data_dir)
    
    for batch in batches:
        dfs = []
        files = os.listdir(data_dir+batch+'/')
        sws = list(range(1, len(files)+1))
        for sw in sws:
            df = pandas.read_csv(data_dir+batch+'/'+'SW_'+str(sw)+'_raw_counts.tsv', sep='\t', index_col=0)
            dfs.append(df)
        df = pandas.concat(dfs, axis=1, ignore_index=False)
        filename = output_dir+'raw_counts_'+str(batch)+'.tsv'
        df.to_csv(filename, sep='\t')
        
        metadata = pandas.DataFrame(data={'sample':df.columns.tolist()})
        metadata.loc[:,'condition'] = metadata.loc[:,'sample'].apply(lambda x: x.split('_sample')[0])
        filename = output_dir+'metadata_'+str(batch)+'.tsv'
        metadata.to_csv(filename, sep='\t')