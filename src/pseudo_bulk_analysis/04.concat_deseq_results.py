#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas
import numpy
import os
import operator



if __name__ == "__main__":
    
    data_dir = '../../results/pseudo_bulk_analysis/deseq_outputs/'
        
    batches = [int(i) for i in os.listdir(data_dir)]
    batches = sorted(batches)
    batches = [str(i) for i in batches]
    
    
    comparisons = {}
    for batch in batches:
        files = os.listdir(data_dir+batch+'/')
        for f in files:
            comparison = f.split('.')[0]
            deseq_results = pandas.read_csv(data_dir+batch+'/'+f, sep=' ')
            deseq_results = deseq_results.loc[:, ['gene_symbol', 'log2FoldChange', 'padj']]
            deseq_results.set_index('gene_symbol', inplace=True)
            comparisons.setdefault(comparison, []).append(deseq_results)
            
    final_results = []
    
    keys = [(int(s.split('_')[1]),s) for s in comparisons.keys()]
    keys = sorted(keys, key=operator.itemgetter(0), reverse=False)
    keys = [c[1] for c in keys]
    
    for comparison in keys:
        
        dfs = comparisons[comparison]
        df = pandas.concat(dfs, axis=1, ignore_index=False)
        
        padj_df = df.loc[:, 'padj']
        padj_df.fillna(1, inplace=True)
        padj_df = -numpy.log10(padj_df)
        padj_df = padj_df.mean(axis=1).to_frame('padj')
        padj_df.loc[:,'padj'] = numpy.power(10,-padj_df.loc[:,'padj'])
        
        logfc_df = df.loc[:, 'log2FoldChange']
        logfc_df.fillna(0, inplace=True)
        logfc_df = logfc_df.mean(axis=1).to_frame('log2fc')
        
        sws = comparison.split('_vs_')
        sws[0] = sws[0].replace('_','')
        sws[1] = sws[1].replace('_','')
        key = '_vs_'.join(sws)
        statistics_df = pandas.concat([logfc_df, padj_df], axis=1, ignore_index=False)
        statistics_df.columns = [key+'_log2FoldChange', key+'_padj']
        
        final_results.append(statistics_df)
        
    df = pandas.concat(final_results, axis=1)
    df.loc[:,'external_gene_name'] = df.index.tolist()
    df = df.reset_index(drop=True)
    df.to_csv('../../results/pseudo_bulk_analysis/deseq_results_SW_vs_previous.csv', sep=',')