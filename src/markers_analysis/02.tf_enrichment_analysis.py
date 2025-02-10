#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas
import igraph
import numpy
import json
from scipy.stats import fisher_exact
from statsmodels.stats import multitest



def tf_enrichment(sample):

    stats_df = []
    for tf, regulon in collect_tri_mapping.items():

        M = len(collect_tri_background)
        n = len(regulon)
        N = len(sample)
        regulon_in_sample = set(regulon).intersection(sample)
        x = len(regulon_in_sample)

        table = numpy.array([[x, n-x],[N-x, M-(n+N)+x]])
        pvalue = fisher_exact(table, alternative='greater').pvalue
        # = fisher_exact(table, alternative='greater').statistic
        stats_df.append([tf, pvalue])

    stats_df = pandas.DataFrame(stats_df)
    stats_df.sort_values(by=1, ascending=True, inplace=True)
    adjusted_pvalues = multitest.multipletests(stats_df.iloc[:,1], method='fdr_bh')[1]
    stats_df.loc[:,'fdr'] = adjusted_pvalues
    stats_df.columns = ['tf', 'pvalue', 'fdr']
    
    return stats_df



if __name__ == '__main__':
    
    
    main_dir = '../../results/networks/network_analysis/'
    output_dir = main_dir+'24_markers/'
    
    # Network
    network_df = pandas.read_csv('../../results/networks/MASLD_unified_directed_network_with_semantics.tsv', sep='\t')
    MASLD_network = igraph.Graph.DataFrame(network_df, directed=True, use_vids=False)
    network_nodes = MASLD_network.vs['name']
    
    # Collect tri database
    collect_tri = pandas.read_csv('../../data/collectTRI_network.tsv', sep='\t')
    collect_tri = collect_tri.loc[:, ['tf', 'target', 'mor']]
    
    filtered_collect_tri = collect_tri.loc[collect_tri.target.isin(network_nodes),:]
    filtered_collect_tri = filtered_collect_tri.loc[filtered_collect_tri.tf.isin(network_nodes),:]
    collect_tri_background = filtered_collect_tri.target.unique().tolist()

    collect_tri_mapping = {}
    groups = filtered_collect_tri.groupby(by='tf')
    for tf, tmp_df in groups:
        if tf in network_nodes:
            collect_tri_mapping.update({tf:tmp_df.target.unique().tolist()})
    
    reference_df = pandas.read_csv(output_dir+'24_plasma_markers.tsv',
                                   sep='\t', index_col=0)
    
    output = {}
    for index, row in reference_df.iterrows():
        sample = row['genes'].split(';')
        stats_df = tf_enrichment(sample)
        tfs = stats_df.loc[stats_df.fdr<0.1,'tf'].tolist()
        regulons_df = filtered_collect_tri.loc[(filtered_collect_tri.tf.isin(tfs)) & (filtered_collect_tri.target.isin(sample)), :]
        regulons = [list(a) for a in regulons_df.values]
        d = {'genes': sample, 'tfs':tfs, 'regulons':regulons}
        output.update({row['gene'] +'_'+ row['sw'] +'_'+ row['direction']: d})
            
        
    filename = output_dir+'24_plasma_markers_extended.json'
    with open(filename, 'w') as fobj:
        json.dump(output, fobj)