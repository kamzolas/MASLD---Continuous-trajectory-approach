#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import igraph
import pandas
import os
import numpy


if __name__ == '__main__':
    
    output_dir = '../../data/sc_liver_cell_atlas/filtered_cell_type_markers/'
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        
    network_df = pandas.read_csv('../../results/networks/MASLD_unified_undirected_network.csv', sep=',')
    network = igraph.Graph.DataFrame(network_df, directed=False, use_vids=False)
    network_nodes = network.vs['name']

    files = os.listdir('../../data/sc_liver_cell_atlas/cell_type_markers/')
    for f in files:
        df = pandas.read_csv('../../data/sc_liver_cell_atlas/cell_type_markers/'+f, sep=',', index_col=0)
        step = 0.1
        logfc = numpy.log2(1.5)
        df = pandas.read_csv('../../data/sc_liver_cell_atlas/cell_type_markers/'+f, sep=',', index_col=0)
        df = df.loc[df.gene.isin(network_nodes),:]
        df = df.loc[df.gene.apply(lambda x: not x.startswith('RPL')),:]
        df = df.loc[df.gene.apply(lambda x: not x.startswith('RPS')),:]
        df = df.loc[df.gene.apply(lambda x: not x.startswith('RPN')),:]
        df = df.loc[df.gene.apply(lambda x: not x.startswith('MT')),:]
        df = df.loc[df.gene.apply(lambda x: not x.startswith('MRPL')),:]
        df = df.loc[df.gene.apply(lambda x: not x.startswith('MRPS')),:]
        df = df.loc[df.loc[:, 'avg_log2FC'] > logfc,]
        while df.shape[0] > 700:
            logfc += step
            df = df.loc[df.loc[:, 'avg_log2FC'] > logfc,]
        cell_type = f.split('.')[0]
        df = df.loc[:,['gene', 'avg_log2FC']]
        df.sort_values(by='avg_log2FC', ascending=False, inplace=True)
        df.reset_index(inplace=True, drop=True)
        df.to_csv(output_dir+cell_type+'.tsv', sep='\t')
    