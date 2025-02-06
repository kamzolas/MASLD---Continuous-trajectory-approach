#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas
import sys

if __name__ == '__main__':
    
    deep_split = 4
    min_size = 180
    
    main_dir = '../../results/networks/'
    nodes_df = pandas.read_csv(main_dir+'MASLD_nodes_without_semantics.csv', 
                               sep=',', index_col=0)
    entrez_df = pandas.read_csv(main_dir+'valid_entrez_ids.csv', 
                                sep=',', index_col=0, dtype=str)
    
    mapping_df = pandas.read_csv('../../data/gene_symbols_to_entrez.tsv',
                                 sep='\t', dtype=str)
    to_keep = mapping_df.loc[~pandas.isnull(mapping_df.gene_symbol)].index.tolist()
    mapping_df = mapping_df.loc[to_keep,:]
    to_keep = mapping_df.loc[~pandas.isnull(mapping_df.entrez_id)].index.tolist()
    mapping_df = mapping_df.loc[to_keep,:]
    
    # First merging to get the entrez ids
    merged_df = pandas.merge(nodes_df, mapping_df, right_on='gene_symbol', 
                             left_on='gene_symbol', how='inner')
    
    # Remove the entries with nan in entrez ids field
    to_keep = merged_df.loc[~pandas.isnull(merged_df.entrez_id)].index.tolist()
    merged_df = merged_df.loc[to_keep,:]
    
    # Second merging to keep the valid entrez ids
    merged_df = pandas.merge(merged_df, entrez_df, right_on='entrez_id', 
                             left_on='entrez_id', how='inner')

    # Remove any remaining duplicated entries    
    to_remove = merged_df.loc[merged_df.duplicated(subset='gene_symbol', keep='last')].index
    to_keep = list(filter(lambda i: not i in to_remove, merged_df.index))
    
    merged_df = merged_df.loc[to_keep,:]
    merged_df.reset_index(drop=True, inplace=True)
    
    merged_df.to_csv(main_dir+'MASLD_nodes_without_semantics_filtered.csv', sep=',')