#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import igraph
import pandas
from reactome_library import OntoGraph

###############################################################################
# Description:
###############################################################################
# This script does the same analysis are the previous one, but instead of
# calculating the scores for MASLD netwrok nodes, it does it for the enriched
# Reactome pathways, summing up the scores of the respecitve gene sets.
# Outputs:
# - A tsv file with all the scores
###############################################################################


if __name__ == "__main__":
    
    main_dir = '../../results/ucam_sanyal/networks/'
    output_dir = main_dir+'network_analysis/'
    
    filename = output_dir+'reactome_enriched_graph.pkl'
    reactome_graph = igraph.Graph.Read_Pickle(filename)
    
    scores_df = pandas.read_csv(output_dir+'final_scores/network_nodes_scores_per_sw.tsv',
                                sep='\t', index_col=0)
    
    selected_terms = []
    for term in reactome_graph.vs:
        if term['plot'] == True:
            selected_terms.append(term['name'])
            
    absolute_dfs = []
    cumsum_dfs = []
    for term in selected_terms:
        genes = reactome_graph.vs.find(name=term)['annotation']
        centralities_df = reactome_graph.vs.find(name=term)['centralities']
        gene_scores_df = scores_df.loc[genes,:]
        term_dfs = []
        for sw in gene_scores_df.columns:
            sw_df = gene_scores_df.loc[:, [sw]]
            sw_df = pandas.merge(sw_df, centralities_df, left_index=True, right_index=True)
            sw_df = (sw_df.iloc[:,0]*sw_df.iloc[:,1]).to_frame(sw)
            term_dfs.append(sw_df)
        term_df = pandas.concat(term_dfs, axis=1).sum(axis=0).to_frame(term)#.cumsum()
        absolute_dfs.append(term_df)
        term_df = pandas.concat(term_dfs, axis=1).sum(axis=0).cumsum().to_frame(term)
        cumsum_dfs.append(term_df)
        
    absolute_df = pandas.concat(absolute_dfs, axis=1).transpose()
    rownames = []
    for term in absolute_df.index:
        rownames.append(reactome_graph.vs.find(name=term)['definition'])
    absolute_df.index = rownames
    
    cumsum_df = pandas.concat(cumsum_dfs, axis=1).transpose()
    rownames = []
    for term in cumsum_df.index:
        rownames.append(reactome_graph.vs.find(name=term)['definition'])
    cumsum_df.index = rownames
    
    cumsum_df.to_csv(output_dir+'final_scores/reactome_cumsum_values.tsv', sep='\t')
    absolute_df.to_csv(output_dir+'final_scores/reactome_absolute_values.tsv', sep='\t')