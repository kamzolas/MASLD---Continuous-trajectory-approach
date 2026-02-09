#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas # 2.1.4
import igraph # 0.11.4
import numpy # 1.26.4
from scipy.stats import fisher_exact # 1.11.4
from statsmodels.stats import multitest # 0.14.1
import sys
import os
sys.path.append('../07.network_analysis/')
from reactome_library import OntoGraph
 
###############################################################################
# Description:
###############################################################################
# This script compares the results from real and bulk data analysis and 
# identifies cell types enriched in de-regulated pathways. The MASLD-related 
# Reactome pathways, as well as their gene sets are retrieved from the 
# enriched_graph object. These pathways describe the differentiated 
# process in MASLD. This differentiation could have been produced either by 
# changes in cell type proportions or actual de-regulation. In order to untangle
# this, genes present in both real and pseudo signatures, whose activation
# score in the pseudo-signature exceeds that in thereal signature, are classified
# as composition-driven and removed from the MASLD pathways. In this way, the 
# remaining gene sets of MASLD-related pathways contain genes which 
# differentiation can be explained only by the deregulation of the respective 
# processes. An enrichment analysis between these filtered gene sets and
# cell type marker lists reveals those cell types which are enriched in the 
# MASLD de-regulated pathways. 
# Outputs:
# - Tsv files with pathways scores for the remaining gene sets in the MASLD-
# related pathways.
# - Tsv file with the enrichemnt analysis results
###############################################################################


if __name__ == '__main__':
    
    ###########################################################################
    # Inputs
    ###########################################################################
    main_dir = '../../results/ucam_sanyal/networks/'
    real_results_dir = '../../results/ucam_sanyal/networks/network_analysis/'
    pseudo_results_dir = '../../results/ucam_sanyal/pseudo_bulk_analysis/networks/network_analysis/'
    
    filename = '../../results/ucam_sanyal/networks/network_analysis/reactome_enriched_graph.pkl'
    reactome_graph = igraph.Graph.Read_Pickle(filename)
    
    real_scores_df = pandas.read_csv(real_results_dir+'final_scores/network_nodes_scores_per_sw.tsv',
                                 sep='\t', index_col=0)
    pseudo_scores_df = pandas.read_csv(pseudo_results_dir+'final_scores/network_nodes_scores_per_sw.tsv',
                                       sep='\t', index_col=0)
    
    
    ###########################################################################
    # 1. Select the pathways for which the signature should be compared and do 
    # the comparisons. Save the new pathway scores from which the pseudo bulk-
    # derived scores have been removed.
    ###########################################################################    
    selected_terms = []
    for term in reactome_graph.vs:
        if term['plot'] == True:
            selected_terms.append(term['name'])
            
    absolute_dfs = []
    cumsum_dfs = []
    genes_dictionary = {}
    for term in selected_terms:
        genes = reactome_graph.vs.find(name=term)['annotation']
        centralities_df = reactome_graph.vs.find(name=term)['centralities']
        real_gene_scores_df = real_scores_df.loc[genes,:]
        pseudo_gene_scores_df = pseudo_scores_df.loc[genes,:]
        term_dfs = []
        for sw in real_gene_scores_df.columns:
    
            real_sw_df = real_gene_scores_df.loc[:, [sw]]
            pseudo_sw_df = pseudo_gene_scores_df.loc[:, [sw]]
    
            comparison_df = pandas.merge(real_sw_df, pseudo_sw_df, left_index=True, 
                                         right_index=True, suffixes=['_real', '_pseudo'])
            term_genes_in_real = comparison_df.loc[comparison_df.loc[:,sw+'_real'] !=0 ,:].index.tolist()[:]
            
            new_col = numpy.sign(comparison_df.iloc[:,0])*numpy.sign(comparison_df.iloc[:,1])
            comparison_df.loc[:, 'direction'] = new_col
            comparison_df.loc[comparison_df.loc[:, 'direction'] == -1, sw+'_pseudo'] = 0
            new_col = (comparison_df.iloc[:,1].abs()>0) & (comparison_df.iloc[:,1].abs() >= comparison_df.iloc[:,0].abs())
            comparison_df.loc[:, 'ct_derived_signal'] = new_col
            comparison_df.loc[comparison_df.loc[:, 'ct_derived_signal'] == True, sw+'_real'] = 0
    
            term_genes_in_pseudo = comparison_df.loc[comparison_df.loc[:,'ct_derived_signal'] == True ,:].index.tolist()[:]
            dereg_genes = set(term_genes_in_real).difference(term_genes_in_pseudo)
            genes_dictionary.setdefault(term, set()).update(dereg_genes)
            
            sw_df = comparison_df.loc[:,[sw+'_real']]
            sw_df.columns = [sw]
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
    
    cumsum_df.to_csv(pseudo_results_dir+'final_scores/reactome_cumsum_values.tsv', sep='\t')
    absolute_df.to_csv(pseudo_results_dir+'final_scores/reactome_absolute_values.tsv', sep='\t')
    
    
    ###########################################################################
    # 2. Enrichment analysis between the remaining gene sets and cell type 
    # markers to identify disease-relevant cell types.
    ###########################################################################
    network_df = pandas.read_csv(main_dir+'MASLD_unified_undirected_network_with_semantics.tsv', sep='\t')
    network = igraph.Graph.DataFrame(network_df, directed=False, use_vids=False)
    background_genes = network.vs['name']
    
    cell_types = ['Circulating_NK_NKT_cells',
              'T_cells',
              'Cholangiocytes',
              'Hepatocytes',
              'Migratory_dendritic_cells',
              'Neutrophils',
              'Fibroblasts',
              'Macrophages',
              'Resident_NK_cells',
              'Endothelial_cells',
              'Basophils',
              'Plasmacytoid_dendritic_cells']
    
    cell_type_markers = {}
    tmp_dir = '../../data/sc_liver_cell_atlas/pscf_cell_type_markers_markers/'
    files = os.listdir(tmp_dir)
    files = list(filter(lambda x: x.split('.')[0] in cell_types, files))
    for f in files:
        cell_type = f.split('.')[0]
        df = pandas.read_csv(tmp_dir+f, sep='\t')
        genes = df.loc[:,'gene_symbol'].tolist()
        cell_type_markers.setdefault(cell_type, set()).update(set(genes))
        
    stats_df = []
    for term, dereg_genes in genes_dictionary.items():
        for cell_type, markers in cell_type_markers.items():
            M = len(background_genes)
            n = len(dereg_genes)
            N = len(markers)
            markers_in_term = set(markers).intersection(dereg_genes)
            x = len(markers_in_term)
            table = numpy.array([[x, n-x],[N-x, M-(n+N)+x]])
            pvalue = fisher_exact(table, alternative='greater').pvalue
            stats_df.append([term, cell_type, pvalue])

    stats_df = pandas.DataFrame(stats_df)
    adjusted_pvalues = multitest.multipletests(stats_df.iloc[:,2], method='fdr_bh')[1]
    stats_df.loc[:,'fdr'] = adjusted_pvalues
    stats_df.columns = ['term', 'cell_type', 'pvalue', 'fdr']        
    
    rownames = list(genes_dictionary.keys())
    colnames = cell_types        
    stats_matrix = numpy.ones((len(rownames),len(colnames)))
    stats_matrix = pandas.DataFrame(stats_matrix, columns=colnames, index=rownames)
    
    for index, row in stats_df.iterrows():
        stats_matrix.loc[row['term'], row['cell_type']] = row['fdr'] 
    stats_matrix = -numpy.log10(stats_matrix)
    stats_matrix_to_write = stats_matrix.copy()
    
    rownames = []
    for term in stats_matrix_to_write.index.tolist():
        definition = reactome_graph.vs(name=term)['definition'][0]
        rownames.append(definition)
    stats_matrix_to_write.index = rownames
    stats_matrix_to_write.to_csv(pseudo_results_dir+'final_scores/cell_type_enrichments.tsv', sep='\t')