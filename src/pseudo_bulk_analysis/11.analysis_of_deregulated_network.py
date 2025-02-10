#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas
import json
import igraph
import numpy
from scipy.stats import fisher_exact
from statsmodels.stats import multitest
import os


def get_pathways_activation_matrix(sw_signatures, pathway_signatures, up_weights, down_weights):
    
    sws = up_weights.columns.tolist()
    matrix = []
    
    for term, term_dict in pathway_signatures.items():
        #genes_set = set(term_dict['genes']).intersection(up_weights.index.tolist())
        genes_set = term_dict['genes']
        vector = [term]
        df = []
        for sw in sws:

            sw_genes = sw_signatures[sw]['up']
            not_in_sw_genes = list(set(genes_set).difference(sw_genes))
            in_sw_genes = list(set(genes_set).intersection(sw_genes))
            genes_set_up_weights = up_weights.loc[genes_set, sw].to_frame()
            genes_set_up_weights.loc[not_in_sw_genes, sw] = 0
            genes_set_up_weights.columns = ['up_'+sw]

            sw_genes = sw_signatures[sw]['down']
            not_in_sw_genes = list(set(genes_set).difference(sw_genes))
            in_sw_genes = list(set(genes_set).intersection(sw_genes))
            genes_set_down_weights = down_weights.loc[genes_set, sw].to_frame()
            genes_set_down_weights.loc[not_in_sw_genes, sw] = 0
            genes_set_down_weights = - genes_set_down_weights
            genes_set_down_weights.columns = ['down_'+sw]

            tmp_df = pandas.concat([genes_set_up_weights, genes_set_down_weights], axis=1)
            sw_weights = tmp_df.sum(axis=1)
            df.append(sw_weights)

        df = pandas.concat(df, axis=1)#.cumsum(axis=1)
        df = pandas.merge(pathway_networks[term], df, left_index=True, right_index=True)
        df = df.iloc[:,1:].multiply(df.iloc[:,0], axis="index")
        vector.extend(df.sum(axis=0).tolist())
        matrix.append(vector)

    return matrix



def tranform_matrix(matrix):
    rownames = []
    generic_terms = []
    for term in matrix.index.tolist():
        definition = pathway_signatures[term]['definition']
        generic_term = pathway_signatures[term]['generic_term']
        rownames.append(definition)
        if generic_term in list(colors_mapping.keys()) or generic_term.startswith('Metabolism'):
            generic_terms.append(generic_term)
        else:
            generic_terms.append('Others')
    output_matrix = matrix.copy()
    output_matrix.index = rownames
    
    sorting_df = pandas.DataFrame({'term':rownames, 'generic_term':generic_terms})
    sorting_df.sort_values(by='generic_term', inplace=True)
    sorted_terms = sorting_df.term.tolist()
    output_matrix = output_matrix.loc[sorted_terms,:].copy(deep=True)
    #matrix = matrix.cumsum(axis=1)

    signs = numpy.sign(output_matrix.values)
    output_log_matrix = numpy.log2(1+numpy.abs(output_matrix.values))
    output_log_matrix = signs*output_log_matrix
    output_log_matrix = pandas.DataFrame(output_log_matrix, columns=output_matrix.columns,
                                        index=output_matrix.index)
    return sorting_df, output_log_matrix



colors_mapping = {
    'Signal Transduction':'#38761d',
    'Immune System':'#cd9d4e',
    'Developmental Biology':'#79a686',
    'Disease':'#9c27b0',
    'Cellular responses to stimuli':'#00bcd4',
    'Extracellular matrix organization':'#be4b00',
    'Metabolism':'#660000', #red
    'Gene expression (Transcription)':'#3f51b5',
    'Cell Cycle':'#000000',
    'Muscle contraction':'#848999',
    'Transport of small molecules':'#673ab7',
    'Others':'#ff9800'
}




if __name__ == '__main__':
    
    main_dir = '../../results/'
    output_dir = '../../results/pseudo_bulk_analysis/networks/pathways_deconvolution/'
    if not os.path.exists(output_dir):
        os.makedirs(output_dir, exist_ok=True)
        
        
    filename = main_dir+'networks/network_analysis/consensus_propagation_signatures/random_sw-lapl_norm_weight.json'
    with open(filename, 'r') as fobj:
        real_signatures = json.load(fobj)
    
    filename = main_dir+'pseudo_bulk_analysis/networks/network_analysis/consensus_propagation_signatures/random_sw-lapl_norm_weight.json'
    with open(filename, 'r') as fobj:
        pseudo_signatures = json.load(fobj)
        
    
    
    
    ct_part = {}
    dereg_part = {}
    for sw, real_sw_dict in real_signatures.items():
        pseudo_sw_dict = pseudo_signatures[sw]
        
        real_up = set(real_sw_dict['up'])
        real_down = set(real_sw_dict['down'])
        
        pseudo_up = set(pseudo_sw_dict['up'])
        pseudo_down = set(pseudo_sw_dict['down'])
        
        ct_part.setdefault(sw, {}).update({'up': list(real_up.intersection(pseudo_up))})
        ct_part.setdefault(sw, {}).update({'down': list(real_down.intersection(pseudo_down))})
        
        dereg_part.setdefault(sw, {}).update({'up': list(real_up.difference(pseudo_up))})
        dereg_part.setdefault(sw, {}).update({'down': list(real_down.difference(pseudo_down))})
    
    filename = main_dir+'pseudo_bulk_analysis/networks/pathways_deconvolution/dereg_part.json'
    with open(filename, 'w') as fobj:
        json.dump(dereg_part, fobj)    
    
    filename = main_dir+'pseudo_bulk_analysis/networks/pathways_deconvolution/ct_part.json'
    with open(filename, 'w') as fobj:
        json.dump(ct_part, fobj)   
    
    
    # Reactome terms - signatures to plot
    filename =  '../../results/networks/network_analysis/networks_of_final_Reactome_pathways.json'
    with open(filename, 'r') as f:
        pathway_signatures = json.load(f)    
    
    # Network
    network_df = pandas.read_csv('../../results/networks/MASLD_unified_undirected_network.csv', sep=',')
    network = igraph.Graph.DataFrame(network_df, directed=False, use_vids=False)
    
    # Calcualte the centralities - weights of genes in each signature 
    pathway_networks = {}
    for term, term_dict in pathway_signatures.items():
        genes_set = term_dict['genes']
        term_network = network.induced_subgraph(genes_set, implementation='create_from_scratch')
        centralities = term_network.pagerank(damping=0.85, weights='weight')
        term_df = pandas.DataFrame(data={'centralities':centralities})
        term_df.index = term_network.vs['name']
        pathway_networks.update({term:term_df})    
        
        
    # Network node weights
    up_weights = pandas.read_csv(main_dir+'networks/network_analysis/initial_weights/up.tsv', sep='\t')
    down_weights = pandas.read_csv(main_dir+'networks/network_analysis/initial_weights/down.tsv', sep='\t')
    
    real_matrix = get_pathways_activation_matrix(real_signatures, pathway_signatures, 
                                                 up_weights, down_weights)    
    sws = up_weights.columns.tolist()
    real_matrix = pandas.DataFrame(real_matrix)
    real_matrix.columns = ['name'] + sws
    real_matrix.set_index('name', drop=True, inplace=True)
        
    
    dereg_matrix = get_pathways_activation_matrix(dereg_part, pathway_signatures, 
                                                  up_weights, down_weights)
    sws = up_weights.columns.tolist()
    dereg_matrix = pandas.DataFrame(dereg_matrix)
    dereg_matrix.columns = ['name'] + sws
    dereg_matrix.set_index('name', drop=True, inplace=True)    
    
    #ct_matrix = get_pathways_activation_matrix(ct_part, pathway_signatures, 
    #                                           up_weights, down_weights)
    #sws = up_weights.columns.tolist()
    #ct_matrix = pandas.DataFrame(ct_matrix)
    #ct_matrix.columns = ['name'] + sws
    #ct_matrix.set_index('name', drop=True, inplace=True)    
    #sorting_df, ct_log_matrix = tranform_matrix(ct_matrix)

    sorting_df, dereg_log_matrix = tranform_matrix(dereg_matrix)

    network_nodes = network.vs['name']
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
        
        
    rownames = list(pathway_signatures.keys())
    colnames = cell_types        
    stats_matrix = numpy.ones((len(rownames),len(colnames)))
    stats_matrix = pandas.DataFrame(stats_matrix, columns=colnames, index=rownames)
    
    ct_part_genes = set()
    for sw, tmp_dict in ct_part.items():
        ct_part_genes.update(tmp_dict['up'])
        ct_part_genes.update(tmp_dict['down'])
    
    dereg_part_genes = set()
    for sw, tmp_dict in dereg_part.items():
        dereg_part_genes.update(tmp_dict['up'])
        dereg_part_genes.update(tmp_dict['down'])
        
    background_genes = network_nodes[:]
    stats_df = []
    
    for pathway in rownames:    
        pathway_genes = set(pathway_signatures[pathway]['genes']).difference(ct_part_genes)
        for cell_type in cell_types:
            M = len(background_genes)
            n = len(pathway_genes)
            markers = cell_type_markers[cell_type]
            N = len(markers)
            markers_in_pathway = set(markers).intersection(pathway_genes)
            x = len(markers_in_pathway)
            table = numpy.array([[x, n-x],[N-x, M-(n+N)+x]])
            pvalue = fisher_exact(table, alternative='greater').pvalue
            stats_df.append([pathway, cell_type, pvalue])    
                
    stats_df = pandas.DataFrame(stats_df)
    adjusted_pvalues = multitest.multipletests(stats_df.iloc[:,2], method='fdr_bh')[1]
    stats_df.loc[:,'fdr'] = adjusted_pvalues
    stats_df.columns = ['pathway', 'cell_type', 'pvalue', 'fdr']        
    for index, row in stats_df.iterrows():
        stats_matrix.loc[row['pathway'], row['cell_type']] = row['fdr']        
        
    log_matrix = -numpy.log10(stats_matrix)
    rownames = []
    generic_terms = []
    for term in log_matrix.index.tolist():
        definition = pathway_signatures[term]['definition']
        generic_term = pathway_signatures[term]['generic_term']
        rownames.append(definition)
        if generic_term in list(colors_mapping.keys()) or generic_term.startswith('Metabolism'):
            generic_terms.append(generic_term)
        else:
            generic_terms.append('Others')
    log_matrix.index = rownames
            
    sorting_df = pandas.DataFrame({'term':rownames, 'generic_term':generic_terms})
    sorting_df.sort_values(by='generic_term', inplace=True)
    sorted_terms = sorting_df.term.tolist()
    log_matrix = log_matrix.loc[sorted_terms,:].copy(deep=True)        
        
        
    log_matrix.rename({'Circulating_NK_NKT_cells':'Circulating_NK/NKT_cells'}, axis=1, inplace=True)        
    log_matrix.columns = list(map(lambda x: x.replace('_', ' '), log_matrix.columns.tolist()))
    #log_matrix = log_matrix.loc[(log_matrix > 1).any(axis=1), :]
    #sorting_df = sorting_df.loc[sorting_df.term.isin(log_matrix.index),:]
    filename = main_dir+'pseudo_bulk_analysis/networks/network_analysis/cell_types_enrichment_in_dereg_part.tsv'
    log_matrix.to_csv(filename, sep='\t')
    

        
        
        
        
        
        
        