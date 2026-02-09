#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas # 2.1.4
import os
import igraph # 0.11.4
from reactome_library import OntoGraph


###############################################################################
# Description:
###############################################################################
# This scirpt uses the slim version of Reactome and the results of enrichment
# analysis on MASLD nodes and network signatures in order to further filter the 
# Reactome graph and generate a MASLD-related version of it, keeping the nodes 
# which have been found strongly enriched to the derived network signatures.
# Outputs:
# - A pkl object for the enriched Reactome graph.
###############################################################################


if __name__ == '__main__':
    
    ###########################################################################
    # Inputs
    ###########################################################################
    main_dir = '../../results/ucam_sanyal/networks/'
    output_dir = main_dir+'network_analysis/'
    
    filename = main_dir+'MASLD_unified_undirected_network_with_semantics.tsv'
    network_df = pandas.read_csv(filename, sep='\t')
    masld_network = igraph.Graph.DataFrame(network_df, directed=False, use_vids=False)
    
    reference_reactome_graph = igraph.Graph.Read_Pickle('../../data/reactome_graph/reactome_graph_pruned.pkl')
    reactome_graph = igraph.Graph.Read_Pickle('../../data/reactome_graph/reactome_graph_with_semantics.pkl')
    
    ###########################################################################
    # 1. Creation of the MASLD network-related Reactome graph. Enrichment
    # analysis of network nodes revealed a group of significantly enriched terms.
    # These terms are used to create a subset of the pruned Reactome graph as
    # the MASLD-related one.
    ###########################################################################
    root_node = 'R-HSA-Reactome Pathway'
    # Retrieve the enriched terms for the whole network
    tmp_df = pandas.read_csv(main_dir+'MASLD_nodes_Reactome.tsv', sep='\t', index_col=0)
    tmp = list(map(lambda x: x.split(' ')[-1], tmp_df.index.tolist()))
    tmp_df.loc[:,'term_id'] = tmp
    tmp_df = tmp_df.loc[~pandas.isnull(tmp_df.genes),:]
    tmp_df = tmp_df.loc[tmp_df.loc[:, 'adj_p.value'] < 0.05,:]
    masld_network_terms = tmp_df.term_id.tolist() # len(network_terms) = 434
    
    # Retrieve all the paths to the Reactome root of the enriched terms, to construct the MASLD-related Reactome graph
    paths = reactome_graph.get_all_simple_paths(v=root_node, to=masld_network_terms)
    terms = [t for p in paths for t in p]
    masld_network_subgraph = reactome_graph.induced_subgraph(vertices=set(terms))
    
    # Subset the MASLD-related graph to keep only those terms which are in the pruned Reactome graph
    intersection = list(set(masld_network_subgraph.vs['name']).intersection(reference_reactome_graph.vs['name']))
    masld_network_subgraph = masld_network_subgraph.induced_subgraph(vertices=intersection, 
                                                                     implementation='create_from_scratch')
    
    ###########################################################################
    # 2. Further filtering, keeping only those terms which have been found
    # as enriched in the signatures from network analysis.
    ###########################################################################
    tmp_dir = main_dir+'network_analysis/Reactome_for_consensus_propagation_signatures/'
    files = os.listdir(tmp_dir)
    enriched_terms_dict = {} # mapping of enriched terms and genes
    terms_dict = {} # general mapping of terms and genes (even not the enriched ones)
    for filename in files:
        tmp_df = pandas.read_csv(tmp_dir+filename, sep='\t', index_col=0)
        tmp = list(map(lambda x: x.split(' ')[-1], tmp_df.index.tolist()))
        tmp_df.loc[:,'term_id'] = tmp
        tmp_df = tmp_df.loc[~pandas.isnull(tmp_df.genes),:]
        # get all the anotation
        tmp = zip(tmp_df.term_id.tolist(), tmp_df.genes.tolist())    
        for entry in tmp:
            terms_dict.setdefault(entry[0], set()).update(entry[1].split(';'))    
        # get only the annotation of enriched terms
        tmp_df = tmp_df.loc[tmp_df.loc[:, 'adj_p.value'] < 0.05,:]
        tmp = zip(tmp_df.term_id.tolist(), tmp_df.genes.tolist())
        for entry in tmp:
            if entry[0] in masld_network_subgraph.vs['name']:
                enriched_terms_dict.setdefault(entry[0], set()).update(entry[1].split(';'))

    # Filter the enriched terms and keep only those that are included in the MASLD-related graph 
    enriched_terms_in_reference = list(set(masld_network_subgraph.vs['name']).intersection(enriched_terms_dict.keys()))
    terms_for_enriched_graph = enriched_terms_in_reference[:]
    for term in enriched_terms_in_reference:
        ancs = [a['name'] for a in reference_reactome_graph.get_vertex_ancestors(term)]
        terms_for_enriched_graph.extend(ancs)
    terms_for_enriched_graph = list(set(terms_for_enriched_graph))
    enriched_graph = masld_network_subgraph.induced_subgraph(terms_for_enriched_graph,
                                                             implementation='create_from_scratch')
    
    ###########################################################################
    # 3. Two last filterings:
    # - Leaves which are generic terms (root's children)
    # - Terms with less than 10 genes in the de-activated signatures
    ###########################################################################
    generic_terms = enriched_graph.get_vertex_children(root_node)
    generic_terms = [t['name'] for t in generic_terms]
    leaves = []
    for term in enriched_graph.vs:
        out_degree = enriched_graph.degree(term['name'], mode='out')
        if out_degree == 0:
            leaves.append(term['name'])
    to_remove = list(set(leaves).intersection(generic_terms))
  
    for node in enriched_graph.vs:
        term = node['name']
        reference_annotation = reference_reactome_graph.vs.find(name=term)['annotation']
        network_annotation = list(set(masld_network.vs['name']).intersection(reference_annotation))
        term_network = masld_network.induced_subgraph(network_annotation, implementation='create_from_scratch')
        centralities = term_network.pagerank(damping=0.85, weights='weight')
        term_df = pandas.DataFrame(data={'centralities':centralities})
        term_df.index = term_network.vs['name']
        node['annotation'] = network_annotation
        node['centralities'] = term_df
        try:
            genes_in_signatures = enriched_terms_dict[term]
            if len(genes_in_signatures) < 10:
                to_remove.append(term)
        except KeyError:
            pass
        
    leaves = list(set(leaves).difference(to_remove))
    enriched_graph.delete_vertices(to_remove)
    
    # Annotate the terms that will be used in the heatmap
    leaves = []
    for term in enriched_graph.vs:
        out_degree = enriched_graph.degree(term['name'], mode='out')
        if out_degree == 0:
            leaves.append(term['name'])
    leaf_parents = []
    for leaf in leaves:
        parents = enriched_graph.get_vertex_parents(leaf)
        for parent in parents:
            leaf_parents.extend([parent['name']])
    leaf_parents = list(set(leaf_parents))
    examined_terms = leaves[:]
    for parent in leaf_parents:
        parent_node = enriched_graph.vs.find(name=parent)
        children = enriched_graph.get_vertex_children(parent)
        children_annotation = [c['annotation'] for c in children]
        children_annotation = set([g for l in children_annotation for g in l])
        parent_annotation = parent_node['annotation']
        covered_ratio = round(len(children_annotation)/len(parent_annotation),2)
        if covered_ratio < 0.8 and not parent in generic_terms:
            examined_terms.append(parent)
    
    for term in examined_terms:
        node = enriched_graph.vs.find(name=term)
        node['plot'] = True
        ancestors = [a['name'] for a in enriched_graph.get_vertex_ancestors(term)]
        ancestors = list(set(ancestors).intersection(generic_terms))
        node['generic_term'] = ancestors
    
    filename = output_dir+'reactome_enriched_graph.pkl'
    enriched_graph.write_pickle(filename)
    
    
    