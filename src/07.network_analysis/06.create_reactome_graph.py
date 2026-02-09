#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas # 2.1.4
import os
import numpy # 1.26.4
import json
from reactome_library import OntoGraph


###############################################################################
# Description:
###############################################################################
# This scirpt creates a directed acyclic graph (DAG) for Reactome database and
# then filters it in order to generate a slim version, keeping the most generic
# but still informative nodes.
# Outputs:
# - pkl objects for the whole and pruned Reactome graphs.
###############################################################################

    
def create_reactome_graph():
    
    filename = '../../data/annotation_databases/raw_data/reactome/ReactomePathways.txt'
    reactome_base = pandas.read_csv(filename, sep='\t', header=None)
    reactome_base = reactome_base.loc[reactome_base.iloc[:,2] == 'Homo sapiens',]
    reactome_base.reset_index(drop=True, inplace=True)
    reactome_base.columns = ['term_id', 'definition', 'species']
    reactome_base.loc[:,'graph_id'] = reactome_base.index.tolist()
    
    filename = '../../data/annotation_databases/raw_data/reactome/ReactomePathwaysRelation.txt'
    reactome_relations = pandas.read_csv(filename, sep='\t', header=None)
    bool1 = reactome_relations.iloc[:,0].isin(reactome_base.iloc[:,0].tolist())
    bool2 = reactome_relations.iloc[:,1].isin(reactome_base.iloc[:,0].tolist())
    reactome_relations = reactome_relations.loc[bool1 | bool2,:]
    reactome_relations.reset_index(drop=True, inplace=True)
    reactome_relations.columns = ['parent', 'child']
    
    source = reactome_relations.apply(lambda x: reactome_base.loc[reactome_base.term_id == x.parent,'graph_id'].values[0], 
                                      axis=1)
    target = reactome_relations.apply(lambda x: reactome_base.loc[reactome_base.term_id == x.child,'graph_id'].values[0], 
                                      axis=1)
    reactome_relations.loc[:,'source'] = source
    reactome_relations.loc[:,'target'] = target
    reactome_relations = reactome_relations.loc[:, ['source', 'target']]
    reactome_base.sort_values(by='graph_id', inplace=True)
    
    reactome_graph = OntoGraph().DataFrame(reactome_relations)
    reactome_graph.vs['name'] = reactome_base.term_id.tolist()
    reactome_graph.vs['term_id'] = reactome_base.term_id.tolist()
    reactome_graph.vs['definition'] = reactome_base.definition.tolist()
    
    root_node_name = 'R-HSA-Reactome Pathway'
    root_node = reactome_graph.add_vertex(root_node_name)
    root_node = reactome_graph.vs.find(name=root_node_name)
    reactome_graph.vs[root_node.index]['definition'] = root_node_name
    reactome_graph.vs[root_node.index]['term_id'] = root_node_name
        
    degrees = numpy.array(reactome_graph.degree(mode='in'))
    first_level_terms = list(numpy.where(degrees == 0)[0])
    extra_edges = []
    for term_index in first_level_terms:
        if term_index == root_node.index:
            pass
        else:
            extra_edges.append((root_node.index, term_index))
    reactome_graph.add_edges(extra_edges)
    
    return reactome_graph



if __name__ == '__main__':

    ###########################################################################
    # 1. Create the graph and load the annotation
    ###########################################################################    
    reactome_graph = create_reactome_graph()
    root_node = reactome_graph.get_root()

    with open('../../data/annotation_databases/raw_data/reactome/reactome_2024.json') as fobj:
        reactome_annotation = json.load(fobj)        
    all_genes = set()
    for term, genes in reactome_annotation.items():
        reactome_graph.vs.find(name=term)['annotation'] = genes
        all_genes.update(genes)
    reactome_graph.vs.find(name='R-HSA-Reactome Pathway')['annotation'] = list(all_genes)
    
    output_dir = '../../data/reactome_graph/'
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    reactome_graph.write_pickle(output_dir+'reactome_graph_global.pkl')

    ###########################################################################
    # 2. Remove nodes without annotation
    ###########################################################################
    indices_to_remove = numpy.where(pandas.isnull(reactome_graph.get_vertex_dataframe().annotation))
    indices_to_remove = list(indices_to_remove[0])
    reactome_graph.delete_vertices(indices_to_remove)
    
    ###########################################################################
    # 3. Calculate the semantic attributes for graph nodes (IC and SV)
    ###########################################################################
    p = 1/len(all_genes)
    max_ic_annot = -numpy.log2(p)
    p = 1/len(reactome_graph.vs)
    max_ic_graph = -numpy.log2(p)
    
    for node in reactome_graph.vs:
        p = len(node['annotation'])/len(all_genes)
        ic = -numpy.log2(p)
        node['IC_annot'] = ic/max_ic_annot
        p = (1+len(reactome_graph.get_vertex_descendants(node['name'])))/len(reactome_graph.vs)
        ic = -numpy.log2(p)
        node['IC_graph'] = ic/max_ic_graph
        if ic > 0:
            node['SW'] = 1/(1+numpy.exp(-max_ic_graph/ic))
        else:
            node['SW'] = 1
    
    for node in reactome_graph.vs:
        ancs = reactome_graph.get_vertex_ancestors(node['name'])
        node['SV'] = sum([anc['SW'] for anc in ancs]) + node['SW']

    for node in reactome_graph.vs:
        children = reactome_graph.get_vertex_children(node['name'])
        if len(children) != 0:
            diffs = [c['SV']-node['SV'] for c in children]
            diff = numpy.mean(diffs)
            node['SV_diff'] = diff
            diffs = [c['IC_annot']-node['IC_annot'] for c in children]
            diff = numpy.mean(diffs)
            node['IC_annot_diff'] = diff
        else:
            node['SV_diff'] = numpy.nan
            node['IC_annot_diff'] = numpy.nan

    reactome_graph.write_pickle(output_dir+'reactome_graph_with_semantics.pkl')

    ###########################################################################
    # 4. First pruning of graph, keeping the most generic terms.
    # The most generic terms are defined using the IC and SV values. All terms 
    # with values higher than the respective 25th percentile (or condition) 
    # are filtered out. These terms can be considered as extremely specific.
    ###########################################################################
    # Prune the graph
    filtered_reactome_graph = reactome_graph.copy()
    
    # Filtering no1: IC and SV thresholds
    ic_threshold = numpy.percentile(filtered_reactome_graph.vs['IC_annot'], 25)
    sv_threshold = numpy.percentile(filtered_reactome_graph.vs['SV'], 25)
    
    to_remove = [t['name'] for t in filtered_reactome_graph.vs if t['IC_annot']>=ic_threshold or t['SV']>=sv_threshold]
    filtered_reactome_graph.delete_vertices(to_remove)

    ###########################################################################
    # 5. Second pruning of graph, by examining the semantic distances between
    # parent and childs. If these distances are lower that the respective
    # median values (and condition), then the descendants terms are filtered out
    # because they are semantically close to their parent.
    ###########################################################################
    leaves = []
    for term in filtered_reactome_graph.vs:
        out_degree = filtered_reactome_graph.degree(term['name'], mode='out')
        if out_degree == 0:
            leaves.append(term['name'])
    leaf_parents = []
    for leaf in leaves:
        parents = filtered_reactome_graph.get_vertex_parents(leaf)
        for parent in parents:
            leaf_parents.extend([parent['name']])
    leaf_parents = list(set(leaf_parents))
    
    corrected_leaf_parents = []
    for term in leaf_parents:
        children = [c['name'] for c in filtered_reactome_graph.get_vertex_children(term)]
        if len(set(leaf_parents).intersection(children)) > 0:
            pass
        else:
            corrected_leaf_parents.append(term)    

    sv_diff_thr = numpy.nanpercentile(filtered_reactome_graph.vs['SV_diff'], 50)
    ic_annot_diff_thr = numpy.nanpercentile(filtered_reactome_graph.vs['IC_annot_diff'], 50)
    generic_terms = filtered_reactome_graph.get_vertex_children('R-HSA-Reactome Pathway')
    generic_terms = [t['name'] for t in generic_terms]
    print(ic_annot_diff_thr, sv_diff_thr)
    to_remove = []
    for term in corrected_leaf_parents:
        term_sv_diff = filtered_reactome_graph.vs.find(name=term)['SV_diff']
        term_ic_annot_diff = filtered_reactome_graph.vs.find(name=term)['IC_annot_diff']
        if (term_sv_diff >= sv_diff_thr and term_ic_annot_diff >= ic_annot_diff_thr) or term in generic_terms:
            pass
        else:
            to_remove.extend([c['name'] for c in filtered_reactome_graph.get_vertex_descendants(term)])
    to_remove = list(set(to_remove))    
    filtered_reactome_graph.delete_vertices(to_remove)
    filtered_reactome_graph.write_pickle(output_dir+'reactome_graph_pruned.pkl')
