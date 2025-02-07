#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas
import os
import igraph
import numpy
from scipy import spatial
import json


class OntoGraph(igraph.Graph):

    def __init__(self,  *args, **kwargs):
        super().__init__(*args, **kwargs)

    def get_vertex_parents(self, node):
        if type(node) == int:
            in_edges = self.incident(node, mode='in')
        else:
            node = self.vs.find(name=node)
            in_edges = self.incident(node.index, mode='in')
        parents = []
        for e in in_edges:
            parents.append(self.vs[self.es[e].source])
        return parents
    
    def get_vertex_ancectors(self, node):
        if type(node) == int:
            ancestors = self.bfs(node.index, mode='in')[0][1:]
        else:
            node = self.vs.find(name=node)
            ancestors = self.bfs(node.index, mode='in')[0][1:]
        return self.vs[ancestors]

    def get_vertex_children(self, node):
        if type(node) == int:
            out_edges = self.incident(node, mode='out')
        else:
            node = self.vs.find(name=node)
            out_edges = self.incident(node.index, mode='out')
        children = []
        for e in out_edges:
            children.append(self.vs[self.es[e].target])
        return children
    
    def get_vertex_descendants(self, node):
        if type(node) == int:
            descendants = self.bfs(node, mode='out')[0][1:]
        else:
            node = self.vs.find(name=node)
            descendants = self.bfs(node.index, mode='out')[0][1:]
        return self.vs[descendants]
    
    def get_distance_from_root(self, node):
        degrees = zip(self.vs['name'], self.degree(mode='in'))
        root_node = list(filter(lambda x: x[1]==0, degrees))[0][0]
        distance = self.distances(source=root_node, target=node)[0][0]
        return distance

    def get_root(self):
        degrees = zip(self.vs['name'], self.degree(mode='in'))
        root_node = list(filter(lambda x: x[1]==0, degrees))[0][0]
        return root_node
    


def create_reactome_graph():
    
    filename = '../../data/annotation_databases/reactome/ReactomePathways.txt'
    reactome_base = pandas.read_csv(filename, sep='\t', header=None)
    reactome_base = reactome_base.loc[reactome_base.iloc[:,2] == 'Homo sapiens',]
    reactome_base.reset_index(drop=True, inplace=True)
    reactome_base.columns = ['term_id', 'definition', 'species']
    reactome_base.loc[:,'graph_id'] = reactome_base.index.tolist()
    
    filename = '../../data/annotation_databases/reactome/ReactomePathwaysRelation.txt'
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
    
    main_dir = '../../results/networks/network_analysis/'
    output_dir = main_dir+'Reactome_pathways_clustering/'
    os.makedirs(output_dir, exist_ok=True)
    
    reactome_graph = create_reactome_graph()
    root_node = reactome_graph.get_root()
    
    ###########################################################################
    #
    # Creation of the reactome graph which contains 
    #
    ###########################################################################
    
    # Retrieve enrichment analysis results for the whole network
    tmp_df = pandas.read_csv('../../results/networks/MASLD_nodes_Reactome.tsv',
                             sep='\t', index_col=0)
    tmp = list(map(lambda x: x.split(' ')[-1], tmp_df.index.tolist()))
    tmp_df.loc[:,'term_id'] = tmp
    tmp_df = tmp_df.loc[~pandas.isnull(tmp_df.genes),:]
    tmp_df = tmp_df.loc[tmp_df.loc[:, 'adj_p.value'] < 0.01,:]
    network_terms = tmp_df.term_id.tolist() # len(network_terms) = 434
    
    # Get the paths that connect network nodes with root_node and induce the
    # correspoding part of the reactome graph
    paths = reactome_graph.get_all_simple_paths(v=root_node, to=network_terms)
    terms = [t for p in paths for t in p]
    network_subgraph = reactome_graph.induced_subgraph(vertices=set(terms))
    # len(network_subgraph.vs['name']) = 502
    
    # Find the generic terms in the graph and split them to enriched and
    # non-enriched. Descendants of non-enriched terms will be removed from the
    # graph.
    generic_terms = network_subgraph.get_vertex_children(root_node)
    generic_terms = [n['name'] for n in generic_terms]
    enriched_generic_terms = tmp_df.loc[tmp_df.term_id.isin(generic_terms),'term_id'].tolist()
    
    non_enriched_generic_terms = set(generic_terms).difference(enriched_generic_terms)
    network_subgraph.delete_vertices(non_enriched_generic_terms)
    main_component = network_subgraph.components(mode='weak')[0]
    network_subgraph = network_subgraph.induced_subgraph(vertices=main_component, implementation='create_from_scratch')
    # len(network_subgraph.vs['name']) 465
    
    
    ###########################################################################
    #
    # Create the enriched ractome graph, using the enrichment analysis results 
    # for each sw
    #
    ###########################################################################
    
    tmp_dir = main_dir+'Reactome_for_consensus_propagation_signatures/'
    files = os.listdir(tmp_dir)
    
    enriched_terms_data = {} # mapping of terms and genes (terms enriched in the whole network)
    terms_data = {} # general mapping of terms and genes (even not the enriched ones)
    for f in files:
        tmp_df = pandas.read_csv(tmp_dir+f, sep='\t', index_col=0)
        tmp = list(map(lambda x: x.split(' ')[-1], tmp_df.index.tolist()))
        tmp_df.loc[:,'term_id'] = tmp
        tmp_df = tmp_df.loc[~pandas.isnull(tmp_df.genes),:]
        # get all the anotation
        tmp = zip(tmp_df.term_id.tolist(), tmp_df.genes.tolist())    
        for entry in tmp:
            terms_data.setdefault(entry[0], set()).update(entry[1].split(';'))    
        # get only the annotation of enriched terms
        tmp_df = tmp_df.loc[tmp_df.loc[:, 'adj_p.value'] < 0.05,:]
        tmp = zip(tmp_df.term_id.tolist(), tmp_df.genes.tolist())
        for entry in tmp:
            if entry[0] in network_terms:
                enriched_terms_data.setdefault(entry[0], set()).update(entry[1].split(';'))

    # Get only the eriched terms, connect them up to the root and find the 
    # respective subgraph of network_subgraph
    terms_in_subgraph = set(list(enriched_terms_data.keys())).intersection(network_subgraph.vs['name'])
    paths = network_subgraph.get_all_simple_paths(v=root_node, to=list(terms_in_subgraph))
    all_terms = set([t for p in paths for t in p])
    enriched_subgraph = network_subgraph.induced_subgraph(vertices=all_terms)
    
    # Annotate the terms with genes using the terms_data
    for term in enriched_subgraph.vs['name']:
        try:
            genes = list(set(terms_data[term]))
        except KeyError:
            genes = []
        terms_data.update({term:genes})
        enriched_subgraph.vs.find(name=term)['genes'] = genes

    # Genes in the enriched graph - not all the deregulated genes are annotated
    # with enriched terms
    genes_in_the_graph = set([g for l in enriched_subgraph.vs['genes'] for g in l])
    F = open(main_dir+'deregulated_network_genes.txt', 'w')
    F.write('\n'.join(genes_in_the_graph))
    F.close()


    tmp_json = {}
    for term_dict in enriched_subgraph.vs:
        d = {'genes':term_dict['genes'], 'definition':term_dict['definition']}
        tmp_json.update({term_dict['name']: d})
        
    with open(main_dir+'deregulated_pathways.json', 'w') as fobj:
        json.dump(tmp_json, fobj)
    
    
    ###########################################################################
    #
    # Pruning step 1: remove all terms with distance > 4 from the root
    #
    ###########################################################################
    to_write = set()
    to_remove = set()
    for term in enriched_subgraph.vs['name']:
        distance = enriched_subgraph.get_distance_from_root(term)
        if distance > 4:
            to_remove.update([term])    
    enriched_subgraph.delete_vertices(to_remove)
    to_write.update(to_remove)
    
    degrees = list(zip(enriched_subgraph.vs['name'], enriched_subgraph.degree()))
    to_remove = set()
    for term, degree in degrees:
        if degree == 0:
            to_remove.update([term])    
    enriched_subgraph.delete_vertices(to_remove)
    to_write.update(to_remove)

    F = open(output_dir+'terms_with_dist4.txt', 'w')
    F.write('\n'.join(list(to_write)))
    F.close()

    ###########################################################################
    #
    # Pruning step 2: examine the branches of terms with distance=3 from the 
    # root
    #
    ###########################################################################
    selected_terms = []
    level = 3
    for term in enriched_subgraph.vs['name']:
        distance = enriched_subgraph.get_distance_from_root(term)
        if distance == 3:
            selected_terms.append(term)
    
    to_remove = set()
    to_keep = set()
    to_write = []
    new_edges = list()
    for level_term in selected_terms:
        children = [t['name'] for t in enriched_subgraph.get_vertex_children(level_term)]
        if len(children) > 1:
            to_remove.update(children)
            for c in children:
                to_write.append([level_term, c])
        elif len(children) == 1:
            c_genes = enriched_subgraph.vs.find(name=children[0])['genes']
            p_genes = enriched_subgraph.vs.find(name=level_term)['genes']
            jac_sim = len(c_genes)/len(p_genes)
            if jac_sim >= 0.8:
                ancs = [t['name'] for t in enriched_subgraph.get_vertex_parents(level_term)]
                for anc in ancs:
                    new_edges.append([anc, children[0]])
                to_remove.update([level_term])
                to_keep.update([children[0]])
                to_write.append([children[0],level_term])
            else:
                to_remove.update([children[0]])
                to_write.append([level_term, children[0]])
        else:
            pass    
    to_remove = to_remove.difference(to_keep)
    
    enriched_subgraph.delete_vertices(list(to_remove))
    for edge in new_edges:
        enriched_subgraph.add_edge(source=edge[0], target=edge[1])
    
    degrees = list(zip(enriched_subgraph.vs['name'], enriched_subgraph.degree()))
    to_remove = set()
    for term, degree in degrees:
        if degree == 0:
            to_remove.update([term])
    enriched_subgraph.delete_vertices(list(to_remove))





    ###########################################################################
    #
    # Pruning step 3: examine the branches of terms with distance=2 from the 
    # root
    #
    ###########################################################################
    selected_terms = []
    level = 2
    for term in enriched_subgraph.vs['name']:
        distance = enriched_subgraph.get_distance_from_root(term)
        if distance == 2:
            selected_terms.append(term)
            
    to_remove = set()
    to_keep = set()
    new_edges = list()
    for level_term in selected_terms:
        children = [c['name'] for c in enriched_subgraph.get_vertex_children(level_term)]
        if len(children) == 1:
            c_genes = enriched_subgraph.vs.find(name=children[0])['genes']
            p_genes = enriched_subgraph.vs.find(name=level_term)['genes']
            jac_sim = len(c_genes)/len(p_genes)
            if jac_sim >= 0.8:
                ancs = [t['name'] for t in enriched_subgraph.get_vertex_parents(level_term)]
                for anc in ancs:
                    new_edges.append([anc, children[0]])
                to_remove.update([level_term])
                to_keep.update([children[0]])
                to_write.append([children[0],level_term])
            else:
                to_remove.update([children[0]])
                to_write.append([level_term,children[0]])
        elif len(children) > 1:
            genes_mapping = dict((c, enriched_subgraph.vs.find(name=c)['genes']) for c in children)
            genes_pool = set([g for t,l in genes_mapping.items() for g in l])
            matrix = numpy.zeros((len(genes_mapping.keys()), len(genes_pool)))
            matrix = pandas.DataFrame(matrix)
            matrix.index = genes_mapping.keys()
            matrix.columns = genes_pool
            for t, l in genes_mapping.items():
                for g in l:
                    matrix.loc[t,g] = 1
            tmp_distances = spatial.distance.pdist(matrix.values, metric='jaccard')
            #sim_matrix = spatial.distance.squareform(1-tmp_distances)
            tmp_sims = 1 - tmp_distances
            if numpy.median(tmp_sims) >= 0.3:
                to_remove.update(set(genes_mapping.keys()))
                for c in set(genes_mapping.keys()):
                    to_write.append([level_term, c])
        else:
            pass    
    
    enriched_subgraph.delete_vertices(list(to_remove))
    for edge in new_edges:
        enriched_subgraph.add_edge(source=edge[0], target=edge[1])    


    tmp_df = pandas.DataFrame(to_write, columns=['keep', 'delete'])
    tmp_df.to_csv(output_dir+'term_substitutions.csv')
    ###########################################################################
    #
    # Define the final terms - leaves of the remaining graph
    #
    ###########################################################################
    generic_terms = enriched_subgraph.get_vertex_children(root_node)
    generic_terms = [t['name'] for t in generic_terms]
    data = []
    for generic_term in generic_terms:
        for term_dict in enriched_subgraph.get_vertex_descendants(generic_term):
            d = enriched_subgraph.get_distance_from_root(term_dict['name'])
            out_degree = enriched_subgraph.degree(term_dict['name'], mode='out')
            data.append([generic_term, term_dict['name'], d, out_degree==0])
    data = pandas.DataFrame(data)
    data.columns = ['generic_term_id', 'term_id', 'distance', 'is_leaf']
    data.to_csv(output_dir+'terms_after_clustering.csv')
    tmp_df = data.loc[data.loc[:, 'is_leaf'],:]
    final_terms = tmp_df.term_id.tolist()
    

    counts = []
    for term in final_terms:
        counts.append(len(enriched_subgraph.vs.find(name=term)['genes']))
    
    median = numpy.median(counts)
    deviations = numpy.abs(numpy.median(counts) - numpy.array(counts))
    thr = median - numpy.median(deviations)
    
    to_remove = set()
    for term in final_terms:
        if len(enriched_subgraph.vs.find(name=term)['genes']) <= thr:
            to_remove.update([term])

    final_terms = list(set(final_terms).difference(to_remove))

    F = open(output_dir+'leaf_terms_with_low_genes_count.txt', 'w')
    F.write('\n'.join(list(to_remove)))
    F.close()


    ###########################################################################
    #
    # Correct the final terms by adding some generic ones, due to the absent
    # of genes from the already defined final terms
    #
    ###########################################################################
    final_terms_signatures = {}
    generic_terms_counts = set()
    for term in final_terms:
        term_genes = enriched_subgraph.vs.find(name=term)['genes']
        definition = enriched_subgraph.vs.find(name=term)['definition']
        ancestors = [t['name'] for t in enriched_subgraph.get_vertex_ancectors(term)]
        ancestor = list(set(ancestors).intersection(generic_terms))[0]
        anc_definition = enriched_subgraph.vs.find(name=ancestor)['definition']
        final_terms_signatures.update({term:{'genes':term_genes, 
                                             'definition':definition,
                                             'generic_term':anc_definition}})
        generic_terms_counts.update([anc_definition])    
    
    for generic_term in generic_terms:
        generic_term_genes = enriched_subgraph.vs.find(name=generic_term)['genes']
        definition = enriched_subgraph.vs.find(name=generic_term)['definition']
        descendants = [t['name'] for t in enriched_subgraph.get_vertex_descendants(generic_term)]
        # find the leaves
        leaves = []
        for desc in descendants:
            out_degree = enriched_subgraph.degree(desc, mode='out')
            if out_degree == 0:
                leaves.append(desc)
        # get leaves' genes
        leaves_genes = set()
        for leaf in leaves:
            leaves_genes.update(enriched_subgraph.vs.find(name=leaf)['genes'])
        # check the similarity with the generic term gene set
        jac_sim =  len(leaves_genes)/len(generic_term_genes)
        if jac_sim >= 0.8:
            pass
        else:
            new_term = 'Subset of '+generic_term
            if definition in generic_terms_counts:
                new_definition = 'Other genes related to '+definition.lower()
                new_definition = new_definition.replace('dna', 'DNA').replace('rna', 'RNA')
            else:
                new_definition = 'Genes related to '+definition.lower()
            new_gene_set = list(set(generic_term_genes).difference(leaves_genes))
            final_terms_signatures.update({new_term:{'genes':new_gene_set,
                                                     'definition':new_definition,
                                                     'generic_term':definition}})
            
    print(len(final_terms_signatures))
    filename = main_dir+'networks_of_final_Reactome_pathways.json'
    with open(filename,'w') as f:
        json.dump(final_terms_signatures, f)
        
        
        
        
        
        
        
        
        
        