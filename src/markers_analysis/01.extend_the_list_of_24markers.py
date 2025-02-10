#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas
import numpy
import os
import igraph
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
    output_dir = main_dir+'24_markers/'
    if not os.path.exists(output_dir):
        os.makedirs(output_dir, exist_ok=True)

    F = open('../../data/24_plasma_markers.txt', 'r')
    genes24 = [line.replace('\n', '') for line in F.readlines()]
    F.close()
    
    # Results from network propagation
    filename =  main_dir+'consensus_propagation_signatures/random_sw-lapl_norm_weight.json'
    with open(filename, 'r') as f:
        sw_signatures = json.load(f)
    
    # Network node weights
    up_weights = pandas.read_csv(main_dir+'initial_weights/up.tsv', sep='\t')
    down_weights = pandas.read_csv(main_dir+'initial_weights/down.tsv', sep='\t')
    
    # Network
    network_df = pandas.read_csv('../../results/networks/MASLD_unified_directed_network_with_semantics.tsv', sep='\t')
    MASLD_network = igraph.Graph.DataFrame(network_df, directed=True, use_vids=False)
    network_nodes = MASLD_network.vs['name']
    
    # Network scores and subset for the 24 markers
    scores_df = pandas.read_csv(main_dir+'network_scores_per_sw.tsv', sep='\t',
                                index_col=0)
    scores_df = scores_df.loc[(scores_df != 0).any(axis=1),:]
    genes24_scores_df = scores_df.loc[genes24,:]

    sws = scores_df.columns.tolist()
    events = {}
    for gene in genes24:
        absolute_values = numpy.abs(genes24_scores_df.loc[gene,:].values)
        values = numpy.where(absolute_values > 0.2)[0]
        for v in values:
            sw = sws[v]
            direction = 'up' if genes24_scores_df.loc[gene,sw] > 0 else 'down'
            events.setdefault(gene, []).append([sw, direction])


    tmp_dir = main_dir+'Reactome_for_consensus_propagation_signatures/'
    files = os.listdir(tmp_dir)
    pathways_data = {}
    for f in files:
        sw = f.split('_')[0]
        direction = f.split('_')[1].split('.')[0]
        tmp_df = pandas.read_csv(tmp_dir+f, sep='\t', index_col=0)
        tmp = list(map(lambda x: x.split(' ')[-1], tmp_df.index.tolist()))
        tmp_df.loc[:,'term_id'] = tmp
        tmp_df = tmp_df.loc[tmp_df.loc[:, 'adj_p.value'] < 0.05,:]
        tmp = zip(tmp_df.term_id.tolist(), tmp_df.genes.tolist())
        for entry in tmp:
            pathways_data.setdefault(sw, {}).setdefault(direction, {}).update({entry[0]:entry[1].split(';')})

    
    reactome_graph = create_reactome_graph()
    root_node = reactome_graph.get_root()

    genes24_mapping_to_terms = {}
    genes24_mapping_to_genes = {}
    for gene, gene_events in events.items():
        for e in gene_events:
            e_sw = e[0]
            e_direction = e[1]
            try:
                reactome_dict = pathways_data[e_sw][e_direction]
            except KeyError:
                continue
            df = [[t,g] for t in reactome_dict.keys() for g in reactome_dict[t]]
            df = pandas.DataFrame(df, columns=['term', 'gene_symbol'])
            df = df.loc[df.gene_symbol == gene,:]
            terms_list = []
            for term in df.term.tolist():
                descs = reactome_graph.get_vertex_descendants(term)
                descs = [d['name'] for d in descs]
                if len(set(descs).intersection(df.term.tolist())):
                    pass
                else:
                    terms_list.append(term)
            genes24_mapping_to_terms.setdefault(gene, {}).update({(e_sw, e_direction):terms_list})
            for term in terms_list:
                genes24_mapping_to_genes.setdefault(gene, {}).setdefault((e_sw, e_direction), set()).update(reactome_dict[term])

    reference_df = []
    for gene, event in genes24_mapping_to_genes.items():
        for e, e_genes in event.items():
            terms_list = genes24_mapping_to_terms[gene][e]
            reference_df.append([gene, e[0], 
                                 int(e[0].replace('SW', '')), e[1], 
                                 ';'.join(list(e_genes)), 
                                 ';'.join(terms_list)])
            
    columns = ['gene', 'sw', 'sw_int', 'direction', 'genes', 'terms']
    reference_df = pandas.DataFrame(reference_df, columns=columns)    
    reference_df.sort_values(by='sw_int', ascending=True, inplace=True)
    reference_df.reset_index(inplace=True, drop=True)
    reference_df.to_csv(output_dir+'24_plasma_markers.tsv', sep='\t')
