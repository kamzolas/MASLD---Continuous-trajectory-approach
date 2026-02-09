#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import json
import pandas # 2.1.4
import igraph # 0.11.4
import os


###############################################################################
# Description:
###############################################################################
# This scirpt reads the results of network propagation and generates a score 
# matrix which quantifies the level of differentiation of each node in the 
# MASLD network in each sliding window. The score of a node for a specific 
# sliding window is calculated as the difference between the respective 
# log-transformed pvalues in the up and down-regulated signatures.
# Outputs:
# - A tsv file with all the scores
###############################################################################


if __name__ == '__main__':
    
    ###########################################################################
    # Inputs
    ###########################################################################
    main_dir = '../../results/ucam_sanyal/networks/'
    output_dir = main_dir+'network_analysis/final_scores/'
    os.makedirs(output_dir, exist_ok=True)
    
    # Signatures per sw    
    filename =  main_dir+'network_analysis/consensus_propagation_signatures/random_sw-lapl_norm_weight.json'
    with open(filename, 'r') as f:
        sw_signatures = json.load(f)
    
    # log-transformed p-values (up & down)
    filename = main_dir+'network_analysis/consensus_propagation_signatures/random_sw-lapl_norm_weight_up.tsv'
    scores_up = pandas.read_csv(filename, sep='\t')
    filename = main_dir+'/network_analysis/consensus_propagation_signatures/random_sw-lapl_norm_weight_down.tsv'
    scores_down = pandas.read_csv(filename, sep='\t')
    
    # Network nodes
    network_df = pandas.read_csv(main_dir+'MASLD_unified_undirected_network_with_semantics.tsv', sep='\t')
    masld_network = igraph.Graph.DataFrame(network_df, directed=True, use_vids=False)
    network_nodes = masld_network.vs['name']

    # Open data frames to store the results
    network_nodes_up = masld_network.vs['name']
    network_nodes_up = pandas.DataFrame(data={'weight':[0.0]*len(network_nodes_up)}, index=network_nodes_up)
    network_nodes_down = masld_network.vs['name']
    network_nodes_down = pandas.DataFrame(data={'weight':[0.0]*len(network_nodes_down)}, index=network_nodes_down)
    
    for sw, signature in sw_signatures.items():
        new_col = [0.0]*len(network_nodes)
        network_nodes_up.loc[:,sw] = new_col
        for node in signature['up']:
            network_nodes_up.loc[node, sw] = scores_up.loc[node, sw]
        
    for sw, signature in sw_signatures.items():
        new_col = [0.0]*len(network_nodes)
        network_nodes_down.loc[:,sw] = new_col
        for node in signature['down']:
            network_nodes_down.loc[node, sw] = -scores_down.loc[node, sw]
    
    del network_nodes_up['weight']
    del network_nodes_down['weight']
    scores_df = network_nodes_up + network_nodes_down        
    scores_df.to_csv(output_dir+'network_nodes_scores_per_sw.tsv', sep='\t')