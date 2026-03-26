#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas # 2.1.4
import json #  2.0.9
from chembl_webresource_client.new_client import new_client # 0.10.9 (Python 3.8 environment)


if __name__ == '__main__':
    
    # Inputs
    main_dir = '../../results/ucam_sanyal/networks/network_analysis/'
    filename1 =  main_dir+'consensus_propagation_signatures/random_sw-lapl_norm_weight.json'
    filename2 = '../../data/194_proteo_transcriptomic_signature.txt'
    output_dir = '../../results/ucam_sanyal/drugs/'

    # Results from network analysis - selection of deregulated genes
    with open(filename1, 'r') as f:
        sw_signatures = json.load(f)
    all_genes = set()
    for sw, sw_dict in sw_signatures.items():
        for direction, gene_list in sw_dict.items():
            all_genes.update(gene_list)
    
    # Pre-defined set of plasma markers
    F = open(filename2, 'r')
    plasma_markers = []
    for entry in F.readlines():
        gene = entry.strip().split(';')[0]
        plasma_markers.append(gene)
    plasma_markers = list(set(plasma_markers))
    
    # Intersection of plasma markers and deregulated genes in the MASLD network
    selected_markers = list((set(plasma_markers).intersection(all_genes)))
    
    target = new_client.target
    attrs = ['target_chembl_id', 'pref_name', 'target_type']
    targets_chembl_ids = {}
    for gene in selected_markers:    
        res = target.filter(target_synonym__iexact=gene, organism='Homo sapiens').only(attrs)
        for r in res:
            targets_chembl_ids.setdefault(gene, []).append(r['target_chembl_id'])
    
    filename = output_dir+'plasma_markers_chembl_ids.json'
    with open(filename, 'w') as fobj:
        json.dump(targets_chembl_ids, fobj)
