#!/usr/bin/env python3
# -*- coding: utf-8 -*-


from chembl_webresource_client.new_client import new_client
import json


def retrieve_target_ids(tmp_dict):
    genes = tmp_dict['genes']
    tfs = tmp_dict['tfs']
    all_genes = list(set(genes).union(tfs))
    
    attrs = ['target_chembl_id', 'pref_name', 'target_type']
    targets_raw_data = {}
    targets_chembl_ids = {}
    for gene in all_genes:    
        res = target.filter(target_synonym__iexact=gene, organism='Homo sapiens').only(attrs)
        targets_raw_data.update({gene:res})
        for r in res:
            targets_chembl_ids.setdefault(gene, []).append(r['target_chembl_id'])
    
    return targets_chembl_ids


if __name__ == '__main__':
    
    output_dir ='../../results/networks/network_analysis/24_markers/'
    filename = output_dir+'24_plasma_markers_extended.json'
    with open(filename, 'r') as fobj:
        markers_dict = json.load(fobj)
    
    target = new_client.target
    target_ids = {}
    for key, tmp_dict in markers_dict.items():
        tmp_dict = retrieve_target_ids(tmp_dict)
        target_ids.update({key: tmp_dict})
    
    filename = '24_plasma_markers_chembl_ids.json'
    with open(output_dir+filename, 'w') as fobj:
        json.dump(target_ids, fobj)