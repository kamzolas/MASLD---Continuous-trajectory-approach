#!/usr/bin/env python3
# -*- coding: utf-8 -*-


from chembl_webresource_client.new_client import new_client
import json
import pandas


if __name__ == "__main__":
    
    output_dir = '../../results/networks/network_analysis/24_markers/'
    
    filename = output_dir+'24_plasma_markers_chembl_ids.json'
    with open(filename, 'r') as fobj:
        markers_dict = json.load(fobj)
        
    #filename = output_dir+'chembl_approved_drugs.tsv'
    #approved_drugs_df = pandas.read_csv(filename, sep='\t', index_col=0)
    
    mechanism = new_client.mechanism
    
    gene_chembl_activities = []
    examined_keys = []
    examined_genes = []
    
    #attrs = ['molecule_chembl_id', 'action_type', 'activity_comment', 'assay_description', 'target_chembl_id']
    for key, targets_dict in list(markers_dict.items()):
        print(key)
        if key in examined_keys:
            continue
        for gene, chembl_ids in targets_dict.items():
            print(gene)
            if gene in examined_genes:
                continue
            for chembl_id in chembl_ids:
                res = mechanism.filter(target_chembl_id=chembl_id)
                for entry in res:
                    entry.update({'gene':gene})
                    gene_chembl_activities.append(entry)
            examined_genes.append(gene)
        examined_keys.append(key)

    gene_chembl_activities = pandas.DataFrame(gene_chembl_activities)
    gene_chembl_activities.to_csv(output_dir+'gene_activities.tsv', sep='\t')
    
    
    
    
    
    
    
    
    