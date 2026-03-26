#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas # 2.1.4
import json #  2.0.9
from chembl_webresource_client.new_client import new_client # 0.10.9 (Python 3.8 environment)


if __name__ == "__main__":
    
    output_dir = '../../results/ucam_sanyal/drugs/'
    filename = output_dir+'plasma_markers_chembl_ids.json'
    with open(filename, 'r') as fobj:
        markers_dict = json.load(fobj)
        
    mechanism = new_client.mechanism
    gene_chembl_activities = []
    
    #attrs = ['molecule_chembl_id', 'action_type', 'activity_comment', 'assay_description', 'target_chembl_id']
    for gene, chembl_ids in markers_dict.items():
        for chembl_id in chembl_ids:
            res = mechanism.filter(target_chembl_id=chembl_id)
            for entry in res:
                entry.update({'gene':gene})
                gene_chembl_activities.append(entry)

    gene_chembl_activities = pandas.DataFrame(gene_chembl_activities)
    gene_chembl_activities.to_csv(output_dir+'gene_activities.tsv', sep='\t')
    
    
    
    
    
    
    
    