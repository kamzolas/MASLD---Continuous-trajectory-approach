#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import pandas # 2.1.4
import json #  2.0.9
from chembl_webresource_client.new_client import new_client # 0.10.9 (Python 3.8 environment)


if __name__ == '__main__':
    
    output_dir ='../../results/ucam_sanyal/drugs/'
    drug_indications = new_client.drug_indication
        
    df = pandas.read_csv(output_dir+'gene_activities.tsv', sep='\t')
    ids = df.molecule_chembl_id.unique().tolist()
    s = []
    for _id in ids:
        res = drug_indications.filter(molecule_chembl_id=_id)
        for r in res:
            s.append(r)
    df = pandas.DataFrame(s)
    df.to_csv(output_dir+'drug_indications.csv', sep=',')
