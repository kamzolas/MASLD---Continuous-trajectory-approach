#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 17 12:28:27 2024

@author: thodoris
"""

from chembl_webresource_client.new_client import new_client
import pandas


if __name__ == '__main__':

    output_dir = '../../results/networks/network_analysis/24_markers/'
    molecule = new_client.molecule
    attrs = ['molecule_chembl_id', 'pref_name', 'molecule_type']
    approved_drugs4 = molecule.filter(max_phase=4).only(attrs)
    approved_drugs_df4 = pandas.DataFrame(approved_drugs4)
    
    attrs = ['molecule_chembl_id', 'pref_name', 'molecule_type']
    approved_drugs3 = molecule.filter(max_phase=3).only(attrs)
    approved_drugs_df3 = pandas.DataFrame(approved_drugs3)    
    
    attrs = ['molecule_chembl_id', 'pref_name', 'molecule_type']
    approved_drugs2 = molecule.filter(max_phase=2).only(attrs)
    approved_drugs_df2 = pandas.DataFrame(approved_drugs2)      
    
    approved_drugs_df = pandas.concat([approved_drugs_df4,
                                       approved_drugs_df3,
                                       approved_drugs_df2])
    
    approved_drugs_df.to_csv(output_dir+'chembl_drugs.tsv', sep='\t')
    