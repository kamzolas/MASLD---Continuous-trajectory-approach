#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas # 2.1.4
import json #  2.0.9
import numpy # 1.26.4

    
def worker_for_down_regulated_gene_sets(gene):
    mapping = {}
    # Genes as targets - force their upregulation
    genes_mechanisms_df = gene_chembl_activities_df.loc[(gene_chembl_activities_df.gene == gene) & 
                                                        gene_chembl_activities_df.action_type.isin(positive_keywords),:]
    groups = genes_mechanisms_df.groupby(by='gene')
    for gene, gene_df in groups:
        tmp = list(zip(gene_df.molecule_chembl_id, gene_df.action_type, gene_df.max_phase))
        mapping.update({gene:tmp})

    return mapping


def worker_for_up_regulated_gene_sets(gene):
    mapping = {}
    # Genes as targets - force their downregulation
    genes_mechanisms_df = gene_chembl_activities_df.loc[(gene_chembl_activities_df.gene == gene) & 
                                                        gene_chembl_activities_df.action_type.isin(negative_keywords),:]
    groups = genes_mechanisms_df.groupby(by='gene')
    for gene, gene_df in groups:
        tmp = list(zip(gene_df.molecule_chembl_id, gene_df.action_type, gene_df.max_phase))
        mapping.update({gene:tmp})
    
    return mapping


if __name__ == '__main__':
    
    output_dir = '../../results/ucam_sanyal/drugs/'
    network_analysis_dir = '../../results/ucam_sanyal/networks/network_analysis/'
    
    filename = output_dir+'gene_activities.tsv'
    gene_chembl_activities_df = pandas.read_csv(filename, sep='\t', index_col=0, low_memory=False)
    
    filename = output_dir+'chembl_drugs.tsv'
    chembl_gene_ids = pandas.read_csv(filename, sep='\t', index_col=0)
    
    filename = output_dir+'drug_indications.csv'
    drug_indications_df = pandas.read_csv(filename, sep=',', index_col=0)
    
    ensembl_ids = pandas.read_csv('../../data/ensembl_description.tsv', sep='\t')
    ensembl_ids = ensembl_ids.loc[:, ['external_gene_name', 'description']]
    ensembl_ids.drop_duplicates(inplace=True)
    
    positive_keywords = ['ACTIVATOR', 'AGONIST', 'STABILISER']
    negative_keywords = ['INHIBITOR', 'ANTAGONIST', 'NEGATIVE ALLOSTERIC MODULATOR', 'HYDROLYTIC ENZYME']
    keywords = list(set(positive_keywords).union(negative_keywords))
    
    drug_targets = gene_chembl_activities_df.gene.unique().tolist()
    
    filename = output_dir+'plasma_markers_chembl_ids.json'
    with open(filename, 'r') as fobj:
        drug_targets_chembl_ids = json.load(fobj)
    
    scores_df = pandas.read_csv(network_analysis_dir+'final_scores/network_nodes_scores_per_sw.tsv', sep='\t', index_col=0)
    drug_targets_scores_df = scores_df.loc[drug_targets,:]
    
    sws = drug_targets_scores_df.columns.tolist()
    events = {}
    for gene in drug_targets:
        absolute_values = numpy.abs(drug_targets_scores_df.loc[gene,:].values)
        values = numpy.where(absolute_values > 0)[0]
        for v in values:
            sw = sws[v]
            direction = 'up' if drug_targets_scores_df.loc[gene,sw] > 0 else 'down'
            events.setdefault(gene, []).append([sw, direction])


    drug_targets_pairs = {}
    for gene, gene_events in events.items():
        try:
            chembl_ids = drug_targets_chembl_ids[gene]
        except KeyError:
            continue
        for event in gene_events:
            sw = event[0]
            direction = event[1]
            key = gene + '_' + sw + '_' + direction
            if direction == 'down':
                mapping = worker_for_down_regulated_gene_sets(gene)
            else:
                mapping = worker_for_up_regulated_gene_sets(gene)
            drug_targets_pairs.update({key:mapping})
            

    data = []
    for key, tmp_dict in drug_targets_pairs.items():
        for gene, gene_entries in tmp_dict.items():
            for entry in gene_entries:
                data.append([key, gene, entry[0], entry[1], entry[2]])
                
    data = pandas.DataFrame(data)
    data.columns = ['event', 'target', 'chembl_drug_id', 'action_type', 'max_phase']
    data.target = data.target.apply(lambda x: 'L ribosomal proteins' if x.startswith('RPL') else x)
    data.target = data.target.apply(lambda x: 'S ribosomal proteins' if x.startswith('RPS') else x)
    data.drop_duplicates(inplace=True)
    data = pandas.merge(data, chembl_gene_ids, left_on='chembl_drug_id', right_on='molecule_chembl_id', how='left')
    data = pandas.merge(data, ensembl_ids, left_on='target', right_on='external_gene_name', how='left')
    data = data.loc[:, ['event', 'target', 'description', 'chembl_drug_id', 
                        'molecule_type', 'pref_name', 'action_type', 'max_phase']]
    data.description = data.description.apply(lambda x: x if pandas.isnull(x) else x.split(' [Source')[0])
    
    
    filename = output_dir+'drug_indications.csv'
    indications_df = pandas.read_csv(filename, sep=',', index_col=0)
    indications_df.drop('indication_refs', axis=1, inplace=True)
    df = pandas.merge(data, indications_df, left_on='chembl_drug_id', right_on='molecule_chembl_id', how='left')
    df.loc[:,'disease_stage_(sw)'] = df.event.apply(lambda x: int(x.split('_')[1].replace('SW','')) )
    
    df = df.loc[:, ['disease_stage_(sw)', 'event', 'target', 'description', 'chembl_drug_id', 'molecule_type',
       'pref_name', 'action_type', 'max_phase', 'efo_term', 'mesh_heading']]
    df.sort_values('disease_stage_(sw)', inplace=True, ascending=True)
    df.reset_index(inplace=True, drop=True)
    
    df.sort_values(by='target', inplace=True)
    
    df.to_csv(output_dir+'drug_repurposing_for_plasma_markers.tsv', sep='\t')
    
    
    
    
    
    
    
    
    
    
    