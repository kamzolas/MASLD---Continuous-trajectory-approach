#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas # 2.1.4
import numpy # 1.26.4
import collections
import os
from scipy.stats import gaussian_kde # 1.11.4
from scipy.io import mmread # 1.11.4

###############################################################################
# Description:
###############################################################################
# This scirpt is used to generate pseudo busk samples given the cell type 
# proportions of real samples, estimated through cell type deconvolution, and
# the reference counts matrix of single cell data. A random sample is generated
# given that it belongs to a specific sliding window, while its library size
# is being randomly selected within [20e6, 50e6]. For this sample, all the 
# respective (real) cell type proportions are retrieved to create a kernel for
# each cell type and finally produce a random vector of cell type proportions. 
# This vector is used as input to randomly select cells from the reference 
# counts matrix and create the pseudo bulk sample by adding cells gradually.
# The workflow generates a random number of pseudo bulk samples for each sliding 
# window, while the whole process runs 30 times. 
# Outputs:
# - 30 folders with count matrices for the pseudo bulk samples.
###############################################################################


def kernel_construction(samples_data, condition, cell_type, weighted=True, 
                        control_factor=0.1):

    tmp_values = samples_data.loc[samples_stratification[condition], cell_type].values
    tmp_mean = tmp_values.mean()
    
    global_values = samples_data.loc[:, cell_type].values
    global_mean = global_values.mean()
    
    if weighted == True:
        weights = global_mean*control_factor + numpy.abs((global_mean - tmp_values) + (global_mean - tmp_mean))
    else:
        weights = numpy.ones(len(tmp_values))
    
    kde = gaussian_kde(tmp_values, weights=weights)
    return kde


if __name__ == "__main__":

    data_dir = '../../data/'
    weighted_kernels = False
    #control_factor = 0.1
    Nr = 10000
    Ntrials = 30
    low_counts = 20e6
    high_counts = 50e6
    output_dir = '../../results/ucam_sanyal/pseudo_bulk_analysis/'
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        
    # A. Single cell dataset
    counts_matrix = mmread(data_dir+'sc_liver_cell_atlas/MASLD_selected_qc_counts_matrix.mtx')
    counts_matrix = counts_matrix.tocsr()
    
    # B. QC cells
    qc_cells_df = pandas.read_csv(data_dir+'sc_liver_cell_atlas/selected_qc_cells.tsv',
                                  sep='\t', header=None)
    qc_cells_df.columns = ['cell']
    cells_annotation_df = pandas.read_csv(data_dir+'sc_liver_cell_atlas/cells_annotation.tsv', 
                                          sep='\t', index_col=0)
    qc_cells_df = qc_cells_df.merge(cells_annotation_df, left_on='cell', 
                                    right_on='cell', how='inner')
    new_col = qc_cells_df.loc[:,'final_annotation'].apply(lambda x: x.replace(' ', '_'))
    qc_cells_df.loc[:,'final_annotation'] = new_col
    
    # C. QC genes
    qc_genes_df = pandas.read_csv(data_dir+'sc_liver_cell_atlas/MASLD_qc_genes.tsv', 
                                  sep='\t', header=None)
    qc_genes = qc_genes_df.iloc[:,0].tolist()
    
    # D. Bulk RNA-seq - metadata
    metadata = pandas.read_csv(data_dir+'metadata.csv', sep=',', index_col=0)
    metadata.set_index('Sample name', inplace=True)
    
    # E. Bulk RNA-seq samples - cell proportions
    filename = '../../results/ucam_sanyal/cell_type_deconvolution/cell_proportions_per_sample.tsv'
    cell_type_proportions = pandas.read_csv(filename, sep='\t', index_col=0)
    cell_type_proportions = cell_type_proportions.pivot_table(index='sample', 
                                                              columns='cell_type', 
                                                              values='consensus')
    
    # F. Change the names of cell types in cell_type_proportions df to be 
    # consistent with the cell types annotation 
    tmp_mapping = {'B_cell': 'B_cells', 'Basophils': 'Basophils',
                   'Cholangiocytes': 'Cholangiocytes',
                   'Circulating_NK_NKT': 'Circulating_NK/NKT',
                   'Endothelial': 'Endothelial_cells',
                   'Fibroblasts': 'Fibroblasts', 'Hepatocytes': 'Hepatocytes',
                   'Macrophages': 'Macrophages', 'Mig_cDCs': 'Mig_cDCs',
                   'Monocytes': 'Mono+mono_derived_cells', 
                   'Neutrophils': 'Neutrophils', 'Plasma_cell': 'Plasma_cells',
                   'Resident_NK': 'Resident_NK', 'T_cell': 'T_cells', 
                   'cDC1s': 'cDC1s', 'cDC2s': 'cDC2s', 'pDCs': 'pDCs'}
    cell_type_proportions.columns = [tmp_mapping[a] for a in cell_type_proportions.columns]
    cell_type_proportions.drop('Mig_cDCs', axis=1, inplace=True)
    
    # G. Organize the samples in SW groups, so as to create the cell 
    # type proportion kernels for each SW
    samples_data = cell_type_proportions.join(metadata, how='inner')
    samples_stratification_df = pandas.read_csv(data_dir+'sw_samples.csv', sep=',')
    samples_stratification = {}
    for index, row in samples_stratification_df.iterrows():
        samples = row.samples.split(';')
        samples_stratification.update({row.sw: samples})
    
    # Kernels construction
    numpy.random.seed(1234)
    kernels = {}
    distributions = {}
    for condition, tmp_samples in samples_stratification.items():
        for cell_type in cell_type_proportions.columns:
            kde = kernel_construction(samples_data, condition, cell_type, weighted_kernels)
            random_points = kde.resample(size=Nr)
            random_points = random_points[random_points >= 0]
            kernels.setdefault(condition, {}).update({cell_type:kde})
            distributions.setdefault(condition, {}).update({cell_type:random_points})
        
    for n in range(Ntrials):
        tmp_output_dir = output_dir+'pseudo_bulk_RNAseq/'+str(n)+'/'
        os.makedirs(tmp_output_dir, exist_ok=True)
        for sw in distributions.keys():
            dataset = []
            # Number of pseudo samples that will be generated
            Nsamples = numpy.random.randint(6,20,1)[0]
            for s in range(Nsamples): # for each random sample
                prop_df = []
                for cell_type, random_points in distributions[sw].items():
                    tmp_prop = numpy.random.choice(random_points, 1)[0] # random cell type proportion
                    prop_df.append([cell_type, tmp_prop])
                prop_df = pandas.DataFrame(prop_df, columns=['cell_type', 'proportion'])
                prop_df.proportion = prop_df.proportion/prop_df.proportion.sum() # scaling
                prop_df.set_index('cell_type', drop=False, inplace=True)
                prop_df.drop('cell_type', axis=1, inplace=True)
    
                # random library size: [20e6, 50e50] counts
                lib_size_thr = numpy.random.uniform(low=low_counts, high=high_counts, size=1).astype(int)[0]
                lib_size = 0
                batch_size = 100
    
                sample_cell_types = prop_df.index.tolist()
                sample_prop = prop_df.proportion.tolist()    
                sample_cells = []
                while lib_size < lib_size_thr:
                    # randomly select 'batch_times' cell type values given the sample cell type proportions
                    batch_cell_types = numpy.random.choice(sample_cell_types, size=batch_size, p=sample_prop)
                    batch_cell_types = list(batch_cell_types)
                    # count how many times each cell type is in the batch_cell_types list
                    batch_counts = collections.Counter(batch_cell_types)
                    batch_cells = []
                    for cell_type, number in batch_counts.items():
                        # for each cell type pick randomly 'number' cells and save them to the list 'batch_cells'
                        cell_type_cells = qc_cells_df.loc[qc_cells_df.final_annotation == cell_type,:].index.tolist()
                        random_cells = list(numpy.random.choice(cell_type_cells, size=number, replace=True))
                        batch_cells.extend(random_cells)
                    lib_size = lib_size + counts_matrix[:,batch_cells].toarray().sum()
                    # all the batches of cells are concatenated in the 'sample_cells' list
                    sample_cells.extend(batch_cells)
                # use the 'sample_cells' list to subset the initial counts matrix and create the pseudo sample
                raw_counts = counts_matrix[:,sample_cells].toarray().sum(axis=1).reshape(1,-1)
                dataset.append(raw_counts)
    
            final_dataset = numpy.concatenate(dataset)
            final_dataset = pandas.DataFrame(final_dataset.transpose())
            final_dataset.index = qc_genes
            final_dataset.columns = [sw+'_sample_'+str(i) for i in range(Nsamples)]
            filename = tmp_output_dir+sw+'_raw_counts.tsv'
            final_dataset.to_csv(filename, sep='\t')
    
            print('pseudo samples generation', str(n), sw, str(Nsamples))    
        
    
    
    
    
    