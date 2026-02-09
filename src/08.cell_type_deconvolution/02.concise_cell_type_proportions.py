#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov  7 15:24:59 2025

@author: tkoutsandreas
"""

import pandas
import numpy
from matplotlib import pyplot as plt


if __name__ == "__main__":
    
    main_dir = '../../results/ucam_sanyal/cell_type_deconvolution/'
    data_dir = '../../data/ucam_sanyal/'
    
    filename = '../../data/cell_types_grouping/cell_types_names_and_features.tsv'
    cell_type_groups_df = pandas.read_csv(filename, sep='\t', index_col=0)
    main_classes = sorted(cell_type_groups_df.main_class.unique())
    groups_dict = {}
    for main_class in main_classes:
        tmp_df = cell_type_groups_df.loc[cell_type_groups_df.loc[:,'main_class'] == main_class,:]
        cell_type = tmp_df.cell_type.tolist()
        groups_dict.update({main_class:cell_type})
        
        
    filename = main_dir+'cell_proportions_per_sample.tsv'
    prop_per_sample_df = pandas.read_csv(filename, sep='\t', index_col=0)
    prop_per_sample_df = prop_per_sample_df.pivot_table(columns='cell_type',
                                                        values='consensus', 
                                                        index='sample')
    new_prop_per_sample_df = []
    for main_class in main_classes:
        cell_types = groups_dict[main_class]
        tmp = prop_per_sample_df.loc[:,cell_types].sum(axis=1).to_frame(main_class)
        new_prop_per_sample_df.append(tmp)
    new_prop_per_sample_df = pandas.concat(new_prop_per_sample_df, axis=1)
    new_prop_per_sample_df.to_csv(main_dir+'concise_cell_proportions_per_sample.tsv', sep='\t')
    
    
    samples_per_sw = {}
    F = open(data_dir+'sw_samples.csv', 'r')
    next(F)
    for line in F:
        fields = line.replace('\n', '').split(',')
        sw = fields[0]
        samples = fields[1].split(';')
        samples_per_sw.update({sw:samples})
    F.close()
    
    dfs = []
    sw_ints = list(range(1, len(samples_per_sw.keys())+1))
    for i in sw_ints:
        sw = 'SW_' + str(i)
        samples = samples_per_sw[sw]
        tmp_df = new_prop_per_sample_df.loc[new_prop_per_sample_df.index.isin(samples),:]
        tmp_df = tmp_df.mean(axis=0).to_frame(i)
        dfs.append(tmp_df)
    prop_per_sw_df = pandas.concat(dfs, axis=1).transpose()
    

    groups = [['Hepatocytes'],
              ['Endothelial cells', 'Fibroblasts', 'Cholangiocytes', 
               'Lymphoid cells', 'Macrophages', 'Myeloid cells']]
    
    colors_dict = dict(zip(cell_type_groups_df.main_class.tolist(),
                           cell_type_groups_df.main_color.tolist()))

    fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(16,5))
    
    axes[0].plot(prop_per_sw_df.index.tolist(), prop_per_sw_df.loc[:,groups[0]],
                 linewidth=3.5, label=groups[0][0], color=colors_dict[groups[0][0]])
    axes[0].set_ylim(0.45, 0.75)
    axes[0].tick_params(labelsize=12)
    axes[0].set_yticks(numpy.arange(0.45, 0.76, 0.05))
    axes[0].set_xticks(prop_per_sw_df.index.tolist(),labels=prop_per_sw_df.index.tolist(), size=12)
    
    axes[0].grid(alpha=0.25, linestyle='--')
    leg = axes[0].legend(fontsize=14, frameon=False, handlelength=1.5)
    for line in leg.get_lines():
        line.set_linewidth(5.0)
    axes[0].set_xlabel('Sliding window', fontsize=16)
    axes[0].set_ylabel('Proportion', fontsize=16)
    
    for index, cell_type in enumerate(groups[1]):    
        axes[1].plot(prop_per_sw_df.index.tolist(), prop_per_sw_df.loc[:,cell_type],
                     linewidth=3.5, label=cell_type, color=colors_dict[cell_type])
    axes[1].set_ylim(0, 0.2)
    axes[1].tick_params(labelsize=12)
    axes[1].set_yticks(numpy.arange(0, 0.21, 0.04))
    axes[1].set_xticks(prop_per_sw_df.index.tolist(),labels=prop_per_sw_df.index.tolist(), size=12)

    axes[1].grid(alpha=0.25, linestyle='--')
    leg = axes[1].legend(ncol=2, fontsize=14, frameon=False, handlelength=1.5)
    for line in leg.get_lines():
        line.set_linewidth(5.0)
    axes[1].set_xlabel('Sliding window', fontsize=16)
    axes[1].set_ylabel('Proportion', fontsize=16)
    plt.savefig(main_dir+'concise_cell_type_proportions.pdf')