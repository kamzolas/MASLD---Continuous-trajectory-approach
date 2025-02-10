#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May  3 16:17:40 2024

@author: thodoris
"""

import pandas
import json



if __name__ == "__main__":
    
    output_dir = '../../results/finding_of_optimal_SWs_sequence/'
    pc_sorting_df = pandas.read_csv('../../data/PC1_sorted_samples.csv', sep=',')
    pc_sorting_df.columns = ['Sample name', 'PC1_value', 'PC2_value']
    pc_sorting_df.set_index('Sample name', inplace=True)
    
    metadata = pandas.read_csv('../../data/metadata.csv', sep=',', index_col=0)
    metadata.set_index('Sample name', inplace=True, drop=True)
    features = ['PC1_value', 'Steatosis', 'Ballooning', 'Fibrosis', 'NAS', 'Inflammation']
    metadata = metadata.join(pc_sorting_df, on='Sample name')
    
    with open(output_dir+'final_paths_after_filtering.json', 'r') as f:
        paths = json.load(f)

    
    # Calculate the correlations between PC1 and the other NAFLD variables    
    correlations_df = []
    indexes = []
    for param, paths_dict in paths.items():
        mean_feature_values = []
        for sw, samples in paths_dict.items():
            tmp = metadata.loc[samples, features].mean().to_frame()
            tmp.columns = [sw]
            mean_feature_values.append(tmp)
    
        mean_feature_values = pandas.concat(mean_feature_values, axis=1).transpose()
        correlations = mean_feature_values.corr(method='spearman')
        row = -correlations.loc['PC1_value',:].to_frame().transpose().iloc[:,1:]
        #row.loc[:,'N'] = len(paths_dict)

        row.index = [param]
        correlations_df.append(row)
    
    correlations_df = pandas.concat(correlations_df, axis=0)
    correlations_df.loc[:,'average_cor'] = correlations_df.loc[:,features[1:]].mean(axis=1)
    correlations_df.loc[:,'median_cor'] = correlations_df.loc[:,features[1:]].median(axis=1)
    correlations_df.sort_values(by='median_cor', ascending=False, inplace=True)
    correlations_df.to_csv(output_dir+'final_paths_spearman_correlations_with_MASLD_traits.tsv', sep='\t')
    
    
    
    # Load stats, join them with correlations and calculate the final factors
    stats_df = pandas.read_csv(output_dir+'final_paths_stats.tsv', sep='\t', index_col=0)
    #final_df = stats_df.join(correlations_df)
    final_df = stats_df.copy()
    final_df.loc[:,'inv_cv'] = 1/final_df.loc[:,'cv']
    final_df.loc[:,'pass_lower_thr'] = final_df.loc[:,'pass_lower_thr']/final_df.loc[:,'N']
    
    
    #final_df.loc[:,'f1_median_corr'] = final_df.loc[:,'median_cor']/final_df.loc[:,'median_cor'].max()
    final_df.loc[:,'f2_inv_cv'] = final_df.loc[:,'inv_cv']/final_df.loc[:,'inv_cv'].max()
    final_df.loc[:,'f3_median_degs'] = final_df.loc[:,'median']/final_df.loc[:,'median'].max()
    final_df.loc[:,'f4_pass_lower_thr'] = final_df.loc[:,'pass_lower_thr']/final_df.loc[:,'pass_lower_thr'].max()
    
    #features = ['f1_median_corr', 'f2_inv_cv', 'f3_median_degs', 'f4_pass_lower_thr']
    features = ['f2_inv_cv', 'f3_median_degs', 'f4_pass_lower_thr']
    final_df = final_df.loc[:, features]
    final_df.loc[:,'score'] = final_df.product(axis=1)
    final_df.sort_values(by='score', ascending=False, inplace=True)
    final_df.to_csv(output_dir+'final_paths_ranking_scores.tsv', sep='\t')
    
    
    
    
    
    
    