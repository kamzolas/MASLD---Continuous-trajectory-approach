#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas


###############################################################################
# Description:
###############################################################################
# The final paths in sw_graph aew ranked using the product of three max-scaled
# criteria: (i) median DEG count, (ii) number of edges exceeding the lower DEG 
# threshold, and (iii) coefficient of variation of DEG counts along the path.
# - Outputs:
# A tsv file which contains the scores for all the potential SW paths.
###############################################################################


if __name__ == "__main__":
    
    output_dir = '../../results/ucam_sanyal/finding_of_optimal_SWs_sequence/'
    output_filename = 'final_paths_ranking_scores.tsv'
    
    # Load stats, join them with correlations and calculate the final factors
    stats_df = pandas.read_csv(output_dir+'final_paths_stats.tsv', sep='\t', index_col=0)
    filtered_indexes = list(filter(lambda x: x.startswith('sw_graph_8_24_03'), stats_df.index.tolist()))
    stats_df = stats_df.loc[stats_df.index.isin(filtered_indexes),:]
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
    final_df.to_csv(output_dir+output_filename, sep='\t')
    
    
    
    
    
    
    