#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas
import os


if __name__ == '__main__':
    
    
    ###########################################################################
    #
    # CALCULATE THE CONSENSUS CELL TYPE PROPORTIONS PER SAMPLE
    #
    ###########################################################################
    
    # Get patient ids and organise them regarding their SW
    samples_per_sw = {}
    F = open('../../data/sw_samples.csv', 'r')
    next(F)
    for line in F:
        fields = line.replace('\n', '').split(',')
        sw = fields[0]
        samples = fields[1].split(';')
        samples_per_sw.update({sw:samples})
    F.close()
    ##samples_per_sw = dict()
    #for col in tmp_df.columns:
    #    samples = tmp_df.loc[:, col]
    #    samples = list(filter(lambda x: not pandas.isnull(x), samples))
    #    samples_per_sw.update({col:samples})
    
    
    # Load the matrices of cell type deconvolution (results per patient)
    files = os.listdir('../../results/cell_type_deconvolution/matrices/')
    matrices = {}
    for f in files:
        key = f.split('.')[0]
        matrix = pandas.read_csv('../../results/cell_type_deconvolution/matrices/'+f,
                                 sep='\t', index_col=0)
        matrices.update({key:matrix})
    
    
    # Collect the deconvolution results fron the different tools and create a 
    # consesus score for each cell type and patient
    tools = ['DWLS', 'bseqsc', 'EpiDISH', 'CIBERSORT', 'BayesPrism', 'FARDEEP']
    dfs = []
    for tool in tools:
        matrix = matrices[tool]
        tmp_df = matrix.melt(value_vars=matrix.columns, 
                             var_name='cell_type', 
                             value_name=tool, 
                             ignore_index=False)
        tmp_df.loc[:, 'sample'] = tmp_df.index.tolist()
        multi_index = pandas.MultiIndex.from_frame(tmp_df.loc[:,['sample','cell_type']],
                                                   names=['sample','cell_type'])
        tmp_df.set_index(multi_index, inplace=True)
        tmp_df = tmp_df.loc[:,[tool]]
        dfs.append(tmp_df)
    tmp_df = pandas.concat(dfs, axis=1, ignore_index=False)
    tmp_df.to_csv('../../results/cell_type_deconvolution/cell_proportions_raw.tsv', sep='\t')
    '''
     		DWLS 	bseqsc 	EpiDISH 	CIBERSORT 	BayesPrism
    sample 	cell_type 					
    Sample.1 	B_cell 	0.010641 	0.0 	0.000000 	0.000000 	0.002648
    Sample.2 	B_cell 	0.011627 	0.0 	0.000000 	0.027076 	0.001087
    Sample.3 	B_cell 	0.038191 	0.0 	0.000000 	0.000000 	0.003574
    Sample.4 	B_cell 	0.005401 	0.0 	0.000000 	0.000000 	0.001243
        
    '''
    tmp_df = tmp_df.mean(axis=1).to_frame()
    tmp_df.columns = ['consensus']
    tmp_df.loc[:, 'cell_type'] = tmp_df.index.get_level_values(1)
    tmp_df = tmp_df.droplevel(level=1)
    '''
        consensus 	cell_type
    sample 		
    Sample.1 	0.002658 	B_cell
    Sample.2 	0.007958 	B_cell
    Sample.3 	0.008353 	B_cell
    Sample.4 	0.001329 	B_cell
    '''
    
    # Normalize the consensus scores per patient to have sum equal to one
    dfs = []
    for sample, sample_df in tmp_df.groupby(level=0):
        sample_df.loc[:, 'consensus'] = sample_df.consensus/sample_df.consensus.sum()
        sample_df.reset_index(inplace=True, drop=False)
        dfs.append(sample_df)
    prop_per_sample_df = pandas.concat(dfs, axis=0, ignore_index=True)
    prop_per_sample_df.pivot_table(values='consensus', index='sample', columns='cell_type')
    prop_per_sample_df.loc[:,'sample'] = prop_per_sample_df.loc[:,'sample'].map(lambda x: x.replace('.',' '))
    prop_per_sample_df.to_csv('../../results/cell_type_deconvolution/cell_proportions_per_sample.tsv', sep='\t')
    


    ###########################################################################
    #
    # CALCULATE THE CONSENSUS CELL TYPE PROPORTIONS PER SW
    #
    ###########################################################################
    dfs = []
    sw_ints = list(range(1, len(samples_per_sw.keys())+1))
    for i in sw_ints:
        sw = 'SW_' + str(i)
        samples = samples_per_sw[sw]
        #if i < 13:
        #    next_samples = samples_per_sw['SW_' + str(i+1)]
        #    samples = list(set(samples).difference(next_samples))
        tmp_df = prop_per_sample_df.loc[prop_per_sample_df.loc[:,'sample'].isin(samples),:]
        tmp_df = tmp_df.groupby(by='cell_type').consensus.mean().to_frame().reset_index()
        tmp_df.loc[:, 'sw'] = [sw]*tmp_df.shape[0]
        dfs.append(tmp_df)
    
    prop_per_sw_df = pandas.concat(dfs, axis=0, ignore_index=True)
    prop_per_sw_df = prop_per_sw_df.pivot_table(values='consensus', index='sw', columns='cell_type')
    prop_per_sw_df = prop_per_sw_df.loc[:, prop_per_sw_df.max(axis=0) > 0.005]
    prop_per_sw_df.index = prop_per_sw_df.index.map(lambda x: int(x.split('_')[1]))
    prop_per_sw_df = prop_per_sw_df.transpose()
    prop_per_sw_df.iloc[1:,:] = prop_per_sw_df.iloc[1:,:]/prop_per_sw_df.iloc[1:,:].sum(axis=0)
    prop_per_sw_df = prop_per_sw_df.transpose()    
    prop_per_sw_df.sort_values(by='sw', inplace=True)
    prop_per_sw_df.to_csv('../../results/cell_type_deconvolution/cell_type_proportions_per_sw.tsv', sep='\t')


