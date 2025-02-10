#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import pandas
import os


if __name__ == '__main__':
    
    
    main_dir = '../../results/networks/'
    
    mapping_df = pandas.read_csv(main_dir+'MASLD_nodes_without_semantics_filtered.csv', 
                                 sep=',', dtype=str, index_col=0)
    
    files = os.listdir(main_dir+'semantic_similarities/')
    files = [main_dir+'semantic_similarities/'+f for f in files]
    
    for f in files:
        matrix = pandas.read_csv(f, index_col=0, sep=',')
        str_index = [str(i) for i in matrix.index.tolist()]
        tmp_df = pandas.DataFrame(data={'str_index':str_index, 
                                        'columns':matrix.columns.tolist()})
        
        s = sum(tmp_df.loc[:,'str_index'] == tmp_df.loc[:,'columns'])
        print(matrix.shape, s)
        
        tmp_df = pandas.merge(tmp_df, mapping_df, left_on='str_index', 
                              right_on='entrez_id', how='inner')
        
        new_matrix = matrix.loc[:, tmp_df.entrez_id.tolist()]
        new_matrix.columns = tmp_df.gene_symbol.tolist()
        new_matrix = new_matrix.transpose()
        new_matrix = new_matrix.loc[:, [int(i) for i in tmp_df.entrez_id.tolist()]]
        new_matrix.columns = tmp_df.gene_symbol.tolist()
        
        #matrix.index = tmp_df.gene_symbol.tolist()
        #matrix.columns = tmp_df.gene_symbol.tolist()
        #print(new_matrix.index[1:5], new_matrix.columns[1:5])
        
        os.remove(f)
        f = f.replace('tsv', 'csv')
        new_matrix.to_csv(f, sep=',')
        