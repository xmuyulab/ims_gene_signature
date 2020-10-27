import pandas as pd
import numpy as np
import json
import os, sys


def load_IO_data():
    """Load sample clinical data and RNA-seq data, RNA-seq data is after performance 
    of quantile normalization and log2 transformation
    
    Returns
    -------
    data:
        pandas DataFrame, clinical information and gene expression.
    """
    data_path = '../data/'
    cli_data, exp_data = [], []

    for i in os.listdir(data_path):
        
        #if i in ['mel_van_exp_data.csv','mel_van_cli_data.csv']:
        #    continue
        if ('mel' in i) or ('gas' in i) :
            if 'cli_data.csv' in i:
                cli_data.append(pd.read_csv(os.path.join(data_path,i),index_col=0))
                #print(i,pd.read_csv(data_path+i,index_col=0).shape)
            if 'exp_data.csv' in i:
                # log2 
                exp = pd.read_csv(os.path.join(data_path,i), index_col=0).T
                if exp.max().max() > 500:
                    exp = np.log2(exp+1)   
                exp_data.append(exp)
    
    clinical = pd.concat(cli_data)
    expression = pd.concat(exp_data)
    # we focus on pre-treatment
    return pd.merge(clinical[clinical.treatment=='pre'],expression, left_index=True, right_index=True)


# load cancer immune-related genes curated in Nanostring’s IO 360 panel
with open('../data/NanoString_gene_list.txt','r') as f:
    nanostring_sig_gene = list(f.read().split('\n'))

survival_data = {'Hugo16':pd.read_csv('../data/mel_hugo16_survival_data.csv',index_col=0),
                 'Riaz17':pd.read_csv('../data/mel_bms038_survival_data.csv',index_col=0),
                 'Gide19': pd.read_csv('../data/mel_gide19_survival_data.csv',index_col=0),
                 'VanAllen': pd.read_csv('../data/mel_van_survival_data.csv',index_col=0),
                 'Liu19':pd.read_csv('../data/mel_liu_survival_data.csv',index_col=0)}
