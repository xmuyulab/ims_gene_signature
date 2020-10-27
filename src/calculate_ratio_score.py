from warnings import filterwarnings
filterwarnings('ignore')
import pandas as pd
import numpy as np
import os
from scipy.stats import ttest_ind

from load_data import load_IO_data, nanostring_sig_gene


def data_preprocess(data:pd.DataFrame, ifn_set:list):
    
    data['ifn_set'] = data[ifn_set].mean(axis=1)
    pre_data = data[[i for i in data.columns if i not in['response','data_set','treatment','ifn_set']]].apply(lambda x:x-data['ifn_set']).copy()
    pre_data['ifn_set'] = data['ifn_set']
    clinical = data[['response','data_set','treatment']]
    return pre_data, clinical


def ttest_pvalues(pre_data,clinical, data_set=False,single_sided=True):

    data = pd.merge(pre_data,clinical,left_index=True,right_index=True)
    response_sample = data[data.response.isin([1,-1])].copy()
    
    nonresponse_sample = data[data.response==0].copy()
    
    if data_set:
        
        ds_result = []
        for d in set(clinical['data_set']):
            result = {}
            for i in list(pre_data.columns):
                _,pvalues = ttest_ind(response_sample[response_sample.data_set==d][i],
                                      nonresponse_sample[nonresponse_sample.data_set==d][i])
                result[i] = pvalues

            result = pd.DataFrame.from_dict(result,orient='index',columns=['{}_p_value'.format(d)])
            result['response_mean_{}'.format(d)] = response_sample[response_sample.data_set==d][pre_data.columns].mean()
            result['non-response_mean_{}'.format(d)] = nonresponse_sample[nonresponse_sample.data_set==d][pre_data.columns].mean()
            result['DE(r-nr)_{}'.format(d)] = result['response_mean_{}'.format(d)]-result['non-response_mean_{}'.format(d)]
            
            if single_sided==True:
                for i in result.index:
                    if result['DE(r-nr)_{}'.format(d)].loc[i]>0:
                        result['{}_p_value'.format(d)].loc[i]=1-0.5*result['{}_p_value'.format(d)].loc[i]
                    else:
                        result['{}_p_value'.format(d)].loc[i]=0.5*result['{}_p_value'.format(d)].loc[i]
            
            result['response_variance_{}'.format(d)] = response_sample[response_sample.data_set==d][pre_data.columns].var()
            result['non-response_variance_{}'.format(d)] = nonresponse_sample[nonresponse_sample.data_set==d][pre_data.columns].var()
            ds_result.append(result.T)
        return ds_result
    else:
        result = {}
        for i in list(pre_data.columns):
            _,pvalues = ttest_ind(response_sample[i],nonresponse_sample[i])
            result[i] = pvalues

        result = pd.DataFrame.from_dict(result,orient='index',columns=['p_value'])
        result['response_mean'] = response_sample[pre_data.columns].mean()
        result['non-response_mean'] = nonresponse_sample[pre_data.columns].mean()
        result['DE(r-nr)'] = result['response_mean']-result['non-response_mean']
        result['response_variance'] = response_sample[pre_data.columns].var()
        result['non-response_variance'] = nonresponse_sample[pre_data.columns].var()
    
        return result.sort_values('p_value',ascending=True)


def combined_p_pearson(order_table:pd.DataFrame):

    order_table['DE_mean'] = order_table[[i for i in order_table.columns if 'DE(r-nr)' in i]].min(axis=1)
    order_table['pearson_p'] = order_table[[i for i in order_table.columns if 'p_value' in i]].apply(lambda x:-1*sum([np.log10(1-i) for i in x.values]),axis=1)
    order_table['p_max'] = order_table[[i for i in order_table.columns if 'p_value' in i]].max(axis=1)

    return order_table


# Identification of the IMS genes list
def Identified_IMS_gene(train_set=['Riaz17', 'Gide19', 'Hugo16'], IFN_set=['IFNG', 'STAT1', 'CCR5', 'CXCL9', 'CXCL10', 'CXCL11', 'IDO1', 'PRF1', 'GZMA', 'HLA-DRA'], n=18):
    
    IFN_GAMMA_SET = IFN_set
    data = load_IO_data()
    # find IMS gene set in discovery dataset
    data = data[data.data_set.isin(train_set)].T.dropna().T
    
    # Skip Stable Disease (SD) sample
    data = data[data.response!=-1]
    
    # Normalize by IFN gamma gene set
    de_pre_data, de_clinical = data_preprocess(data,IFN_GAMMA_SET)
    
    candidate_gene = list(set(nanostring_sig_gene)&set(data.columns))
    print('candidate gene: {}'.format(len(candidate_gene)))

    #single sided t-test DE analysis
    ttest_result = ttest_pvalues(de_pre_data[candidate_gene],
                                 de_clinical[de_clinical.data_set.isin(train_set)],
                                 data_set=True, 
                                 single_sided=True)
    
    ttest_result = pd.concat(ttest_result).T
    
    #combined p-value using Pearson’s method
    order_table = combined_p_pearson(ttest_result)
    IMS_set = list(order_table.sort_values('pearson_p',ascending=True)[:n].index)
    #print(IMS_set)

    return IMS_set


def calculate_ratio_score(expr:pd.DataFrame, find_IMS:bool=False):

    if find_IMS:
        # Identification of the IMS genes
        IMS_set = Identified_IMS_gene()
    else:
        # result of find_IMS_set
        IMS_set = ['CCL8', 'VCAN', 'CCL2', 'BCAT1', 'ISG15', 'CD163', 'AXL', 'CCL13', 'COL6A3', 
                    'SIGLEC1', 'PDGFRB', 'IL10', 'STC1', 'ADAM12', 'OLFML2B', 'FAP', 'TWIST2', 'INHBA']
        
    gene_set ={'IFN_GAMMA_SET': ['IFNG', 'STAT1', 'CCR5', 'CXCL9', 'CXCL10', 'CXCL11', 'IDO1', 'PRF1', 'GZMA', 'HLA-DRA'],
               'House_keeping': ['ABCF1', 'DNAJC14', 'ERCC3', 'G6PD', 'GUSB', 'MRPL19', 'OAZ1', 'POLR2A', 'PSMC4',
                                 'PUM1', 'SDHA', 'SF3A1', 'STK11IP', 'TBC1D10B', 'TBP', 'TFRC', 'TLK2', 'TMUB2', 'UBB'],
               'IMS_SET': IMS_set}
    # Normalize by House keeping gene set
    expr['IFN_GAMMA_SET_score'] = expr[gene_set['IFN_GAMMA_SET']].mean(axis=1)-expr[gene_set['House_keeping']].mean(axis=1)
    expr['IMS_SET_score'] = expr[gene_set['IMS_SET']].mean(axis=1)-expr[gene_set['House_keeping']].mean(axis=1)
    expr['Ratio_score'] = expr['IFN_GAMMA_SET_score']-expr['IMS_SET_score']
    
    return expr


if  __name__ == "__main__":
    
    IMS_set = Identified_IMS_gene()
    print(IMS_set)
