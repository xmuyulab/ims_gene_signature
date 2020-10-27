from warnings import filterwarnings
filterwarnings('ignore')
import pandas as pd
import numpy as np
import os
from sklearn.metrics import roc_auc_score, roc_curve
from sklearn.preprocessing import MinMaxScaler


def calculate_auc(data,var='Ratio_score',sd=1,drop_sd=False):
    """Compute Area Under the Receiver Operating Characteristic Curve (ROC AUC)
    from GEP signatures.
    
    Parameters
    ----------
    data:
        pandas DataFrame, response and score must in data.columns.

    var:
        str, optional, default is None. Target scores.
    sd: 
        int, optional, default is 1, SD(stable disease) samples as responders, if 0 SD 
        samples as Non-responders.
    drop_sd:
        boolean, optional, default: False. If True remove SD sample when calculate AUC

    Returns
    -------
    auc:
        float, auc score
    """
    
    if drop_sd:
        data = data[data.response!=-1].copy()
    else:
        data['response'] = data['response'].replace(-1,sd)
    
    if var:
        assert var in data.columns
        scaler = MinMaxScaler(feature_range=(0, 1), copy=True)
        scaler.fit(data[[var]].astype(np.float64).values)
        score = pd.DataFrame(index=data.index, columns=[var], data=scaler.transform(data[[var]].values))
  
        auc_p = pd.merge(score,
                         data[['response']].copy(),
                         left_index=True, 
                         right_index=True)
        
        return roc_auc_score(auc_p['response'].values.astype(int), auc_p[var].values),\
               roc_curve(auc_p['response'].values.astype(int),auc_p[var].values)
    
    return None
