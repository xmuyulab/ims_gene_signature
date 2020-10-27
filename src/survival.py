import pandas as pd
import numpy as np
import os
from matplotlib import pyplot as plt

from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test
from lifelines.plotting import add_at_risk_counts
from sklearn.metrics import roc_auc_score, roc_curve

from calculate_ratio_score import calculate_ratio_score
from load_data import survival_data, load_IO_data

def cutoff_youdens_j(data, var='Ratio_score', sd=1, rm_sd=True):
    """Calculate Youden's J statistic for selecting the optimum cut-off point of Ratio score in
    dataset.
        J = sensitivity + specificity -1
    
    Parameters
    ----------
    data:
        pandas DataFrame, score and response must in columns.
    sd:
        int, optional, default is 1, SD(stable disease) samples as responders, if 0 SD 
        samples as Non-responders.
    rm_sd:
        boolean, optional, default: True. If True remove SD sample when calculate Youden's J
    
    Returns
    -------
    j_ordered:
        pandas DataFrame, Youden's J statistic.
    """
    
    assert var in data.columns
    assert 'response' in data.columns
    score = data[['response',var]].copy()
    
    if rm_sd:
        score = score[score.response!=-1].copy()
    
    else:
        if sd==1:
            score['response'] = score.response.replace(-1,1)

        elif sd==0:
            score['response'] = score.response.replace(-1,0)

        else:
            print('error')
    
    fpr, tpr, threshold = roc_curve(score['response'].values.astype(int), score[var].values)
    j_scores = tpr-fpr
    j_ordered = sorted(zip(j_scores,threshold))
    j_ordered = pd.DataFrame(columns=['youden_j','threshold'],data=j_ordered)
    return j_ordered


def sub_draw(inp, threshold, ds,var='Ratio_score', time_type='OS'):
    
    # Youdens J score
    threshold = threshold-1e-5
    
    data = inp.copy()
    print(data.shape)
    data['label'] = data[var].apply(lambda x:[0,1][x>=threshold])
    fig = plt.figure(figsize=(10,8))
    ax = plt.subplot(111)
    kmf0 = KaplanMeierFitter()

    kmf0.fit(data[data['label']==0][time_type],
           data[data['label']==0]['status'],
           label=0)
    kmf0.plot(ax=ax,
             show_censors=True,
             at_risk_counts=False,
             ci_show=False)
    
    kmf1 = KaplanMeierFitter()
    kmf1.fit(data[data['label']==1][time_type],
           data[data['label']==1]['status'],
           label=1)
    kmf1.plot(ax=ax,
             at_risk_counts=False,
             show_censors=True, 
             ci_show=False)
    
    # There are equivalent
    add_at_risk_counts(kmf0, kmf1)

    T1 = data[data['label']==1][time_type]
    T2 = data[data['label']==0][time_type]

    E1 = data[data['label']==1]['status']
    E2 = data[data['label']==0]['status']

    results = logrank_test(T1, T2, event_observed_A=E1, event_observed_B=E2)
    results.print_summary()

    if not os.path.exists('../result/'):
        os.makedirs('../result/')
    plt.savefig('../result/{}_survival_{}.pdf'.format(ds,time_type))
    return data


if  __name__ == "__main__":
    IO_data = load_IO_data()
    IO_data = calculate_ratio_score(IO_data)
    for  d in survival_data.keys():
        tmp = cutoff_youdens_j(IO_data[IO_data.data_set==d], 
                               rm_sd=True,
                               var='Ratio_score')
        threshold = tmp.sort_values('youden_j',ascending=False).values[0][1]
        surv_data = pd.merge(IO_data[IO_data.data_set==d][['Ratio_score']],
                             survival_data[d],
                             left_index=True,
                             right_index=True)
        surv_data_ = sub_draw(surv_data,threshold=threshold, ds=d)
