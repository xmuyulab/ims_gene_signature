from warnings import filterwarnings
filterwarnings('ignore')
import pandas as pd
import numpy as np
import os

from calculate_ratio_score import calculate_ratio_score
from load_data import load_IO_data
from calculate_performance import calculate_auc


if __name__ == "__main__":

    IO_data = load_IO_data()
    IO_data = calculate_ratio_score(IO_data)
    
    # calculate AUC of Ratio score in datasets
    for ds in set(IO_data.data_set):
        auc, _ = calculate_auc(IO_data[IO_data.data_set==ds],drop_sd=True)
        print('data_set:{}, AUC: {}'.format(ds, auc))
