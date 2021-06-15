"""Creates Residue Sensitivity Logos for top 25 predictive sequences in ISB & NIH/NIAID data sets for the cohort-trained model."""

import pandas as pd
import numpy as np
from DeepTCR.DeepTCR import DeepTCR_WF
import pickle
import matplotlib.pyplot as plt


with open('cohort_model_seq_preds.pkl', 'rb') as f:
    predicted, beta_sequences,lb = pickle.load(f)

df = pd.DataFrame()
df['beta_sequences'] = beta_sequences
df[lb.classes_[0]] = predicted[:,0]
df[lb.classes_[1]] = predicted[:,1]

sel = 0
cl = lb.classes_[sel]
df.sort_values(by=cl,inplace=True,ascending=False)
DTCR = DeepTCR_WF('cohort_model')
fig,ax = DTCR.Residue_Sensitivity_Logo(beta_sequences=np.array(df['beta_sequences'])[0:25],
                              class_sel=cl,figsize=(5,10),background_color='black',Load_Prev_Data=False,
                              min_size=0.25)
fig.savefig('figures/cohort_rsl_isb.png',dpi=1200,facecolor='black')

sel = 1
cl = lb.classes_[sel]
df.sort_values(by=cl,inplace=True,ascending=False)
DTCR = DeepTCR_WF('cohort_model')
fig,ax = DTCR.Residue_Sensitivity_Logo(beta_sequences=np.array(df['beta_sequences'])[0:25],
                              class_sel=cl,figsize=(5,10),background_color='black',Load_Prev_Data=False,
                              min_size=0.25)
fig.savefig('figures/cohort_rsl_niaid.png',dpi=1200,facecolor='black')





