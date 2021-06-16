"""Creates Residue Sensitivity Logos for top 25 predictive sequences in ISB & NIH/NIAID data sets for the cohort-trained model."""

import pandas as pd
import numpy as np
from DeepTCR.DeepTCR import DeepTCR_WF
import pickle
import matplotlib.pyplot as plt

with open('cohort_model_ep_seq_preds.pkl', 'rb') as f:
    predicted, beta_sequences,lb = pickle.load(f)

df_isb = pd.read_csv('isb_covid_seq.csv')
df_niaid = pd.read_csv('niaid_covid_seq.csv')
df = pd.concat([df_isb,df_niaid])
df = df.groupby(['peptide','sample','cohort','cell_type']).agg({'freq':'sum'}).reset_index()


df[lb.classes_[0]] = predicted[:,0]
df[lb.classes_[1]] = predicted[:,1]

sel = 0
cl = lb.classes_[sel]
df.sort_values(by=cl,inplace=True,ascending=False)
cell_type = 'CD8'
df = df[df['cell_type']==cell_type]
df = df[df['cohort']=='COVID-19-ISB']
DTCR = DeepTCR_WF('cohort_model_epitope')
fig,ax = DTCR.Residue_Sensitivity_Logo(beta_sequences=np.array(df['peptide'])[0:25],
                              class_sel=cl,figsize=(5,10),background_color='black',Load_Prev_Data=False,
                              min_size=0.25)
fig.savefig('figures_ant/cohort_ep_'+cell_type+'_rsl_isb.png',dpi=1200,facecolor='black')

sel = 1
cl = lb.classes_[sel]
df.sort_values(by=cl,inplace=True,ascending=False)
cell_type = 'CD8'
df = df[df['cell_type']==cell_type]
df = df[df['cohort']=='COVID-19-NIH/NIAID']
DTCR = DeepTCR_WF('cohort_model_epitope')
fig,ax = DTCR.Residue_Sensitivity_Logo(beta_sequences=np.array(df['beta_sequences'])[0:25],
                              class_sel=cl,figsize=(5,10),background_color='black',Load_Prev_Data=False,
                              min_size=0.25)
fig.savefig('figures/cohort_ep_'+cell_type+'_rsl_niaid.png',dpi=1200,facecolor='black')





