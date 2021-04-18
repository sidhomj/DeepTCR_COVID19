"""Creates Residue Sensitivity Logos for top 25 predictive sequences in ISB & NIH/NIAID data sets."""
import pandas as pd
import numpy as np
from DeepTCR.DeepTCR import DeepTCR_WF
import pickle
import matplotlib.pyplot as plt

#isb
with open('../supervised/isb_ft_pred.pkl','rb') as f:
    _,_,df_isb = pickle.load(f)
df_isb['cohort'] = 'COVID-19-ISB'

df_isb.sort_values(by='pred',inplace=True,ascending=False)
DTCR = DeepTCR_WF('isb_model')
models = ['model_'+str(x) for x in np.random.choice(100,25,replace=False)]
fig,ax = DTCR.Residue_Sensitivity_Logo(beta_sequences=np.array(df_isb['beta_sequences'])[0:25],models=models,
                              class_sel='severe',figsize=(5,10),background_color='black',Load_Prev_Data=False,
                              min_size=0.25)
fig.savefig('figures/isb_rsl.png',dpi=1200,facecolor='black')


#niaid
with open('../supervised/niaid_ft_pred.pkl','rb') as f:
    _,_,df_niaid = pickle.load(f)
df_niaid['cohort'] = 'COVID-19-NIH/NIAID'

df_niaid.sort_values(by='pred',inplace=True,ascending=False)
DTCR = DeepTCR_WF('niaid_model')
models = ['model_'+str(x) for x in np.random.choice(100,25,replace=False)]
fig,ax = DTCR.Residue_Sensitivity_Logo(beta_sequences=np.array(df_niaid['beta_sequences'])[0:25],models=models,
                              class_sel=True,figsize=(5,10),background_color='black',Load_Prev_Data=False,
                              min_size=0.25)
fig.savefig('figures/niaid_rsl.png',dpi=1200,facecolor='black')

