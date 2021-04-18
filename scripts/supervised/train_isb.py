"""Fits repertoire classifier on COVID=19-ISB data set."""
import pandas as pd
import numpy as np
from DeepTCR.DeepTCR import DeepTCR_WF
import pickle

model = 'isb'

DTCR = DeepTCR_WF('data_1000')
classes = ['COVID-19-ISB','COVID-19-NIH-NIAID','COVID-19-HUniv12Oct']
DTCR.Get_Data('../../Data/ImmuneCODE/repertoires_org',
              Load_Prev_Data=True,
              aa_column_beta=3,
              v_beta_column=49,
              d_beta_column=35,
              j_beta_column=42,
              count_column=4,
              data_cut=1000,
              type_of_data_cut='Num_Seq',
              classes=classes)

beta_sequences = DTCR.beta_sequences
counts = DTCR.counts
sample_id = DTCR.sample_id
v_beta = DTCR.v_beta
d_beta = DTCR.d_beta
j_beta = DTCR.j_beta

#get labels
df = pd.read_csv('../../Data/ImmuneCODE/ImmuneCODE-Repertoire-Tags-002.2.tsv',sep='\t')
df['who_ordinal_scale_bin'] = None
df['who_ordinal_scale_bin'][df['who_ordinal_scale'] <=4] = 'mild'
df['who_ordinal_scale_bin'][df['who_ordinal_scale'] > 4] = 'severe'
df['who_ordinal_scale'] = df['who_ordinal_scale_bin']
df['Age'] = df['Age'].str.split(' ',expand=True)[0]
df['Age'] = df['Age'].astype(float)
df['sample_name'] = df['sample_name']+'.tsv'

if model == 'niaid':
    ds = ['COVID-19-NIH/NIAID']
elif model == 'isb':
    ds = ['COVID-19-ISB']

if ds is not None:
    df = df[df['Dataset'].isin(ds)]

#drop any columns with all nans
df = df.dropna(axis=1,how='all')
cols = np.array(list(df.columns))

if model == 'niaid':
    label_sel = 'icu_admit'
elif model == 'isb':
    label_sel = 'who_ordinal_scale'

df_sel = df.dropna(subset=[label_sel])

label_dict = dict(zip(list(df_sel[cols[2]]),list(df_sel[label_sel])))
keep = np.array(list(df_sel['sample_name']))

idx = np.isin(sample_id,keep)

sample_id = sample_id[idx]
beta_sequences = beta_sequences[idx]
counts = counts[idx]
v_beta = v_beta[idx]
d_beta = d_beta[idx]
j_beta = j_beta[idx]

label_id = np.array(list(map(label_dict.get,sample_id)))

DTCR = DeepTCR_WF(model+'_model',device=3)
DTCR.Load_Data(beta_sequences=beta_sequences,counts=counts,sample_labels=sample_id,class_labels=label_id)

folds = 100
epochs_min = 25
size_of_net = 'small'
num_concepts=12
hinge_loss_t = 0.1
train_loss_min=0.1
seeds = np.array(range(folds))
graph_seed = 0
#
DTCR.Monte_Carlo_CrossVal(folds=folds,epochs_min=epochs_min,size_of_net=size_of_net,num_concepts=num_concepts,
                          train_loss_min=train_loss_min,combine_train_valid=True,hinge_loss_t=hinge_loss_t,
                          seeds=seeds,graph_seed=graph_seed,subsample=100)

df_preds = pd.DataFrame()
df_preds['sample'] =  DTCR.sample_id
df_preds['beta_sequences'] = DTCR.beta_sequences
df_preds['pred'] = DTCR.predicted[:,1]
df_preds['freq'] = DTCR.freq

if model == 'niaid':
    DTCR.DFs_pred[True].to_csv('niaid_preds.csv', index=False)
    with open('niaid_ft_pred.pkl', 'wb') as f:
        pickle.dump([DTCR.features, DTCR.predicted,df_preds], f, protocol=4)
elif model =='isb':
    DTCR.DFs_pred['severe'].to_csv('isb_preds.csv',index=False)
    with open('isb_ft_pred.pkl','wb') as f:
        pickle.dump([DTCR.features,DTCR.predicted,df_preds],f,protocol=4)





