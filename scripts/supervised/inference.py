import pandas as pd
import numpy as np
from DeepTCR.DeepTCR import DeepTCR_WF
import pickle
from sklearn.metrics import roc_auc_score
import seaborn as sns
from copy import deepcopy
from scipy.stats import spearmanr

model = 'isb'
model = 'niaid'

data = 'isb'
# data = 'niaid'
# data = 'huniv'

DTCR = DeepTCR_WF('data_1000')
DTCR.Get_Data('../../Data/ImmuneCODE/repertoires',
              Load_Prev_Data=True,
              aa_column_beta=3,
              v_beta_column=49,
              d_beta_column=35,
              j_beta_column=42,
              count_column=4,
              type_of_data_cut='Num_Seq',
              data_cut=500)
#out = np.unique(DTCR.sample_id,return_counts=True)

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

#select cohort (optional)
if data == 'niaid':
    ds = ['COVID-19-NIH/NIAID']
elif data == 'huniv':
    ds = ['COVID-19-HUniv12Oct']
elif data == 'isb':
    ds = ['COVID-19-ISB']
else:
    ds = None

#ds = ['COVID-19-NIH/NIAID','COVID-19-HUniv12Oct']
if ds is not None:
    df = df[df['Dataset'].isin(ds)]

#drop any columns with all nans
df = df.dropna(axis=1,how='all')

cols = np.array(list(df.columns))

if data == 'isb':
    label_sel = 'who_ordinal_scale_bin'
else:
    label_sel = 'icu_admit'
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

DTCR = DeepTCR_WF(model+'_model',device=4)
models = ['model_'+str(x) for x in np.random.choice(50,25,replace=False)]
models = None
DTCR.Sample_Inference(sample_labels=sample_id,
                      beta_sequences=beta_sequences,
                      counts=counts,batch_size=25,
                      models=models)
if model == 'isb':
    df_preds = deepcopy(DTCR.Inference_Pred_Dict['severe'])
    df_preds['label'] = df_preds['Samples'].map(label_dict)
    ds_dict = dict(zip(df_sel['sample_name'],df_sel['Dataset']))
    df_preds['ds'] = df_preds['Samples'].map(ds_dict)

    print(roc_auc_score(df_preds['label'],df_preds['Pred']))
    sns.violinplot(data=df_preds,x='label',y='Pred',cut=0)
else:
    df_preds = deepcopy(DTCR.Inference_Pred_Dict[True])
    df_preds['label'] = df_preds['Samples'].map(label_dict)
    ds_dict = dict(zip(df_sel['sample_name'],df_sel['Dataset']))
    df_preds['ds'] = df_preds['Samples'].map(ds_dict)

    print(roc_auc_score(df_preds['label'],df_preds['Pred']))
    sns.violinplot(data=df_preds,x='label',y='Pred',cut=0)


