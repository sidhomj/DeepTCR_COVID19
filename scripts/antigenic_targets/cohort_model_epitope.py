import pandas as pd
import numpy as np
from DeepTCR.DeepTCR import DeepTCR_WF
import pickle

df_isb = pd.read_csv('isb_covid_seq.csv')
df_niaid = pd.read_csv('niaid_covid_seq.csv')
df = pd.concat([df_isb,df_niaid])
df = df.groupby(['peptide','sample','cohort','cell_type']).agg({'freq':'sum'}).reset_index()

seq = np.array(df['peptide'])
freq = np.array(df['freq'])
sample_id = np.array(df['sample'])
label_id = np.array(df['cohort'])


DTCR = DeepTCR_WF('cohort_model_epitope',device=3)
DTCR.Load_Data(beta_sequences=seq,freq=freq,sample_labels=sample_id,class_labels=label_id)


folds = 10
epochs_min = 25
size_of_net = 'medium'
num_concepts=64
hinge_loss_t = 0.3
train_loss_min=0.2
seeds = np.array(range(folds))
graph_seed = 0
#
DTCR.Monte_Carlo_CrossVal(folds=folds,epochs_min=epochs_min,size_of_net=size_of_net,num_concepts=num_concepts,
                          train_loss_min=train_loss_min,combine_train_valid=True,hinge_loss_t=hinge_loss_t,
                          seeds=seeds,graph_seed=graph_seed,subsample=100,batch_size=25)

with open('cohort_model_ep_preds.pkl','wb') as f:
    pickle.dump(DTCR.DFs_pred,f,protocol=4)

with open('cohort_model_ep_seq_preds.pkl', 'wb') as f:
    pickle.dump([DTCR.predicted,DTCR.beta_sequences,DTCR.lb], f, protocol=4)

df_out = pd.DataFrame()
df_out['seq'] = DTCR.beta_sequences
for ii,c in enumerate(DTCR.lb.classes_,0):
    df_out[c] = DTCR.predicted[:,ii]
df_out['sum'] = np.sum(DTCR.predicted,1)
df_out['cell_type'] = df['cell_type']
df_out = df_out[df_out['sum']>0.0]
df_out = df_out.groupby(['seq','cell_type']).agg({DTCR.lb.classes_[0]:'mean',DTCR.lb.classes_[1]:'mean'}).reset_index()
df_out.sort_values(by=DTCR.lb.classes_[0],inplace=True,ascending=False)
df_out = df_out[df_out['cell_type']=='CD8']

DTCR.Residue_Sensitivity_Logo(beta_sequences=np.array(df_out['seq'])[0:25],class_sel=DTCR.lb.classes_[0],
                              figsize=(5,10),background_color='black',Load_Prev_Data=False,
                              min_size=0.25)



