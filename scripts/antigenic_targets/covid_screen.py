import seaborn as sns
import numpy as np
import pandas as pd
import pickle
import matplotlib.pyplot as plt
from scipy.stats import spearmanr, rankdata,fisher_exact
from utils import *
from sklearn.metrics import roc_auc_score, roc_curve
import matplotlib
plt.style.use('ggplot')
matplotlib.rc('font', family='sans-serif')


#mira
db = 'mira'
cdr3_col = 'beta_sequences'
epitope_col = 'peptide'
categories = ['orf_name','peptide']

file = '../../Data/ImmuneCODE/mira_data_parsed_ci.csv'
df = pd.read_csv(file,dtype=object,engine='python')
df = Process_Seq(df,cdr3_col)
df.drop_duplicates(subset=[cdr3_col,epitope_col],inplace=True)
df['cell_type'] = 'CD8'
df_cd8 = df

file = '../../Data/ImmuneCODE/mira_data_parsed_cii.csv'
df = pd.read_csv(file,dtype=object,engine='python')
df = Process_Seq(df,cdr3_col)
df.drop_duplicates(subset=[cdr3_col,epitope_col],inplace=True)
df['cell_type'] = 'CD4'
df_cd4 = df

df = pd.concat([df_cd8,df_cd4])

#isb
with open('../supervised/isb_ft_pred.pkl','rb') as f:
    _,_,df_isb = pickle.load(f)
df_isb['cohort'] = 'COVID-19-ISB'
df_isb['covid+'] = None
df_isb['covid+'][df_isb['beta_sequences'].isin(df['beta_sequences'])] = 'COVID (+)'
df_isb['covid+'][~df_isb['beta_sequences'].isin(df['beta_sequences'])] = 'COVID (-)'
df_isb = pd.merge(df_isb,df,how='left')

#df_isb = df_isb[df_isb['cell_type']!='CD8']
l = 'covid+'
sns.violinplot(data=df_isb,x=l,y='pred',cut=0,order=['COVID (-)','COVID (+)'])
plt.gca().set_facecolor('white')
plt.gca().spines['left'].set_color('black')
plt.gca().spines['bottom'].set_color('black')
plt.xlabel('')
plt.ylabel('Prediction',fontsize=18)
plt.xticks(fontsize=18)
plt.gca().axhline(.90,linestyle='--',color='black')
plt.tight_layout()
plt.savefig('figures/covid_preds_isb.eps')


l = 'orf_name'
diff = []
for o in np.unique(df_isb[l].dropna()):
    df_temp = df_isb[df_isb[l]==o]
    df_temp = df_temp.groupby(['cell_type']).agg({'pred':'mean'})
    try:
        val_cd8 = df_temp.loc['CD8']['pred']
    except:
        val_cd8 = 0

    try:
        val_cd4 = df_temp.loc['CD4']['pred']
    except:
        val_cd4 = 0

    diff.append(np.abs(val_cd8-val_cd4))

order = np.unique(df_isb[l].dropna())[np.argsort(diff)]
sns.violinplot(data=df_isb,x=l,y='pred',cut=0,order=order,hue='cell_type')
plt.gca().set_facecolor('white')
plt.gca().spines['left'].set_color('black')
plt.gca().spines['bottom'].set_color('black')
plt.xticks(rotation=90)
plt.ylabel('Prediction',fontsize=18)
plt.tight_layout()
plt.xlabel('')
plt.legend(loc='lower right')
labels = [item.get_text() for item in plt.gca().get_xticklabels()]
labels_new = []
for l in labels:
    labels_new.append(l.replace(' ', '\n'))
plt.gca().set_xticklabels(labels_new)
plt.tight_layout()
plt.savefig('figures/covid_preds_isb_cd84.eps')


#niaid
with open('../supervised/niaid_ft_pred.pkl','rb') as f:
    _,_,df_niaid = pickle.load(f)
df_niaid['cohort'] = 'COVID-19-NIH/NIAID'
df_niaid['covid+'] = None
df_niaid['covid+'][df_niaid['beta_sequences'].isin(df['beta_sequences'])] = 'COVID (+)'
df_niaid['covid+'][~df_niaid['beta_sequences'].isin(df['beta_sequences'])] = 'COVID (-)'
df_niaid = pd.merge(df_niaid,df,how='left')

l = 'covid+'
sns.violinplot(data=df_niaid,x=l,y='pred',cut=0,order=['COVID (-)','COVID (+)'])
plt.gca().set_facecolor('white')
plt.gca().spines['left'].set_color('black')
plt.gca().spines['bottom'].set_color('black')
plt.tight_layout()
plt.xlabel('')
plt.ylabel('Prediction',fontsize=18)
plt.xticks(fontsize=18)
plt.gca().axhline(.90,linestyle='--',color='black')
plt.savefig('figures/covid_preds_niaid.eps')


l = 'orf_name'
sns.violinplot(data=df_niaid,x=l,y='pred',cut=0,order=order,hue='cell_type')
plt.gca().set_facecolor('white')
plt.gca().spines['left'].set_color('black')
plt.gca().spines['bottom'].set_color('black')
plt.xticks(rotation=90)
plt.ylabel('Prediction',fontsize=18)
plt.tight_layout()
plt.xlabel('')
plt.legend(loc='lower right')
labels = [item.get_text() for item in plt.gca().get_xticklabels()]
labels_new = []
for l in labels:
    labels_new.append(l.replace(' ', '\n'))
plt.gca().set_xticklabels(labels_new)
plt.tight_layout()
plt.savefig('figures/covid_preds_niaid_cd84.eps')

df_isb['call'] = df_isb['pred'] > 0.90
df_isb['count'] = 1
df_isb_cont = pd.pivot_table(data=df_isb,index='covid+',columns=['call'],values=['count'],aggfunc='sum')
df_isb_cont = df_isb_cont.loc[['COVID (-)','COVID (+)']]
OR, p_val = fisher_exact(df_isb_cont)

df_niaid['call'] = df_niaid['pred'] > 0.90
df_niaid['count'] = 1
df_niaid_cont = pd.pivot_table(data=df_niaid,index='covid+',columns=['call'],values=['count'],aggfunc='sum')
df_niaid_cont = df_niaid_cont.loc[['COVID (-)','COVID (+)']]
OR, p_val = fisher_exact(df_niaid_cont)

df_isb_covid = df_isb[df_isb['covid+']=='COVID (+)']
df_isb_covid.sort_values(by='pred',ascending=False,inplace=True)
df_isb_covid.to_csv('isb_covid_seq.csv',index=False)

df_niaid_covid = df_niaid[df_niaid['covid+']=='COVID (+)']
df_niaid_covid.sort_values(by='pred',ascending=False,inplace=True)
df_niaid_covid.to_csv('niaid_covid_seq.csv',index=False)



