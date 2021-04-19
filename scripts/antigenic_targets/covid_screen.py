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
df_isb['covid+'] = df_isb['beta_sequences'].isin(df['beta_sequences'])
df_isb = pd.merge(df_isb,df,how='left')

#df_isb = df_isb[df_isb['cell_type']!='CD8']
l = 'covid+'
# l = 'orf_name'
sns.violinplot(data=df_isb,x=l,y='pred',cut=0,order=df_isb.groupby([l]).agg({'pred':'mean'}).sort_values(by='pred').index)
plt.xticks(rotation=90)
plt.gca().set_facecolor('white')
plt.gca().spines['left'].set_color('black')
plt.gca().spines['bottom'].set_color('black')
plt.xlabel('')
plt.ylabel('Prediction')
plt.gca().axhline(.90,linestyle='--')
plt.tight_layout()


l = 'orf_name'
sns.violinplot(data=df_isb,x=l,y='pred',cut=0,order=df_isb.groupby([l]).agg({'pred':'mean'}).sort_values(by='pred').index,hue='cell_type')
plt.xticks(rotation=90)
plt.tight_layout()
plt.xlabel('')
plt.ylabel('Prediction')
plt.legend(loc='lower right')


#niaid
with open('../supervised/niaid_ft_pred.pkl','rb') as f:
    _,_,df_niaid = pickle.load(f)
df_niaid['cohort'] = 'COVID-19-NIH/NIAID'
df_niaid['covid+'] = df_niaid['beta_sequences'].isin(df['beta_sequences'])
df_niaid = pd.merge(df_niaid,df,how='left')
l = 'covid+'
# l = 'orf_name'
# sns.violinplot(data=df_niaid,x=l,y='pred',cut=0,order=df_niaid.groupby([l]).agg({'pred':'mean'}).sort_values(by='pred').index)
sns.violinplot(data=df_niaid,x=l,y='pred',cut=0,order=[False,True])
plt.xticks(rotation=90)
plt.gca().set_facecolor('white')
plt.gca().spines['left'].set_color('black')
plt.gca().spines['bottom'].set_color('black')
plt.tight_layout()
plt.xlabel('')
plt.ylabel('Prediction')
plt.gca().axhline(.90,linestyle='--')


l = 'orf_name'
sns.violinplot(data=df_niaid,x=l,y='pred',cut=0,order=df_niaid.groupby([l]).agg({'pred':'mean'}).sort_values(by='pred').index,hue='cell_type')
plt.xticks(rotation=90)
plt.tight_layout()
plt.xlabel('')
plt.ylabel('Prediction')
plt.legend(loc='lower right')


#AUC curves
dfs = [df_isb,df_niaid]
ds = ['ISB','NIH/NIAID']
fig,ax = plt.subplots()
ax.plot([0, 1], [0, 1], color='navy', linestyle='--')
for df,d in zip(dfs,ds):
    score = roc_auc_score(df['covid+'],df['pred'])
    fpr,tpr,_ = roc_curve(df['covid+'],df['pred'])
    ax.plot(fpr,tpr,label=d+' (%0.2f)' % score)
    ax.spines['bottom'].set_color('black')
    ax.spines['left'].set_color('black')
    ax.set_facecolor('white')
ax.legend(prop={'size': 10}, loc='lower right', facecolor='white')
ax.set_xlabel('False Positive Rate', fontsize=14)
ax.set_ylabel('True Positive Rate', fontsize=14)
plt.tight_layout()

df_isb['call'] = df_isb['pred'] > 0.90
df_isb['count'] = 1
df_isb_cont = pd.pivot_table(data=df_isb,index='covid+',columns=['call'],values=['count'],aggfunc='sum')
OR, p_val = fisher_exact(df_isb_cont)

df_niaid['call'] = df_niaid['pred'] > 0.90
df_niaid['count'] = 1
df_niaid_cont = pd.pivot_table(data=df_niaid,index='covid+',columns=['call'],values=['count'],aggfunc='sum')
OR, p_val = fisher_exact(df_niaid_cont)



