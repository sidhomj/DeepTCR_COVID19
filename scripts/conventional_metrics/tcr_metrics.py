"""Plots showing associations between conventional TCR metrics of diversity/abundance and severe illness."""
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests
import matplotlib
plt.style.use('ggplot')
matplotlib.rc('font', family='sans-serif')

df_metrics = pd.read_csv('../../Data/ImmuneCODE/samples.tsv',sep='\t')
df = pd.read_csv('../../Data/ImmuneCODE/ImmuneCODE-Repertoire-Tags-002.2.tsv',sep='\t')
df['who_ordinal_scale_bin'] = None
df['who_ordinal_scale_bin'][df['who_ordinal_scale'] <=4] = 'mild'
df['who_ordinal_scale_bin'][df['who_ordinal_scale'] > 4] = 'severe'
df['who_ordinal_scale'] = df['who_ordinal_scale_bin']
df['Age'] = df['Age'].str.split(' ',expand=True)[0]
df['Age'] = df['Age'].astype(float)
df['Dataset'] = df['Dataset'].str.split('COVID-19-',expand=True)[1]

df['severe'] = None
df['severe'][(df['icu_admit']==True) | (df['who_ordinal_scale']=='severe')] = 'severe'
df['severe'][(df['icu_admit']==False) | (df['who_ordinal_scale']=='mild')] = 'mild'


df_merge = pd.merge(df,df_metrics,on=['sample_name'])

ds = ['NIH/NIAID','ISB']
label_sel = 'severe'

df_merge = df_merge[df_merge['Dataset'].isin(ds)]
df_merge.dropna(subset=[label_sel],inplace=True)

metrics = df_metrics.columns[1:10]
fig,ax = plt.subplots(3,3,figsize=(10,10))
ax = np.ndarray.flatten(ax)
for ii,c in enumerate(range(9),0):
    sns.violinplot(data=df_merge,x='Dataset',y=metrics[ii],cut=0,
                   ax = ax[ii],hue=label_sel)
    ax[ii].set_xlabel('')
    ax[ii].set_ylabel(metrics[ii],fontsize=12)
    ax[ii].spines['bottom'].set_color('black')
    ax[ii].spines['left'].set_color('black')
    ax[ii].set_facecolor('white')
    ax[ii].tick_params(axis='x', labelsize=14)
    if ii !=0:
        ax[ii].get_legend().remove()
    else:
        ax[ii].get_legend().get_frame().set_linewidth(0.0)
        ax[ii].get_legend().set_title(None)
plt.tight_layout()
plt.savefig('figures/tcr_metrics.eps')

#plot comparing cohorts
fig,ax = plt.subplots(3,3,figsize=(10,10))
ax = np.ndarray.flatten(ax)
for ii,c in enumerate(range(9),0):
    sns.violinplot(data=df_merge,x='Dataset',y=metrics[ii],cut=0,
                   ax = ax[ii])
    ax[ii].set_xlabel('')
    ax[ii].set_ylabel(metrics[ii],fontsize=12)
    ax[ii].spines['bottom'].set_color('black')
    ax[ii].spines['left'].set_color('black')
    ax[ii].set_facecolor('white')
    ax[ii].tick_params(axis='x', labelsize=14)
plt.tight_layout()
plt.savefig('figures/tcr_metrics_cohorts.eps')


#stats
metric_list = []
ds_list = []
p_val = []
for m in metrics:
    for d in ds:
        df_temp = df_merge[df_merge['Dataset']==d]
        _,p = mannwhitneyu(df_temp[df_temp[label_sel]=='severe'][m],df_temp[df_temp[label_sel]=='mild'][m])
        metric_list.append(m)
        ds_list.append(d)
        p_val.append(p)

df_stats = pd.DataFrame()
df_stats['metric'] = metric_list
df_stats['ds'] = ds_list
df_stats['p_val'] = p_val
_,df_stats['corr_p_val'],_,_ = multipletests(p_val,method='fdr_bh')
df_stats.sort_values(by='corr_p_val',inplace=True)
df_stats.to_csv('figures/tcr_stats.csv',index=False)

#stats at cohort level
metric_list = []
p_val = []
for m in metrics:
    _,p = mannwhitneyu(df_merge[df_merge['Dataset']=='ISB'][m],df_merge[df_merge['Dataset']=='NIH/NIAID'][m])
    metric_list.append(m)
    p_val.append(p)

df_stats = pd.DataFrame()
df_stats['metric'] = metric_list
df_stats['p_val'] = p_val
_,df_stats['corr_p_val'],_,_ = multipletests(p_val,method='fdr_bh')
df_stats.sort_values(by='corr_p_val',inplace=True)
df_stats.to_csv('figures/tcr_stats_cohort.csv',index=False)