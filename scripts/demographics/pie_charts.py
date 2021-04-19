""" Demographics -Biological Sex"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.colors import LinearSegmentedColormap
import matplotlib
matplotlib.rc('font', family='sans-serif')
plt.style.use('ggplot')
def make_autopct(values):
    def my_autopct(pct):
        total = sum(values)
        val = int(round(pct*total/100.0))
        return '{p:.2f}%  ({v:d})'.format(p=pct,v=val) if pct > 0 else ''
    return my_autopct

#get labels
df = pd.read_csv('../../Data/ImmuneCODE/ImmuneCODE-Repertoire-Tags-002.2.tsv',sep='\t')
df['who_ordinal_scale_bin'] = None
df['who_ordinal_scale_bin'][df['who_ordinal_scale'] <=4] = 'mild'
df['who_ordinal_scale_bin'][df['who_ordinal_scale'] > 4] = 'severe'
df['who_ordinal_scale'] = df['who_ordinal_scale_bin']
df['Age'] = df['Age'].str.split(' ',expand=True)[0]
df['Age'] = df['Age'].astype(float)
df['severity'] = None
df['severity'][(df['icu_admit']==True) | (df['who_ordinal_scale']=='severe')] = 'severe'
df['severity'][(df['icu_admit']==False) | (df['who_ordinal_scale']=='mild')] = 'mild'

ds = ['COVID-19-ISB','COVID-19-NIH/NIAID']
if ds is not None:
    df = df[df['Dataset'].isin(ds)]


labels = ['Biological Sex','Racial Group','ethnicity','severity']
DFs = []
for label_sel in labels:
    c_list = []
    data_list = []
    for c in np.unique(df['Dataset']):
        df_c = df[df['Dataset'] == c]
        df_c.dropna(subset=[label_sel],inplace=True)
        data_list.append(df_c[label_sel])
        c_list.append(len(df_c)*[c])

    df_data = pd.DataFrame()
    df_data['c'] = np.hstack(c_list)
    df_data[label_sel] = np.hstack(data_list)

    #categorical label
    df_data = pd.DataFrame(df_data.groupby(['c',label_sel]).size())
    df_data.reset_index(inplace=True)
    df_data = pd.pivot(df_data,'c',label_sel,values=0)
    df_data.fillna(value=0,inplace=True)
    sum_rows = np.sum(df_data,1)
    DFs.append(df_data)

fig,ax = plt.subplots(len(DFs),2,figsize=(5,10))
for jj,df_data in enumerate(DFs,0):
    labels = list(df_data.columns)

    for ii,d in enumerate(ds):
        values = np.array(df_data.loc[d])
        w,l,p = ax[jj,ii].pie(values, labels=labels, autopct=make_autopct(values), startangle=90,textprops={'fontsize':8},
                              radius=3000)
        ax[jj,ii].axis('equal')
        for zz, t in enumerate(l, 0):
            if values[zz] == 0:
                t.set_text('')
        # ax[jj,ii].set_title(d)
        # for _ in l:
        #     _.set_fontsize(28)
plt.tight_layout()
plt.savefig('figures/sex.eps')