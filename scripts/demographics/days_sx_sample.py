"""Demographics - Days between symptom onset and sample."""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.colors import LinearSegmentedColormap
import matplotlib
plt.style.use('ggplot')
matplotlib.rc('font', family='sans-serif')

#get labels
df = pd.read_csv('../../Data/ImmuneCODE/ImmuneCODE-Repertoire-Tags-002.2.tsv',sep='\t')
df['who_ordinal_scale_bin'] = None
df['who_ordinal_scale_bin'][df['who_ordinal_scale'] <=4] = 'mild'
df['who_ordinal_scale_bin'][df['who_ordinal_scale'] > 4] = 'severe'
df['who_ordinal_scale'] = df['who_ordinal_scale_bin']
df['Age'] = df['Age'].str.split(' ',expand=True)[0]
df['Age'] = df['Age'].astype(float)

ds = ['COVID-19-ISB','COVID-19-NIH/NIAID']
if ds is not None:
    df = df[df['Dataset'].isin(ds)]

label_sel = 'days_from_symptom_onset_to_sample'
label_sel = 'days_from_diagnosis_to_sample'

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

#numeric label
sns.violinplot(data=df_data,y='c',x=label_sel,cut=0,orient='h')
plt.ylabel('')
plt.yticks(fontsize=12)
plt.gca().set_facecolor('white')
plt.tight_layout()
plt.savefig('figures/days_sx_sample.eps')
plt.savefig('figures/days_dx_sample.eps')




