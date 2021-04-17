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

ds = ['COVID-19-HUniv12Oct']
if ds is not None:
    df = df[df['Dataset'].isin(ds)]

label_sel = 'days_from_recovery_to_sample'
df.dropna(subset=['icu_admit',label_sel],inplace=True)
sns.violinplot(data=df,y='icu_admit',x=label_sel,orient='h')
plt.ylabel('icu_admit')
plt.yticks(fontsize=12)
plt.gca().set_facecolor('white')
plt.tight_layout()
plt.savefig('figures/days_recov_sample.eps')




