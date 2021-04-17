import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.colors import LinearSegmentedColormap
import matplotlib
matplotlib.rc('font', family='sans-serif')

#get labels
df = pd.read_csv('../../Data/ImmuneCODE/ImmuneCODE-Repertoire-Tags-002.2.tsv',sep='\t')
df['who_ordinal_scale_bin'] = None
df['who_ordinal_scale_bin'][df['who_ordinal_scale'] <=4] = 'mild'
df['who_ordinal_scale_bin'][df['who_ordinal_scale'] > 4] = 'severe'
df['who_ordinal_scale'] = df['who_ordinal_scale_bin']
df['Age'] = df['Age'].str.split(' ',expand=True)[0]
df['Age'] = df['Age'].astype(float)

df['sample_name'] = df['sample_name']+'.tsv'
cols = np.array(list(df.columns))

#get labels meta
df_tags = pd.read_csv('../../Data/ImmuneCODE/tags_counts.csv')

time = 'days_from_recovery_to_sample'
label = 'icu_admit'
ds = 'COVID-19-HUniv12Oct'
df = df[df['Dataset'] == ds]
df.dropna(subset=[label,time], inplace=True)

#numeric label
sns.violinplot(data=df,y=label,x=time,cut=0,orient='h')
plt.ylabel(label,fontsize=16)
plt.xlabel(time,fontsize=16)
plt.title(ds)
plt.tight_layout()


