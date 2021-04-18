"""Creates figures to visualize deep learning model performance."""
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedKFold
from sklearn.preprocessing import LabelEncoder,OneHotEncoder, StandardScaler, MinMaxScaler
from sklearn.metrics import roc_curve,roc_auc_score
import matplotlib
plt.style.use('ggplot')
matplotlib.rc('font', family='sans-serif')

DFs = []

#niaid
df_dtcr = pd.read_csv('niaid_preds.csv')
df_dtcr = df_dtcr.groupby(['Samples']).agg({'y_pred':'mean','y_test':'mean'}).reset_index()
df_dtcr['model'] = 'dtcr'
df_dtcr['cohort'] = 'NIH/NIAID'
df_dtcr_niaid = df_dtcr
DFs.append(df_dtcr_niaid)

#isb
df_dtcr = pd.read_csv('isb_preds.csv')
df_dtcr = df_dtcr.groupby(['Samples']).agg({'y_pred':'mean','y_test':'mean'}).reset_index()
df_dtcr['model'] = 'dtcr'
df_dtcr['cohort'] = 'ISB'
df_dtcr_isb = df_dtcr
DFs.append(df_dtcr_isb)

ds = ['NIH/NIAID','ISB']
df_plot = DFs

fig,ax = plt.subplots(1,1,figsize=(5,5))
ax.plot([0, 1], [0, 1], color='navy', linestyle='--')
for ii,d in enumerate(ds,0):
    score = roc_auc_score(df_plot[ii]['y_test'],df_plot[ii]['y_pred'])
    fpr,tpr,_ = roc_curve(df_plot[ii]['y_test'],df_plot[ii]['y_pred'])
    ax.plot(fpr,tpr,label=d+' (%0.2f)' % score)
    ax.spines['bottom'].set_color('black')
    ax.spines['left'].set_color('black')
    ax.set_facecolor('white')
ax.legend(prop={'size': 15},loc='lower right',facecolor='white')
ax.set_xlabel('False Positive Rate',fontsize=18)
ax.set_ylabel('True Positive Rate',fontsize=18)
plt.tight_layout()
plt.savefig('figures/auc_dtcr.eps')
