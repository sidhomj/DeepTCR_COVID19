"""Creates figures to visualize deep learning model performance of cohort model."""
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
import pickle

with open('cohort_model_ep_preds.pkl','rb') as f:
    DFs_pred = pickle.load(f)

for df in DFs_pred:
    DFs_pred[df] = DFs_pred[df].groupby(['Samples']).agg({'y_pred': 'mean', 'y_test': 'mean'}).reset_index()

ds = ['COVID-19-NIH/NIAID','COVID-19-ISB']
fig,ax = plt.subplots(1,1,figsize=(5,5))
ax.plot([0, 1], [0, 1], color='navy', linestyle='--')
for ii,d in enumerate(ds,0):
    score = roc_auc_score(DFs_pred[d]['y_test'],DFs_pred[d]['y_pred'])
    fpr,tpr,_ = roc_curve(DFs_pred[d]['y_test'],DFs_pred[d]['y_pred'])
    ax.plot(fpr,tpr,label=d[9:]+' (AUC = %0.2f)' % score)
    ax.spines['bottom'].set_color('black')
    ax.spines['left'].set_color('black')
    ax.set_facecolor('white')
ax.legend(prop={'size': 15},loc='lower right',facecolor='white')
ax.set_xlabel('False Positive Rate',fontsize=18)
ax.set_ylabel('True Positive Rate',fontsize=18)
plt.tight_layout()
plt.savefig('figures_ant/auc_cohort_ep.eps')

