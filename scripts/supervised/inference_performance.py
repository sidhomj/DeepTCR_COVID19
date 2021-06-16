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

files = ['model_isb_data_niaid.csv','model_niaid_data_isb.csv']
names = ['ISB_on_NIH/NIAID','NIH/NIAID_on_ISB']

fig,ax = plt.subplots(1,1,figsize=(5,5))
ax.plot([0, 1], [0, 1], color='navy', linestyle='--')
for ii,(f,n) in enumerate(zip(files,names),0):
    df_preds = pd.read_csv(f)
    if n == 'ISB_on_NIH/NIAID':
        df_preds['y_test'] = df_preds['label']
    else:
        df_preds['y_test'] = df_preds['label']=='severe'

    score = roc_auc_score(df_preds['y_test'],df_preds['Pred'])
    fpr,tpr,_ = roc_curve(df_preds['y_test'],df_preds['Pred'])
    ax.plot(fpr,tpr,label=n+' (AUC = %0.2f)' % score)
    ax.spines['bottom'].set_color('black')
    ax.spines['left'].set_color('black')
    ax.set_facecolor('white')
ax.legend(prop={'size': 12},loc='lower right',facecolor='white')
ax.set_xlabel('False Positive Rate',fontsize=18)
ax.set_ylabel('True Positive Rate',fontsize=18)
plt.tight_layout()
plt.savefig('figures/cross_inference.eps')