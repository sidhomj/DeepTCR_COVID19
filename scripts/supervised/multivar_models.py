"""Compares deep learning models to logistic regression models fit on TCR metrics of diversity and abundance."""
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

#isb
df_dtcr = pd.read_csv('isb_preds.csv')
df_dtcr = df_dtcr.groupby(['Samples']).agg({'y_pred':'mean','y_test':'mean'}).reset_index()
df_dtcr['model'] = 'dtcr'
df_dtcr['cohort'] = 'ISB'
df_dtcr_isb = df_dtcr

#niaid
df_dtcr = pd.read_csv('niaid_preds.csv')
df_dtcr = df_dtcr.groupby(['Samples']).agg({'y_pred':'mean','y_test':'mean'}).reset_index()
df_dtcr['model'] = 'dtcr'
df_dtcr['cohort'] = 'NIH/NIAID'
df_dtcr_niaid = df_dtcr

df_dtcr = pd.concat([df_dtcr_isb,df_dtcr_niaid])


df_log = pd.read_csv('../conventional_metrics/isb_preds_log.csv')
df_log['model'] = 'log'
df_log['cohort'] = 'ISB'
df_log_isb = df_log

df_log = pd.read_csv('../conventional_metrics/niaid_preds_log.csv')
df_log['model'] = 'log'
df_log['cohort'] = 'NIH/NIAID'
df_log_niaid = df_log

df_log = pd.concat([df_log_isb,df_log_niaid])

df = pd.merge(df_dtcr,df_log,on=['Samples','y_test','cohort'])
df = df[['Samples','y_test','cohort','model_x','y_pred_x','model_y','y_pred_y']]

ds = ['NIH/NIAID','ISB']
lr = LogisticRegression()
skf = StratifiedKFold(shuffle=True,n_splits=2)
iterations = 100
DFs = []
DFs_coef = []
for d in ds:
    df_temp = df[df['cohort'] == d].reset_index(drop=True)
    X = np.array(df_temp[['y_pred_x','y_pred_y']])
    Y = np.array(df_temp['y_test'])
    y_pred = []
    y_test = []
    sample_list = []
    coef = []
    for it in range(iterations):
        for train_idx, test_idx in skf.split(X, Y):
            lr.fit(X[train_idx], Y[train_idx])
            pred = lr.predict_proba(X[test_idx])
            y_pred.append(pred[:, 1])
            y_test.append(Y[test_idx])
            sample_list.append(df['Samples'][test_idx])
            coef.append(lr.coef_)

    y_pred = np.hstack(y_pred)
    y_test = np.hstack(y_test)
    sample_list = np.hstack(sample_list)
    coef = np.vstack(coef)

    df_pred = pd.DataFrame()
    df_pred['Samples'] = sample_list
    df_pred['y_test'] = y_test
    df_pred['y_pred'] = y_pred
    DFs.append(df_pred)

    df_coef = pd.DataFrame(coef)
    df_coef.columns = ['dtcr','log']
    DFs_coef.append(df_coef)

agg_dict = {}
agg_dict['y_test'] = 'first'
agg_dict['y_pred'] = 'mean'

DFs_gp = []
for _ in DFs:
    DFs_gp.append(_.groupby(['Samples']).agg(agg_dict).reset_index())

df_plot = DFs_gp

#check if 99 CI crosses 0
print(np.percentile(DFs_coef[1],0.5,axis=0))
print(np.percentile(DFs_coef[1],99.5,axis=0))


fig,ax = plt.subplots(1,1,figsize=(5,5))
ax.plot([0, 1], [0, 1], color='navy', linestyle='--')
for ii,d in enumerate(ds,0):
    score = roc_auc_score(df_plot[ii]['y_test'],df_plot[ii]['y_pred'])
    fpr,tpr,_ = roc_curve(df_plot[ii]['y_test'],df_plot[ii]['y_pred'])
    ax.plot(fpr,tpr,label=d+' (AUC = %0.2f)' % score)
    ax.spines['bottom'].set_color('black')
    ax.spines['left'].set_color('black')
    ax.set_facecolor('white')
ax.legend(prop={'size': 15},loc='lower right',facecolor='white')
ax.set_xlabel('False Positive Rate',fontsize=18)
ax.set_ylabel('True Positive Rate',fontsize=18)
plt.tight_layout()
plt.savefig('figures/auc_logistic_all.eps')

temp = []
for ii,df in enumerate(DFs_coef,0):
    df_temp = pd.melt(df,value_vars=['dtcr','log'])
    df_temp['cohort'] = ds[ii]
    temp.append(df_temp)

df_coef = pd.concat(temp)
sns.violinplot(data=df_coef,x='cohort',y='value',hue='variable')
plt.gca().set_facecolor('white')
plt.gca().spines['bottom'].set_color('black')
plt.gca().spines['left'].set_color('black')
plt.axhline(linestyle='--',color='k')
plt.legend(facecolor='white')
plt.xlabel('')
plt.xticks(fontsize=18)
plt.ylabel('Model Coefficients',fontsize=18)
plt.gca().legend(prop={'size': 15},loc='upper right',facecolor='white')
plt.tight_layout()
plt.savefig('figures/coeff.eps')



