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

df_metrics = pd.read_csv('../../Data/ImmuneCODE/samples.tsv',sep='\t')
metrics = df_metrics.columns[1:10]
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

ss = StandardScaler()
for m in metrics:
    df_merge[m] = ss.fit_transform(np.array(df_merge[m]).reshape(-1,1))


lr = LogisticRegression()
skf = StratifiedKFold(shuffle=True,n_splits=2)
lb = LabelEncoder()
iterations = 100
DFs = []
for d in ds:
    df_temp = df_merge[df_merge['Dataset'] == d].reset_index(drop=True)
    X = np.array(df_temp[metrics])
    Y = np.array(df_temp['severe'])
    Y = lb.fit_transform(Y)
    y_pred = []
    y_test = []
    sample_list = []
    for it in range(iterations):
        for train_idx, test_idx in skf.split(X,Y):
            y_test.append(Y[test_idx])
            sample_list.append(df_temp['sample_name'][test_idx])
            y_pred_ = []
            for ii,m in enumerate(metrics,0):
                lr.fit(X[train_idx,ii].reshape(-1,1),Y[train_idx])
                pred = lr.predict_proba(X[test_idx,ii].reshape(-1,1))
                y_pred_.append(pred[:,1])
            y_pred.append(np.vstack(y_pred_).T)
    y_pred = np.vstack(y_pred)
    y_test = np.hstack(y_test)
    sample_list = np.hstack(sample_list)
    df_pred = pd.DataFrame()
    df_pred['Samples'] = sample_list
    df_pred['y_test'] = y_test
    df_pred = pd.concat([df_pred,pd.DataFrame(y_pred,columns=metrics)],axis=1)
    DFs.append(df_pred)

agg_dict = {}
agg_dict['y_test'] = 'first'
for m in metrics:
    agg_dict[m] = 'mean'

DFs_gp = []
for _ in DFs:
    DFs_gp.append(_.groupby(['Samples']).agg(agg_dict).reset_index())

df_plot = DFs_gp
# df_plot = df_pred

fig,ax = plt.subplots(1,2,figsize=(10,5))
for ii,d in enumerate(ds,0):
    ax[ii].plot([0, 1], [0, 1], color='navy', linestyle='--')
    for m in metrics:
        score = roc_auc_score(df_plot[ii]['y_test'],df_plot[ii][m])
        fpr,tpr,_ = roc_curve(df_plot[ii]['y_test'],df_plot[ii][m])
        ax[ii].plot(fpr,tpr,label=m + '(%0.2f)' % score)
        ax[ii].set_title(d)
        ax[ii].legend(prop={'size': 8},loc='lower right',facecolor='white')
        ax[ii].set_xlabel('False Positive Rate',fontsize=14)
        ax[ii].set_ylabel('True Positive Rate',fontsize=14)
        ax[ii].spines['bottom'].set_color('black')
        ax[ii].spines['left'].set_color('black')
        ax[ii].set_facecolor('white')
plt.tight_layout()
plt.savefig('figures/auc_logistic_ind.eps')






