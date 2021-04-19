from Bio.Alphabet import IUPAC
from DeepTCR.DeepTCR import DeepTCR_WF, DeepTCR_SS
from scipy.stats import spearmanr, rankdata
import numpy as np

def Process_Seq(df,col):
    #Drop null values
    df = df.dropna(subset=[col])

    #strip any white space and remove non-IUPAC characters
    df[col] = df[col].str.strip()
    df = df[~df[col].str.contains(r'[^A-Z]')]
    iupac_c = set((list(IUPAC.IUPACProtein.letters)))
    all_c = set(''.join(list(df[col])))
    searchfor = list(all_c.difference(iupac_c))
    if len(searchfor) != 0:
        df = df[~df[col].str.contains('|'.join(searchfor))]
    return df

def run_models(beta_sequences,df,models=None):
    DTCR = DeepTCR_WF('../supervised/niaid_model', device=2)
    pred = DTCR.Sample_Inference(beta_sequences=beta_sequences, batch_size=100000, models=models)
    df['niaid_pred'] = pred[:, 1]
    df['niaid_rank'] = rankdata(pred[:, 1])

    DTCR = DeepTCR_WF('../supervised/isb_model', device=2)
    pred = DTCR.Sample_Inference(beta_sequences=beta_sequences, batch_size=100000, models=models)
    df['isb_pred'] = pred[:, 1]
    df['isb_rank'] = rankdata(pred[:, 1])

    # DTCR = DeepTCR_WF('../supervised/combined_model', device=2)
    # pred = DTCR.Sample_Inference(beta_sequences=beta_sequences, batch_size=100000, models=models)
    # df['combo_pred'] = pred[:, 1]
    # df['combo_rank'] = rankdata(pred[:, 1])

    df['count'] = 1

def group_df(df,var_group,categories,thresh=1):
    agg_dict = {}
    for ii, c in enumerate(categories, 0):
        if ii < var_group:
            agg_dict[c] = 'first'
    agg_dict.update({'niaid_pred': 'mean', 'isb_pred': 'mean', #'combo_pred': 'mean',
                     'niaid_rank': 'mean', 'isb_rank': 'mean', #'combo_rank': 'mean',
                     'count': 'sum', })
    df_gp = df.groupby(categories[var_group]).agg(agg_dict).reset_index()
    df_gp = df_gp[df_gp['count'] >= thresh]
    return df_gp