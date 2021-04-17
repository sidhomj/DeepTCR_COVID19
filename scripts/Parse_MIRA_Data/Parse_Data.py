""" This script takes the data from Adaptive ImmuneCode and assembles it into a
dataframe that maps TCRs to their respective subject, cohort, ORF. """

import pandas as pd
import numpy as np
import os
from utils import *

mhc = 'cii'
#Import peptide data from Adaptive
file = '../../Data/ImmuneCODE/ImmuneCODE-MIRA-Release002.1/peptide-detail-'+mhc+'.csv'
data = pd.read_csv(file)
cdr3 = data['TCR BioIdentity'].str.split('+',expand=True)
data['beta_sequences'] = np.array(cdr3[0])

#Process TCRs
data = Process_Seq(data,'beta_sequences')

#Map TCRs to peptides
peptides = data['Amino Acids'].str.split(',',expand=True)
peptide_list = []
beta_sequences_list = []
exp_list = []
for ii,pep in enumerate(peptides.iterrows(),0):
    pep = np.array(pep[1])
    pep = pep[pep!=None]
    peptide_list.append(pep)
    beta_sequences_list.append([data['beta_sequences'].iloc[ii]]*len(pep))
    exp_list.append([data['Experiment'].iloc[ii]]*len(pep))

beta_sequences_list = np.hstack(beta_sequences_list)
peptide_list = np.hstack(peptide_list)
exp_list = np.hstack(exp_list)
data = pd.DataFrame()
data['beta_sequences'] = beta_sequences_list
data['peptide'] = peptide_list
data['Experiment'] = exp_list

#Add meta data
meta = pd.read_csv('../../Data/ImmuneCODE/ImmuneCODE-MIRA-Release002.1/subject-metadata.csv',engine='python')
meta['Subject'] = 's_'+meta['Subject'].astype(str)
sample_dict = dict(zip(meta['Experiment'],meta['Subject']))
data['Subject'] = data['Experiment'].map(sample_dict)
covid_dict = dict(zip(meta['Subject'],meta['Cohort']))
data['Cohort'] = data['Subject'].map(covid_dict)

#Add ORF data
df_orf = pd.read_csv('../../Data/ImmuneCODE/ImmuneCODE-MIRA-Release002.1/minigene-hits.csv')
orfs = np.array(df_orf['Amino Acid'])

orf_list = []
for pep in data['peptide']:
    try:
        orf_list.append(orfs[np.where([pep in x for x in orfs])[0][0]])
    except:
        orf_list.append('ORF1ab')

data['orf'] = orf_list
orf_dict = dict(zip(df_orf['Amino Acid'],df_orf['ORF']))
orf_dict['ORF1ab'] = 'ORF1ab'
data['orf_name'] = data['orf'].map(orf_dict)

#Drop duplicates of beta/peptide pairs in patients
data.drop_duplicates(subset=['beta_sequences','peptide','Subject'],inplace=True)
data.to_csv('../../Data/ImmuneCODE/mira_data_parsed_'+mhc+'.csv',index=False)
