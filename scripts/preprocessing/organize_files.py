import pandas as pd
import numpy as np
import os
import glob
import shutil
from multiprocessing import Pool

df = pd.read_csv('../../Data/ImmuneCODE/ImmuneCODE-Repertoire-Tags-002.2.tsv',sep='\t')
df['sample_name'] = df['sample_name']+'.tsv'
df['Dataset'] = df['Dataset'].str.replace('/','-')

#create folders
folder_write = os.path.join('..','..','Data','ImmuneCODE','repertoires_org')
if not os.path.exists(folder_write):
    os.mkdir(folder_write)

folders = np.unique(df['Dataset'])
for f in folders:
    folder_write = os.path.join('..','..','Data','ImmuneCODE','repertoires_org',f)
    if not os.path.exists(folder_write):
        os.mkdir(folder_write)

def copy_file_(file,f):
    try:
        if not os.path.exists(os.path.join('..', '..', 'Data', 'ImmuneCODE', 'repertoires_org', f, file)):
            shutil.copy(os.path.join('..', '..', 'Data', 'ImmuneCODE', 'repertoires', file),
                        os.path.join('..', '..', 'Data', 'ImmuneCODE', 'repertoires_org', f))
        return None
    except:
        return file

p = Pool(80)
#copy files
for f in folders:
    files = np.array(df[df['Dataset'] == f]['sample_name'])
    out = p.starmap(copy_file_,list(zip(files,[f]*len(files))))
p.join()
p.close()

