""" convert file format so na > np.nan to be compatible with DeepTCR"""

import glob
import os
import pandas as pd
import numpy as np
from multiprocessing import Pool

files = '../../Data/ImmuneCODE/repertoires_raw/*.tsv'
files = glob.glob(files)

def reform(file):
    try:
        if not os.path.exists(os.path.join('../../Data/ImmuneCODE/repertoires',os.path.basename(file))):
            df = pd.read_csv(file,sep='\t')
            df['amino_acid'][df['amino_acid']=='na'] = np.nan
            df.to_csv(os.path.join('../../Data/ImmuneCODE/repertoires',os.path.basename(file)),sep='\t',index=False)
        return None
    except:
        return file

p = Pool(80)
out = p.map(reform,files)
p.join()
p.close()

file = '../../Data/ImmuneCODE/repertoires_raw/INCOV067-AC-3_TCRB.tsv'
df = pd.read_csv(file, sep='\t')
df_out = pd.DataFrame(columns=range(0,52))
df_out[3] = df[df.columns[1]]
df_out[4] = df[df.columns[2]]
df_out[49] = df[df.columns[7]]
df_out[35] = df[df.columns[14]]
df_out[42] = df[df.columns[21]]
df_out.to_csv(os.path.join('../../Data/ImmuneCODE/repertoires', os.path.basename(file)), sep='\t', index=False)

out = [x for x in out if x != None]
df = pd.read_csv(out[0], sep='\t')
