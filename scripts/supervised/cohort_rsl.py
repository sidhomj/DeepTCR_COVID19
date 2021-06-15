"""Creates Residue Sensitivity Logos for top 25 predictive sequences in ISB & NIH/NIAID data sets for the cohort-trained model."""

import pandas as pd
import numpy as np
from DeepTCR.DeepTCR import DeepTCR_WF
import pickle
import matplotlib.pyplot as plt

with open('cohort_model_preds.pkl','rb') as f:
    DFs_pred = pickle.load(f)

