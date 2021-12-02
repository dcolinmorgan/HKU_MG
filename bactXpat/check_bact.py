import matplotlib.pyplot as plt
import numpy as np
import os
import glob
import pandas as pd
import seaborn as sns
import scipy
from scipy import stats
from sklearn import metrics
import pybedtools
import coolbox
from coolbox.api import *
import networkx as nx
from joblib.externals.loky import set_loky_pickler
from joblib import parallel_backend
from joblib import Parallel, delayed
from joblib import wrap_non_picklable_objects
os.chdir('/home/dcmorgan')
os.getcwd()

for i in (pd.unique(bact_stat['bact'])):
    for j in (pd.unique(bact_stat['type'])):

        aa=bact_stat[(bact_stat['type']==j)&(bact_stat['bact']==i)]
        print('unqiue: '+str(len(pd.unique(aa['patient']))/len(aa))+', total: '+str(len(aa)))
        tmp=plt.hist(np.sort(aa['q0']),bins=100,log=True)
        tmp=plt.xticks(rotation=90)