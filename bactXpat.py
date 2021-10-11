import matplotlib.pyplot as plt
import numpy as np
import os
import glob
import pandas as pd
import seaborn as sns
import scipy
from sklearn import metrics
import pybedtools
import coolbox
from coolbox.api import *

from joblib.externals.loky import set_loky_pickler
from joblib import parallel_backend
from joblib import Parallel, delayed
from joblib import wrap_non_picklable_objects

traces= glob.glob('data/meta/*bedgraph')
orif=ori='data/OriC_db.bed'

@delayed
@wrap_non_picklable_objects

# for jac in traces:
def oric_calc(jac,orif,out):
    zzz=pd.read_csv(jac,sep='\t',names=['tmp','start','end','ratio'])
    ori=pd.read_csv(orif,sep='\t')
    tmp=ori[ori['Refseq'].str.contains('NC_'+jac.split('/')[2].split('_')[1])]
    ###match bins here### SOMEHOW -- not true false --
    ###just use min "start-end" to round exact start end
    table=pd.DataFrame()
    print('running: '+jac)
    plt.figure(figsize=(20, 5))
    try:
        ss=np.int(np.round(tmp['start'].item(),-1))
        tt=np.int(np.round(tmp['end'].item(),-1))
        plt.bar(data=zzz[(zzz['start']>ss-250000)&(zzz['end']<tt+250000)],x='start',height='ratio',label='align',color='blue')#][zzz['start']=='1665929']
        plt.bar(data=tmp,x='start',height='OriC AT content',label='OriC start',color='orange')
        plt.bar(data=tmp,x='end',height='OriC AT content',label='OriC end',color='orange')

        plt.axes([.2, .6, .6, .2])
        plt.bar(data=zzz[(zzz['start']>ss-100) &(zzz['end']<tt+100)],x='start',height='ratio',label='align',color='blue')#][zzz['start']=='1665929']
        plt.bar(data=tmp,x='start',height='OriC AT content',label='OriC start',color='orange')
        plt.bar(data=tmp,x='end',height='OriC AT content',label='OriC end',color='orange')

        plt.xlabel("OriC ranges from "+str(tmp['start'].item())+" to "+str(tmp['end'].item())+" for "
            +str(tmp['Organism'].item())+" in "+str(os.path.basename(jac).split('_')[2].split('.')[0]),fontsize=14)
        plt.savefig('img/oric/'+str(os.path.basename(jac).split('.')[0])+'.png',dpi=300,bbox_inches = "tight")

        tmpd=zzz[(zzz['start']>ss-100) &(zzz['end']<tt+100)]

        tmp['mean']=np.mean(zzz['ratio'])
        tmp['ori_mean']=np.mean(tmpd['ratio'])

        tmp['median']=np.median(zzz['ratio'])
        tmp['ori_median']=np.median(tmpd['ratio'])
        table.append(tmp)

    except TypeError:
        pass
    except ValueError:
        pass


oric_calc(jac, orif) for jac in (traces))

