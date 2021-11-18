import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import linregress
import glob,os,sys,scipy,coolbox
from tqdm import tqdm
import pandas as pd
import seaborn as sns
# from scipy import stats
# from sklearn import metrics
# import pybedtools
# from coolbox.api import *
# import networkx as nx
from joblib.externals.loky import set_loky_pickler
from joblib import parallel_backend
from joblib import Parallel, delayed
from joblib import wrap_non_picklable_objects
os.chdir('/home/dcmorgan')
os.getcwd()

table=[]
# color=pd.DataFrame()
color0=pd.DataFrame(columns=['species','subj','ratio','slope'])
jac = sys.argv[1]
# ratio = sys.argv[2]
# traces= np.sort(glob.glob('data/meta*/*bedgraph'))
ori=pd.read_csv('data/doric10/OriC_db.bed',sep='\t')
# cc=ori[ori['Organism']==('Escherichia_coli')]
cc=ori
dd=cc[cc['refseq'].str.contains('NC_')]['refseq']

i=0
# for d,c in (enumerate(traces)):
# indices = [i for i, s in enumerate(traces) if c in s]
#     # if not indices.empty:
# if indices !=[]:
# tmp0=ori[ori.refseq=='NC_'+os.path.basename(traces[0]).split('_')[1]]
# jac=traces[d]
ratio=4
e='NC_'+os.path.basename(jac).split('_')[1]
f=os.path.basename(jac).split('_')[2].split('-')[0]
tmp=ori[ori['refseq'].str.contains(e)]
if tmp.empty==False:
    for uu in range(len(tmp)):

        # for dd,cc in tqdm(enumerate(d)):

        # ax = axes[SLICES[i]]
        zzz=pd.read_csv(jac,sep='\t',names=['tmp','start','end','relAb'])


        zzz['len']=zzz['end']-zzz['start']
        jeff3=np.repeat((zzz.relAb), repeats = zzz['len'])
        jeff3=jeff3.reset_index()
        # df0=jeff3.groupby(jeff3.index // 1000).sum()
        # try:

        ppp=tmp.iloc[uu]['oric_start'].item()
        # ppp=tmp['oric_start'].item()
        jj=jeff3.loc[ppp].values[1]
        cc=jeff3.loc[ppp-500:ppp+500+tmp.iloc[uu]['len']]

        # jj=df0.relAb[int(np.round(ppp/1000,-2))]
        # table.append()
        df1=jeff3.sort_values(by='relAb')
        df2=df1[(df1['relAb']>np.median(df1[df1['relAb']!=0])/ratio)&(df1['relAb']<np.median(df1[df1['relAb']!=0])*ratio)]
        df2=df2.reset_index()
        color0['species']='NC_'+os.path.basename(jac).split('_')[1]
        color0['subj']=os.path.basename(jac).split('_')[2].split('-')[0]
        color0['ratio']=np.mean(cc.relAb)/np.mean(jeff3.relAb),
        color0['slope']=linregress((df2.index), (df2.relAb))[0]
        # table.append(jeff22)
        # if i==0:
            # color0.to_csv('data/relAb_ratio5.txt',sep='\t')
            # i=i+1
        # else:
        color0.to_csv('data/relAb500_ratio'+str(ratio)+'.txt',sep='\t',mode='a',header=False,index=False)

        # print(jj/np.mean(df0.relAb))
        # print(jj)
        # except:
        #     pass