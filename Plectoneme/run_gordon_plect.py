# from scipy.integrate import simps
from numpy import trapz
import pandas as pd
import numpy as np
from tqdm import tqdm
from scipy import stats
import sys,glob,os
sys.path.insert(1, './run/oric/Plectoneme/')
import PlectonemeCode
# importlib.reload(sys.modules['PlectonemeCode'])
from PlectonemeCode import PlectonemeCode
from Bio.Seq import Seq
from joblib.externals.loky import set_loky_pickler
from joblib import parallel_backend
from joblib import Parallel, delayed
from joblib import wrap_non_picklable_objects

j=str(sys.argv[1])
h=str(int(sys.argv[2]))

# for j in [0,500,1000]:

IN=pd.read_csv('run/oric/Plectoneme/gordon_fasta/'+h+'/'+h+'_'+j+'_gene.fa',sep='\t',header=None)
# zero=pd.read_csv('data/Plectoneme/'+h+'_none_gene.fa',sep='\t',header=None)
high=IN[1::2]
sign=IN[0::2]
# if high[1::2]
    # 

# zero=zero[1::2]
high['len'] = [len(i) for i in high[0]]
# zero['len'] = [len(i) for i in zero[0]]
# high=high[(high['len']>j)]
# zero=zero[(zero['len']>j)]
A=np.zeros(8)
# B=np.zeros(5)
# B=np.zeros([int(len(zero)),6])

@delayed
@wrap_non_picklable_objects

# for i,Swave in tqdm(enumerate(high[0].tolist())):
def calc_plect(Swave,A,gene,spec,sign):
    # if (sign.item()).find('+')==0:
        # Swave=Seq.reverse_complement(Swave)
    # A=np.zeros(8)
    Sequence_angle_exp_smth,Sequence_angle_exp= PlectonemeCode(Swave)
    A[0]=np.mean(Sequence_angle_exp)
    A[1]=np.median(Sequence_angle_exp)
    A[2]=np.sum(Sequence_angle_exp>=1)
    A[3]=trapz(Sequence_angle_exp, dx=5)
    A[4]=Sequence_angle_exp.argmax()
    A[5]=len(Sequence_angle_exp)
    A[6]=len(Sequence_angle_exp.nonzero())
    A[7]=np.max(Sequence_angle_exp)
    try:
        with open('data/Plectoneme/'+gene+'_plect_'+spec+'.txt','ab') as f:
            np.savetxt(f,[Sequence_angle_exp.nonzero()], fmt='%5d', delimiter='\t')
    except:
        with open('data/Plectoneme/'+gene+'_plect_'+spec+'.txt','ab') as f:
            np.savetxt(f,[np.nan], delimiter='\t')
    # np.savetxt('data/Plectoneme/'+gene+'_summ_plect_'+str(spec)+'.txt', A, delimiter="\t")
    # return A
    
######-----RUN STUFF-----######
Parallel(n_jobs=24) (calc_plect(Swave,A,j,h,sign[i]) for i,Swave in enumerate(high[0].tolist()))

#####--------PLOT--------######
traces=glob.glob("data/Plectoneme/*_plect_*.txt")
AA=pd.DataFrame()
BB=pd.DataFrame()
from statannot import add_stat_annotation
for trace in traces:
    typ=os.path.basename(trace).split('_')[0]
    spec=os.path.basename(trace).split('_')[2].split('.')[0]
    if typ=='high':
        A=pd.read_csv(trace, delimiter="\t",names=['mean','median','>=1','area','max_ind','len','plect_len','max'])
        # A=A[A['area']!=0]
        A['type']=typ
        A['species']=spec
        # A['filt']=str(os.path.basename(trace).split('_')[4].split('.')[0])
        AA=AA.append(A)
    else:
        B=pd.read_csv(trace, delimiter="\t",names=['mean','median','>=1','area','max_ind','len','plect_len','max'])
        # B=B[B['area']!=0]
        B['type']=typ
        B['species']=spec
        # B['filt']=str(os.path.basename(trace).split('_')[4].split('.')[0])
        BB=BB.append(B)
    
C=AA.append(BB)
C=C[C['area']!=0]
# C=C[C['len']>=1000]

# C=C[C['max_ind']!=0]
C['density']=C['area']/C['len']
C['plect_%']=C['>=1']/C['len']
C['plect_loc_%']=C['max_ind']/C['len']
# C['TSS_plect']=(C['max_ind']<=550)*(C['max_ind']>=10)+0
C['TSS>1']=((C['max_ind']<=10)+0)*C['max']
C.dropna(inplace=True)
C=C.reset_index()
ax=sns.violinplot(data=C,y='TSS>1',x='species',hue='type',cut=0,scale='count')#,kind="violin",cut=0,scale="count", inner="quartile")

add_stat_annotation(ax, data=C, x='species', y='TSS>1',hue='type',
                    box_pairs=[(("821", "high"), ("821", "none")),
                                 (("853", "high"), ("853", "none")),
                                 (("2479767", "high"), ("2479767", "none"))
                                ],
                    test='t-test_ind', text_format='star', loc='inside', verbose=2)
plt.legend(loc='upper left', bbox_to_anchor=(1.03, 1))