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

from joblib.externals.loky import set_loky_pickler
from joblib import parallel_backend
from joblib import Parallel, delayed
from joblib import wrap_non_picklable_objects

j=int(sys.argv[1])

# for j in [0,500,1000]:

high=pd.read_csv('data/Plectoneme/high_gene.fa',sep='\t',header=None)
zero=pd.read_csv('data/Plectoneme/zero_gene.fa',sep='\t',header=None)
high=high[1::2]
zero=zero[1::2]
high['len'] = [len(i) for i in high[0]]
zero['len'] = [len(i) for i in zero[0]]
high=high[(high['len']>j)]
zero=zero[(zero['len']>j)]
A=np.zeros(5)
B=np.zeros(5)
# B=np.zeros([int(len(zero)),6])

@delayed
@wrap_non_picklable_objects

# for i,Swave in tqdm(enumerate(high[0].tolist())):
def calc_plect(Swave,A):
    # A=np.zeros([int(len(Swave)),6])
    Sequence_angle_exp_smth,Sequence_angle_exp,CurvePhase,BasepairCovariance= PlectonemeCode(Swave)
    A[0]=np.mean(Sequence_angle_exp)
    A[1]=np.median(Sequence_angle_exp)
    A[2]=np.percentile(Sequence_angle_exp,95)
    # cc=np.histogram(Sequence_angle_exp,bins=10)
    # jeff[i,0,0:10]=cc[0]
    # jeff[i,1]=cc[1]
    A[3]=trapz(Sequence_angle_exp, dx=5)
    # A[i,4]=simps(Sequence_angle_exp, dx=5)
    A[4]=len(Sequence_angle_exp)
    return A
    
######-----RUN STUFF-----######
A=Parallel(n_jobs=4) (calc_plect(Swave,A) for i,Swave in enumerate(high[0].tolist()))
np.savetxt("data/Plectoneme/high_gene_plect"+str(j)+".txt", A, delimiter="\t")

B=Parallel(n_jobs=4) (calc_plect(Swave,B) for i,Swave in enumerate(zero[0].tolist()))
np.savetxt("data/Plectoneme/none_gene_plect"+str(j)+".txt", B, delimiter="\t")


# jeff=[#stats.ttest_ind(A[:,0]/A[:,5],B[:,0]/B[:,5]),
# # stats.ttest_ind(A[:,1]/A[:,5],B[:,1]/B[:,5]),
# # stats.ttest_ind(A[:,2]/A[:,5],B[:,2]/B[:,5]),
# stats.ttest_ind(A[:,3]/A[:,4],B[:,3]/B[:,4],alternative='greater',permutations=10)]
# # stats.ttest_ind(A[:,4]/A[:,5],B[:388,4]/B[:388,5])]
# np.savetxt("data/Plectoneme/ttest_gene_plect"+str(j)+".txt", jeff, delimiter="\t")
