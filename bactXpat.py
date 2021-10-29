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
        zzz['len']=zzz['end']-zzz['start']
        ss=np.int(np.round(tmp['start'].item(),-1))
        tt=np.int(np.round(tmp['end'].item(),-1))
        # plt.bar(data=zzz[(zzz['start']>ss-250000)&(zzz['end']<tt+250000)],x='start',height='ratio',label='align',color='blue')#][zzz['start']=='1665929']
        # plt.bar(data=tmp,x='start',height='OriC AT content',label='OriC start',color='orange')
        # plt.bar(data=tmp,x='end',height='OriC AT content',label='OriC end',color='orange')

        # plt.axes([.2, .6, .6, .2])
        # plt.bar(data=zzz[(zzz['start']>ss-100) &(zzz['end']<tt+100)],x='start',height='ratio',label='align',color='blue')#][zzz['start']=='1665929']
        # plt.bar(data=tmp,x='start',height='OriC AT content',label='OriC start',color='orange')
        # plt.bar(data=tmp,x='end',height='OriC AT content',label='OriC end',color='orange')

        # plt.xlabel("OriC ranges from "+str(tmp['start'].item())+" to "+str(tmp['end'].item())+" for "
        #     +str(tmp['Organism'].item())+" in "+str(os.path.basename(jac).split('_')[2].split('.')[0]),fontsize=14)
        # plt.savefig('img/oric/'+str(os.path.basename(jac).split('.')[0])+'.png',dpi=300,bbox_inches = "tight")

        # tmpd=zzz[(zzz['start']>ss-100) &(zzz['end']<tt+100)]

        # tmp['mean']=np.mean(zzz['ratio'])
        # tmp['ori_mean']=np.mean(tmpd['ratio'])

        # tmp['median']=np.median(zzz['ratio'])
        # tmp['ori_median']=np.median(tmpd['ratio'])
        # table.append(tmp)

        tmpd=zzz[(zzz['start']>ss) &(zzz['end']<tt)]
        tmpdd=zzz[(zzz['start']>ss-250) &(zzz['end']<tt+250)]
        tmpddd=zzz[(zzz['start']>ss-500) &(zzz['end']<tt+500)]

        tmp['mean']=np.sum(zzz['relAb']*zzz['len'])/np.sum(zzz['len'])#np.mean(zzz['score'])
        tmp['ori_mean']=np.sum(tmpd['relAb']*tmpd['len'])/np.sum(tmpd['len'])#np.mean(tmpd['score'])
        tmp['ori_mean_250']=np.sum(tmpdd['relAb']*tmpdd['len'])/np.sum(tmpdd['len'])#np.mean(tmpdd['score'])
        tmp['ori_mean_500']=np.sum(tmpddd['relAb']*tmpddd['len'])/np.sum(tmpddd['len'])#np.mean(tmpddd['score'])

        tmp['trim_mean01']=stats.trim_mean(zzz['relAb'],0.01)
        tmp['ori_trim_mean01']=stats.trim_mean(tmpd['relAb'],0.01)
        tmp['ori_trim_mean01_250']=stats.trim_mean(tmpdd['relAb'],0.01)
        tmp['ori_trim_mean01_500']=stats.trim_mean(tmpddd['relAb'],0.01)

        tmp['trim_mean1']=stats.trim_mean(zzz['relAb'],0.1)
        tmp['ori_trim_mean1']=stats.trim_mean(tmpd['relAb'],0.1)
        tmp['ori_trim_mean1_250']=stats.trim_mean(tmpdd['relAb'],0.1)
        tmp['ori_trim_mean1_500']=stats.trim_mean(tmpddd['relAb'],0.1)

        tmp['trim_mean25']=stats.trim_mean(zzz['relAb'],0.25)
        tmp['ori_trim_mean25']=stats.trim_mean(tmpd['relAb'],0.25)
        tmp['ori_trim_mean25_250']=stats.trim_mean(tmpdd['relAb'],0.25)
        tmp['ori_trim_mean25_500']=stats.trim_mean(tmpddd['relAb'],0.25)

        tmp['pat']=pat


        ccc=np.int(np.round(tmp['start'].item(),0))-50
        ddd=np.int(np.round(tmp['end'].item(),0))+50
        www0=zzz[((zzz['start']<ccc)|(zzz['end']>ddd))]
        # www1=zzz[]
        yyy=zzz[(zzz['start']>ccc)&(zzz['end']<ddd)]
        # print([len(zzz),len(www0),len(yyy)])
        if (np.max(www0.relAb)<np.max(yyy.relAb)):
            # print(str(tmp['Organism'].item()))
            # print([np.max(www0.relAb),np.max(yyy.relAb)])
            plt.figure(figsize=(20,8))
            tmpA=plt.hist(www0.relAb,bins=100,log=True)
            # tmpA=plt.hist(www1.relAb,bins=100,log=True)
            tmpA=plt.hist(yyy.relAb,bins=100,log=True,alpha=0.5)
            # tmpD=(np .max(yyy.relAb))
            tmp['OriC max_relAb']=(np.max(yyy.relAb))
            tmp['background max_relAb']=(np.max(www0.relAb))
            (np.max(www0.relAb))
            plt.xlabel(str(tmp['Organism'].item())+' OriC max_relAb:'+str(tmp['OriC max_relAb'])+', background max_relAb:'+str(tmp['background max_relAb']))
            # pattable=pd.concat([pattable,pat],axis=0)
            table=pd.concat([table,tmp],axis=0)
    except ValueError:
        pass



def run_vRand(dataset,table,n):
    aa=dataset.loc[dataset.index.isin(np.unique(table['subj']))]
    aa['ecoli']=1
    aa=aa.select_dtypes(include=np.number)
    aa=aa.dropna(how='any',axis=1)

    bb=dataset.loc[~dataset.index.isin(np.unique(table['subj']))]
    bb['ecoli']=0
    bb=bb.select_dtypes(include=np.number)
    bb=bb.dropna(how='any',axis=1)

    for i in aa.columns:
        try:
            cc=stats.kruskal(aa[i].dropna(),bb[i].dropna())
            dd=stats.mannwhitneyu(aa[i].dropna(),bb[i].dropna())
            ee=stats.ttest_ind(aa[i].dropna(),bb[i].dropna())
            # ff=stats.wilcoxon(aa[i].dropna(),bb[i].dropna())
            if ((cc.pvalue<0.005) & (dd.pvalue<0.005)&(ee.pvalue<0.005)):
                print (i,cc.pvalue,dd.pvalue,ee.pvalue)

        except:
            pass

    for i in np.arange(n):
        jeff=np.random.permutation(dataset.index)[0:len(np.unique(table['subj']))]
        aa=dataset.loc[dataset.index.isin(jeff)]
        aa=aa.select_dtypes(include=np.number)
        aa=aa.dropna(how='any',axis=1)

        bb=dataset.loc[~dataset.index.isin(jeff)]
        bb=bb.select_dtypes(include=np.number)
        bb=bb.dropna(how='any',axis=1)
        for i in aa.columns:
            try:
                cc=stats.kruskal(aa[i].dropna(),bb[i].dropna())
                dd=stats.mannwhitneyu(aa[i].dropna(),bb[i].dropna())
                ee=stats.ttest_ind(aa[i].dropna(),bb[i].dropna())
                # ff=stats.wilcoxon(aa[i].dropna(),bb[i].dropna())
                if ((cc.pvalue<0.005) & (dd.pvalue<0.005)&(ee.pvalue<0.005)):
                    print (i,cc.pvalue,dd.pvalue,ee.pvalue)

            except:
                pass

######-----RUN STUFF-----######
traces= glob.glob('data/meta/*bedgraph')
orif=ori='data/OriC_db.bed'

oric_calc(jac, orif) for jac in (traces))


primary=pd.read_excel('data/Data Raw - Gut Microbiome Cohort Project Database - 300 Cohort v3.0_280921.xlsx',index_col=0,sheet_name='Primary Data')
diet=pd.read_excel('data/Data Raw - Gut Microbiome Cohort Project Database - 300 Cohort v3.0_280921.xlsx',index_col=0,sheet_name='Diet Data')
blood_stool=pd.read_excel('data/Data Raw - Gut Microbiome Cohort Project Database - 300 Cohort v3.0_280921.xlsx',index_col=0,sheet_name='blood and stool biomarkers')
secondary=pd.read_excel('data/Data Raw - Gut Microbiome Cohort Project Database - 300 Cohort v3.0_280921.xlsx',index_col=0,sheet_name='Secondary Data')
MRI=pd.read_excel('data/Data Raw - Gut Microbiome Cohort Project Database - 300 Cohort v3.0_280921.xlsx',index_col=0,sheet_name='MRI scores')


run_vRand(primary,table,100)
run_vRand(diet,table,100)
run_vRand(blood_stool,table,100)
run_vRand(secondary,table,100)
run_vRand(MRI,table,100)

