traces=glob.glob("data/Plectoneme/*_gene_plect_*.txt")
AA=pd.DataFrame()
BB=pd.DataFrame()
from statannot import add_stat_annotation
for trace in traces:
    type=os.path.basename(trace).split('_')[0]
    if type=='high':
        A=pd.read_csv(trace, delimiter="\t",names=['mean','median','95_max','area','len'])
        # A=A[A['area']!=0]
        A['type']=type
        A['species']=os.path.basename(trace).split('_')[3]
        A['filt']=str(os.path.basename(trace).split('_')[4].split('.')[0])
        AA=AA.append(A)
    else:
        B=pd.read_csv(trace, delimiter="\t",names=['mean','median','95_max','area','len'])
        # B=B[B['area']!=0]
        B['type']=type
        B['species']=os.path.basename(trace).split('_')[3]
        B['filt']=str(os.path.basename(trace).split('_')[4].split('.')[0])
        BB=BB.append(B)
    
C=AA.append(BB)
C=C[C['area']!=0]
C=C[C['filt']=='1000']
C['density']=C['area']/C['len']
C['95_max/len']=C['95_max']/C['len']
C.dropna(inplace=True)
C=C.reset_index()
ax=sns.violinplot(data=C,y='95_max',x='species',hue='type',cut=0,scale='count')#,kind="violin",cut=0,scale="count", inner="quartile")

add_stat_annotation(ax, data=C, x='species', y='95_max',hue='type',
                    box_pairs=[(("821", "high"), ("821", "none")),
                                 (("853", "high"), ("853", "none")),
                                 (("2479767", "high"), ("2479767", "none"))
                                ],
                    test='t-test_ind', text_format='full', loc='inside', verbose=2)
plt.legend(loc='upper left', bbox_to_anchor=(1.03, 1))