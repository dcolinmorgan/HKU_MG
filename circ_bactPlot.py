import matplotlib.pyplot as plt
import numpy as np
import os
import glob
import pandas as pd
import random
# The plot consists of 2 rows and 5 columns (10 categories in total)
import random
# The plot consists of 2 rows and 5 columns (10 categories in total)
fig, axes = plt.subplots(8, 6, figsize=(20, 20), subplot_kw={"projection": "polar"})
fig.patch.set_facecolor("white")
PLUS = 10
table=[]

traces= np.sort(glob.glob('data/meta*/*bedgraph'))
ori=pd.read_csv('data/doric10/OriC_db.bed',sep='\t')
# cc=ori[ori['Organism']!=('Bifidobacterium_longum')]
cc=ori
dd=cc[cc['refseq'].str.contains('NC_')]['refseq']


ori=pd.read_csv('data/doric10/OriC_db.bed',sep='\t')
SLICES = [(i, j) for i in range(8) for j in range(6)]
i=0
sns.palplot(sns.color_palette("Set2", 40))
color=pd.DataFrame()
color0=pd.DataFrame(columns=['color','species'])
# for jac in traces[0::151]:

for d in dd:
    indices = [i for i, s in enumerate(traces) if d in s]
    # if not indices.empty:
    if indices !=[]:
        jac=traces[indices[0]]
        ax = axes[SLICES[i]]
        zzz=pd.read_csv(jac,sep='\t',names=['tmp','start','end','relAb'])
        tmp=ori[ori['refseq'].str.contains('NC_'+jac.split('/')[2].split('_')[1])]

        cc=tmp['Organism'].iloc[0]
        if cc in table:
            colorA=color0[color0['species']==cc]['color']
        else:
            color['color']=["#"+''.join([random.choice('0123456789ABCDEF') for j in range(6)])] #list(np.random.choice(range(256), size=3))
            color['species']=(cc)
            color0=color0.append(color)
            colorA=color0[color0['species']==cc]['color']
            table.append(cc)
        i=i+1
        # aa=zzz[['start','relAb']]
        # aa.rename(columns={"start":"Name","relAb":"Value"},inplace=True)
        # df0=zzz.groupby(zzz.index // 2000).sum()/100

        # TIMES=zzz['relAb'][1::100]
        zzz['len']=zzz['end']-zzz['start']
        jeff3=np.repeat((zzz.relAb), repeats = zzz['len'])
        jeff3=jeff3.reset_index()
        df0=jeff3.groupby(jeff3.index // 1000).sum()/100
        # print(np.max(jeff3.groupby(jeff3.index // 10000).sum()['relAb']))
        index=jeff3.index[0::1000].tolist()

        TIMES=df0['relAb']
        TIME_MAX = np.max(TIMES)
        TIME_MIN = np.min(TIMES)
        # TIMES=(TIMES-TIME_MIN)/(TIME_MAX-TIME_MIN)
        ANGLES = np.linspace(0, 2 * np.pi, len(TIMES), endpoint=False)
        HEIGHTS = np.array(TIMES)

        ax.set_rscale('symlog')
        ax.set_theta_offset(np.pi / 2)
        ax.set_theta_direction(-1)

        ax.vlines(ANGLES, 0 + PLUS, HEIGHTS + PLUS, lw=.9,color=colorA)
        # ax.scatter(ANGLES, HEIGHTS + PLUS*1000, s=(HEIGHTS),color=colorA);

        ax.spines["start"].set_color("none")
        ax.spines["polar"].set_color("none")

        ax.grid(False)
        ax.set_xticks([])
        ax.set_yticklabels([])
        ax.set_rmax(ax.get_rmax() * 2)
        ppp=tmp['oric_start'].iloc[0]
        ax.text(
            x=ppp, y=int(np.round(ppp/10,-4)),s='OriC_loc_at:'+str(ppp),
            ha="right", va="center", ma="center", size=8,
             weight="bold")

        ax.text(
            x=0.5, y=1, s=str(tmp['Organism'].iloc[0]),
            va="center", ha="center", ma="center",
            fontsize=8, fontweight="bold", linespacing=0.87, transform=ax.transAxes)

        # ax.text(
        #     x=0.5, y=0.46, s='OriC_loc_at:'+str(tmp['end'].iloc[0]),
        #      va="center", ha="center",  ma="center",
        #     fontsize=8, linespacing=0.87, transform=ax.transAxes)
    else:
        pass

# fig.subplots_adjust(wspace=0, hspace=0)
fig.savefig("img/circle_bact_unfilt.png",dpi=300,bbox_inches = "tight")
