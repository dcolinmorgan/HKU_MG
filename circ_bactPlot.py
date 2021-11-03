import matplotlib.pyplot as plt
import numpy as np
import os
import glob
import pandas as pd
import random
# The plot consists of 2 rows and 5 columns (10 categories in total)
fig, axes = plt.subplots(5, 5, figsize=(18, 12), subplot_kw={"projection": "polar"})
fig.patch.set_facecolor("white")
PLUS = 1000
table=[]
# Define the slices used to iterate through 'axes'.
# It iterates in a rowwise manner.
# It starts in the first row, and iterates over all the columns of that row
# from left to right, then it goes to the next row and does the same.
# SLICES = [(i, j) for i in range(2) for j in range(5)]
traces= np.sort(glob.glob('data/meta_unfilt/*bedgraph'))
ori=pd.read_csv('data/doric10/OriC_db.bed',sep='\t')
SLICES = [(i, j) for i in range(5) for j in range(5)]
i=0
sns.palplot(sns.color_palette("Set2", 40))
color=pd.DataFrame()
color0=pd.DataFrame(columns=['color','species'])
for jac in traces[0::151]:
    ax = axes[SLICES[i]]

    # jac=traces[1]
    zzz=pd.read_csv(jac,sep='\t',names=['tmp','start','end','relAb'])
    tmp=ori[ori['Refseq'].str.contains('NC_'+jac.split('/')[2].split('_')[1])]
    if not tmp.empty:
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
        TIMES=zzz['relAb'][1::2000]
        TIME_MAX = np.max(TIMES)
        TIME_MIN = np.min(TIMES)
        ANGLES = np.linspace(0, 2 * np.pi, len(TIMES), endpoint=False)
        HEIGHTS = np.array(TIMES)

        ax.set_rscale('symlog')
        ax.set_theta_offset(np.pi / 2)
        ax.set_theta_direction(-1)

        ax.vlines(ANGLES, 0 + PLUS, HEIGHTS + PLUS, lw=.9,color=colorA)
        ax.scatter(ANGLES, HEIGHTS + PLUS, s=(HEIGHTS),color=colorA);

        ax.spines["start"].set_color("none")
        ax.spines["polar"].set_color("none")

        ax.grid(False)
        ax.set_xticks([])
        ax.set_yticklabels([])
        ax.set_rmax(ax.get_rmax() * 2)

        ax.text(
            x=tmp['end'].iloc[0], y=10,s=str(tmp['Refseq'].iloc[0]),
            ha="right", va="center", ma="center", size=8,
             weight="bold")

        ax.text(
            x=0.5, y=0.58, s=str(tmp['Organism'].iloc[0]),
            va="center", ha="center", ma="center",
            fontsize=8, fontweight="bold", linespacing=0.87, transform=ax.transAxes)

        ax.text(
            x=0.5, y=0.46, s='OriC_loc_at:'+str(tmp['end'].iloc[0]),
             va="center", ha="center",  ma="center",
            fontsize=8, linespacing=0.87, transform=ax.transAxes)
    else:
        pass

fig.subplots_adjust(wspace=0, hspace=0)