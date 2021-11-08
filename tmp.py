
#Usage python generate_matrix.py TSS_regions genome_depth outfolder
import os
import sys
import pandas as pd
import numpy as np
import re
from tqdm import tqdm


TSS_regions = pd.read_table(sys.argv[1], names = ['genome','upstream','downstream','id','strand'])
#TSS_regions = pd.read_table("374606/TSS_region_500.bed", names = ['genome','upstream','downstream','id','strand'])
TSS_regions.id = [re.match(r'.*; locus_tag "([^"]*)', x).group(1) for x in TSS_regions.id]
#Remove genes without 500bp sequence around TSS
TSS_regions = TSS_regions[abs(TSS_regions.upstream - TSS_regions.downstream) == 1000]
bp_depth = pd.read_table(sys.argv[2], names = ['genome','start','end','depth'])
#bp_depth = pd.read_table("374606/beds/Aminobacter_msh1.nextseq.bed", names = ['genome','start','end','depth'])

TSS_regions = TSS_regions.to_numpy()
bp_depth_mat = bp_depth[['start','end','depth']].to_numpy()
bp_depth_genome = bp_depth.genome

mat = []
if len(list(set(TSS_regions[:,0])))==1:
    for row in tqdm(TSS_regions):
        left = np.min([int(row[1]),int(row[2])])
        right  without 500bp sequence around TSS
TSS_regions = TSS_regions[abs(TSS_regions.upstream - TSS_regions.downstream) == 1000]
bp_depth = pd.read_table(sys.argv[2], names = ['genome','start','end','depth'])
#bp_depth = pd.read_table("374606/beds/Aminobacter_msh1.nextseq.bed", names = ['genome','start','end','depth'])
TSS_regions = TSS_regions.to_numpy()
bp_depth_mat = bp_depth[['start','end','depth']].to_numpy()
bp_depth_genome = bp_depth.genome
mat = []
if len(list(set(TSS_regions[:,0])))==1:
    for row in tqdm(TSS_regions):
        left = np.min([int(row[1]),int(row[2])])
        right = np.max([int(row[1]),int(row[2])])
        window = []
        for i in range(left,right+1):
            ind = np.count_nonzero(bp_depth_mat[:,1]<i)
            depth = bp_depth_mat[ind][2]
            window.append(depth)
        if row[4]=="-": window.reverse()
        mat.append(window)
else:
    for row in tqdm(TSS_regions):
        chrom = row[0]
        left = np.min([int(row[1]),int(row[2])])
        right = np.max([int(row[1]),int(row[2])])
        window = []
        depth_dat = bp_depth[bp_depth['genome']==chrom][['start','end','depth']].to_numpy()
        last = depth_dat[-1][1]
        for i in range(left,right+1):
            if i > last:
                depth = 0
            else:
                ind = np.count_nonzero(depth_dat[:,1]<i)
                depth = depth_dat[ind][2]
            window.append(depth)
        if row[4]=="-": window.reverse()
        mat.append(window)
gene_names = TSS_regions[:,3]
mat = np.array(mat)
x = np.hstack((gene_names[np.newaxis].T,mat))
#np.savetxt("c.aero/bamCov_matrix_non_overlap.csv",x,fmt='%5s',delimiter=',')
np.savetxt(sys.argv[3],x,fmt='%5s',delimiter=',') np.max([int(row[1]),int(row[2])])
        window = []
        for i in range(left,right+1):
            ind = np.count_nonzero(bp_depth_mat[:,1]<i)
            depth = bp_depth_mat[ind][2]
            window.append(depth)
        if row[4]=="-": window.reverse()
        mat.append(window)
else:
    for row in tqdm(TSS_regions):
        chrom = row[0]
        left = np.min([int(row[1]),int(row[2])])
        right = np.max([int(row[1]),int(row[2])])
        window = []
        depth_dat = bp_depth[bp_depth['genome']==chrom][['start','end','depth']].to_numpy()
        last = depth_dat[-1][1]
        for i in range(left,right+1):
            if i > last:
                depth = 0
            else:
                ind = np.count_nonzero(depth_dat[:,1]<i)
                depth = depth_dat[ind][2]
            window.append(depth)
        if row[4]=="-": window.reverse()
        mat.append(window)

gene_names = TSS_regions[:,3]

mat = np.array(mat)
x = np.hstack((gene_names[np.newaxis].T,mat))
#np.savetxt("c.aero/bamCov_matrix_non_overlap.csv",x,fmt='%5s',delimiter=',')
np.savetxt(sys.argv[3],x,fmt='%5s',delimiter=',')