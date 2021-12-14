#Usage: python get_tpm_groups.py taxid
import pandas as pd
import os, sys, glob, subprocess, re
import numpy as np
os.chdir('run/oric/Plectoneme/gordon_fasta/')
tax_id = sys.argv[1]
# tax_id = 853

def get_tpm(file):
    dat = pd.read_table(file)
    dat = dat[dat['target_id'].str.contains("P_")]
    dat.target_id = [re.match(r'.*_(.P_[^_]*)',x).group(1) for x in dat.target_id]
    dat = dat.drop_duplicates(subset='target_id')
    dat = dat[["target_id","tpm"]]
    return(dat)
    
#Convert Gene ID to Protein ID
cds_file = [name for name in os.listdir(str(tax_id))  if 'genomic.fna' in name][0]
with open(str(tax_id)+'/'+cds_file,'r') as f:
    file = f.readlines()
genes = [line.strip() for line in file if line[0] ==">"]
genes = [g for g in genes if "protein_id" in g]
protein_2_gene = {re.match(r".*\[protein_id=(\S*)\]",line).group(1):re.match(r".*\[locus_tag=(\S*)\]",line).group(1) for line in genes}

#Read transcript abundance
transcript_files = glob.glob(str(tax_id)+"/kallisto_res/*/abundance.tsv")
samp1 = get_tpm(transcript_files[0])
samp2 = get_tpm(transcript_files[1])
samp3 = get_tpm(transcript_files[2])
merged = pd.merge(samp1.set_index('target_id'), samp2.set_index('target_id'), left_index=True, right_index=True)
merged = pd.merge(merged, samp3.set_index('target_id'), left_index=True, right_index=True)
merged.index = merged.index.map(protein_2_gene)
zero = list(merged.index[merged.sum(axis=1)==0])
high = list(merged.index[merged.min(axis=1)>=100])

### Get Sequences ###
fna = glob.glob(str(tax_id)+'/*.fna')
fna = [x for x in fna if 'cds' not in x][0]
gtf = pd.read_table(str(tax_id)+'/TSS_region.gtf', header = None)
gtf = gtf[gtf[2]=='gene']
gtf["ID"] = [re.match(r'.*; locus_tag "([^"]*)', x).group(1) for x in gtf[8]]
def bingenes(gtf,typ,typp):
    zero_genes = gtf[gtf["ID"].isin(typ)]
    end_index=len(zero_genes)-1
    zero_genes=zero_genes.reset_index(drop=True)
    zero_genes.loc[1:,'up']=(zero_genes.loc[1:,3]-50)
    zero_genes.loc[0,'up']=(zero_genes.loc[0,3])
    zero_genes.loc[:end_index-1,'down']=(zero_genes.loc[:end_index-1,4]+50)
    zero_genes.loc[end_index,'down']=(zero_genes.loc[end_index,4])
    zero_genes['down']=zero_genes['down'].astype('int')
    zero_genes['up']=zero_genes['up'].astype('int')
    zero_genes[10]=1 ## required in 5th column for some reason

    zero_genes[[0,'up','down',"ID",10,6]].to_csv(str(tax_id)+"_temp.bed",sep='\t', header = False, index = False)
    subprocess.call("bedtools getfasta -fi "+fna+" -bed "+str(tax_id)+"_temp.bed -s -name > "+str(tax_id)+"/"+str(tax_id)+"_"+str(typp)+"_gene.fa; rm "+str(tax_id)+"_temp.bed", shell = True)
    
bingenes(gtf,high,'high')
bingenes(gtf,zero,'zero')