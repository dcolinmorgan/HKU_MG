#Usage python main_setup_species_gtf.py tax_id
import sys, os, re, csv
import pandas as pd
from Bio import SeqIO
import subprocess
import numpy as np

tax_id = sys.argv[1]
#tax_id = 853 # Faecalibacterium prausnitzii
#tax_id = 821 # Phocaeicola vulgatus
#tax_id = 2479767 #Blautia sp. SC05B48
#tax_id = 4932 # Saccharomyces cerevisiae

os.chdir('run/oric/Plectoneme/gordon_fasta/')
dir_path = os.getcwd()
tax_path = dir_path+'/'+str(tax_id)


refseq = pd.read_csv("../../../../../../groups/cgsd/gordonq/database/assembly_summary_refseq.txt", sep = '\t', header=1)

#x = refseq[(refseq['taxid']==int(tax_id))]
x = refseq[(refseq['taxid']==int(tax_id)) & (refseq['refseq_category']=="representative genome")]
if x.shape[0] == 0:
    x = refseq[(refseq['taxid']==int(tax_id)) & ((refseq['refseq_category']=="representative genome") | (refseq['assembly_level']=="Complete Genome"))]
    if x.shape[0] == 1:
        ftp = x['ftp_path'].item()
    else:
        ftp = x['ftp_path'].iloc[0]
else:
    ftp = x['ftp_path'].item()
ident = re.match(r'^.*/(.*$)',ftp).group(1)
fna = ftp+'/'+ident+'_genomic.fna.gz'
# os.chdir('~/run/oric/Plectoneme/gordon_fasta')

if not os.path.exists(tax_path):
    os.mkdir(tax_path)
subprocess.call('wget '+fna+" -O "+tax_path+'/'+ident+'.fna.gz', shell= True)
subprocess.call('gzip -d '+tax_path+'/'+ident+'.fna.gz', shell= True)
cds_fna = ftp+'/'+ident+'_cds_from_genomic.fna.gz'
subprocess.call('wget '+cds_fna+" -O "+tax_path+'/'+ident+'_cds_from_genomic.fna.gz', shell= True)
subprocess.call('gzip -d '+tax_path+'/'+ident+'_cds_from_genomic.fna.gz', shell= True)
gtf = ftp+'/'+ident+'_genomic.gtf.gz'

gtf_file = pd.read_csv(gtf, sep='\t', header = None, skiprows = 5)
gtf_file.drop(gtf_file.tail(1).index,inplace=True)
gtf_file = gtf_file.to_numpy()

#Filter overlaping genes with 500bp distance from TSS
non_overlap = []
for chrom in list(set(gtf_file[:,0])):
    dat = gtf_file[np.in1d(gtf_file[:, 0], chrom)]
    slide = dat[0][4]
    non_overlap.append(dat[0])
    for row in dat[1:]:
        if row[3] > slide:
            non_overlap.append(row)
            slide = row[4]
non_overlap = pd.DataFrame(non_overlap)

#non_overlap = non_overlap[[0,3,4,8,6]]
non_overlap[3] = non_overlap[3].astype(int)
non_overlap[4] = non_overlap[4].astype(int)
non_overlap.to_csv(str(tax_id)+'/TSS_region.gtf',index=False, sep='\t', header = False, quoting=csv.QUOTE_NONE)
