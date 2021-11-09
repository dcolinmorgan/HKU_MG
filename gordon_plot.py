# Usage python gordon_plot.py pat_id_list tax_id
import sys, glob, os
import pandas as pd
import subprocess

os.chdir('/groups/cgsd/gordonq/TSS_depth/nanopore_data/')

tax_id = 853#sys.argv[1] #853 Faecalibacterium prausnitzii
# fna_file = glob.glob(str(tax_id)+"/*.fna")
# fna_file = list(filter(lambda x:'cds_from_genomic' not in x, fna_file))[0]
# if not os.path.exists('$HOME'+tax_id+'/bams'):
#     os.mkdir('$HOME/'+tax_id+'/bams')
# if not os.path.exists('$HOME'+tax_id+'/beds'):
#     os.mkdir('$HOME/'+tax_id+'/beds')
# if not os.path.exists('$HOME'+tax_id+'/matrices'):
#     os.mkdir('$HOME/'+tax_id+'/matrices')

#pat_id = 'X317822438_RNAlater'
pat_id = sys.argv[2]

# bamCoverage --bam '853/bams/'+pat_id+'.nextseq.sorted.bam -p 48 -o $HOME/'+pat_id+'.nextseq.bigwig  --binSize 1

# #Get read depth of TSS regions
# subprocess.call('bamCoverage --bam '+str(tax_id)+'/bams/'+pat_id+'.nextseq.sorted.bam -p 48 -o $HOME/'+pat_id+'.nextseq.bigwig  --binSize 1', shell = True)
# subprocess.call('bamCoverage --bam '+str(tax_id)+'/bams/'+pat_id+'.nextseq.sorted.bam -p 48 -o $HOME/'+str(tax_id)+'/beds/'+pat_id+'.nextseq.bed -of bedgraph --binSize 1', shell = True) ## use this to plot TSS and OriC per bact
# subprocess.call('computeMatrix scale-regions -S $HOME'+str(tax_id)+'/beds/'+pat_id+'.nextseq.bigwig -R $HOME/'+str(tax_id)+'/beds/'+pat_id+'.nextseq.bed -a 1000 -b 1000 -o $HOME'+str(tax_id)+'/beds/'+pat_id+'matrixA.mat.gz', shell=True)
# subprocess.call('plotHeatmap -m $HOME'+str(tax_id)+'/beds/'+pat_id+'matrixA.mat.gz --kmeans 4 -o $HOME/'+str(tax_id)+'/beds/'+pat_id+'.png',shell=True)

# subprocess.call('bamCoverage --bam '+str(tax_id)+'/bams/'+pat_id+'.nanopore.sorted.bam -p 48 -o $HOME/'+str(tax_id)+'/beds/'+pat_id+'.nanopore.bigwig  --binSize 1', shell = True)
# subprocess.call('bamCoverage --bam '+str(tax_id)+'/bams/'+pat_id+'.nanopore.sorted.bam -p 48 -o $HOME/'+str(tax_id)+'/beds/'+pat_id+'.nanopore.bed -of bedgraph --binSize 1', shell = True)
# subprocess.call('computeMatrix scale-regions -S $HOME'+str(tax_id)+'/beds/'+pat_id+'.nanopore.bigwig -R $HOME'+str(tax_id)+'/beds/'+pat_id+'.nanopore.bed -a 1000 -b 1000 -o $HOME/'+str(tax_id)+'/beds/'+pat_id+'matrixB.mat.gz', shell=True)
# subprocess.call('plotHeatmap -m $HOME'+str(tax_id)+'/beds/'+pat_id+'matrixB.mat.gz --kmeans 4 -o $HOME/'+str(tax_id)+'/beds/'+pat_id+'.png',shell=True)
tax_id=853
pat_id=SRR10666963
subprocess.call('bamCoverage --bam '+str(tax_id)+'/bams/'+pat_id+'.sorted.bam -p 48 -o $HOME/'+str(tax_id)+'/beds/'+pat_id+'.bigwig --binSize 1', shell = True)
subprocess.call('bamCoverage --bam '+str(tax_id)+'/bams/'+pat_id+'.sorted.bam -p 48 -o $HOME/'+str(tax_id)+'/beds/'+pat_id+'.bed -of bedgraph --binSize 1', shell = True)
subprocess.call('computeMatrix scale-regions -S $HOME'+str(tax_id)+'/beds/'+pat_id+'.bigwig -R $HOME/'+str(tax_id)+'/beds/'+pat_id+'.bed -a 1000 -b 1000 -o $HOME'+str(tax_id)+'/beds/'+pat_id+'matrixC.mat.gz', shell=True)
subprocess.call('plotHeatmap -m $HOME'+str(tax_id)+'/beds/'+pat_id+'matrixC.mat.gz --kmeans 4 -o $HOME/'+str(tax_id)+'/beds/'+pat_id+'.png',shell=True)


bamCoverage --bam 853/bams/SRR10666963.sorted.bam -p 48 -o $HOME/853/beds/SRR10666963.bigwig --binSize 1

