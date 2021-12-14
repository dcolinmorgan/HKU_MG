#Usage python main_run_patient.py tax_id pat_id
import sys, glob, os
import pandas as pd
import subprocess
os.chdir('run/oric/Plectoneme/gordon_fasta/')
tax_id = sys.argv[1] #853 Faecalibacterium prausnitzii
fna_file = glob.glob(str(tax_id)+"/*.fna")
fna_file = list(filter(lambda x:'cds_from_genomic' not in x, fna_file))[0]
#if not os.path.exists(tax_id+'/bams'):
#    os.mkdir(tax_id+'/bams')
#if not os.path.exists(tax_id+'/beds'):
#    os.mkdir(tax_id+'/beds')
#if not os.path.exists(tax_id+'/matrices'):
#    os.mkdir(tax_id+'/matrices')
if not os.path.exists(tax_id+'/kallisto_res'):
    os.mkdir(tax_id+'/kallisto_res')

#pat_id = 'X317822438_RNAlater'
# X316701492_RNAlater
# X311245214_RNAlater
pat_id = sys.argv[2]
# Align raw reads to reference genome
# subprocess.call('bwa index '+fna_file, shell=True)
# subprocess.call('bwa mem -t 48 '+fna_file+' DNA_primary/'+pat_id+'_1.fastq DNA_primary/'+pat_id+'_2.fastq |samtools view --threads 48 -b - |samtools sort - -o '+str(tax_id)+'/bams/'+pat_id+'.sorted.bam --threads 48', shell=True)
# subprocess.call('samtools index '+str(tax_id)+'/bams/'+pat_id+'.sorted.bam -@ 48', shell = True)
# subprocess.call('picard MarkDuplicates --REMOVE_DUPLICATES true --INPUT '+str(tax_id)+'/bams/'+pat_id+'.sorted.bam --OUTPUT '+str(tax_id)+'/bams/'+pat_id+'.rm_dup.sorted.bam --METRICS_FILE '+str(tax_id)+'/bams/'+pat_id+'.metrics.txt', shell = True)

# # Get read depth of TSS regions
# subprocess.call('bamCoverage --bam '+str(tax_id)+'/bams/'+pat_id+'.sorted.bam -p 48 --extendReads -o '+str(tax_id)+'/beds/'+pat_id+'.bed -of bedgraph --binSize 1', shell = True)

#Quantify RNA data
cds_fna = glob.glob(str(tax_id)+"/*genomic.fna")[0]
subprocess.call('kallisto index '+cds_fna+' -i '+str(tax_id)+'/'+str(tax_id)+'_kallisto.db', shell = True)

subprocess.call('kallisto quant -t 48 -i '+str(tax_id)+'/'+str(tax_id)+'_kallisto.db -o '+str(tax_id)+'/kallisto_res/'+pat_id+' ../../../../../../../groups/cgsd/gordonq/TSS_depth/matched_public_data/RNA_primary/'+pat_id+'_1.fastq ../../../../../../../groups/cgsd/gordonq/TSS_depth/matched_public_data/RNA_primary/'+pat_id+'_2.fastq', shell = True)

                
# kallisto quant -t 48 -i '+str(tax_id)+'/'+str(tax_id)+'_kallisto.db -o '+str(tax_id)+'/kallisto_res/'+pat_id+' groups/cgsd/gordonq/TSS_depth/matched_public_data/RNA_primary/'+pat_id+'_1.fastq groups/cgsd/gordonq/TSS_depth/matched_public_data/RNA_primary/'+pat_id+'_2.fastq