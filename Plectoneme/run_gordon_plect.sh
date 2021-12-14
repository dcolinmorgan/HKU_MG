#!/bin/bash
#PBS -l nodes=1:ppn=1
#PBS -l mem=24g
#PBS -l walltime=96:00:00
#PBS -N bactXplect
#PBS -q cgsd
#PBS -e .err
#PBS -o .out
##qsub run/oric/Plectoneme/run_gordon_plect.sh

module purge
module add miniconda3
module add parallel
source activate mypy3
# parallel bactXpat ::: $HOME/data/*bact/*/*

chmod +x run/oric/Plectoneme/gordon_fasta/main_setup_species_gtf.py
chmod +x run/oric/Plectoneme/gordon_fasta/main_run_patient.py
chmod +x run/oric/Plectoneme/gordon_fasta/get_tpm_groups.py
chmod +x run/oric/Plectoneme/run_gordon_plect.py

function bactXplect {

bact=$1
python run/oric/Plectoneme/gordon_fasta/main_setup_species_gtf.py $bact 
python run/oric/Plectoneme/gordon_fasta/main_run_patient.py $bact 'X317822438_RNAlater' &
python run/oric/Plectoneme/gordon_fasta/main_run_patient.py $bact 'X316701492_RNAlater' &
python run/oric/Plectoneme/gordon_fasta/main_run_patient.py $bact 'X311245214_RNAlater' 
rm -rf run/oric/Plectoneme/gordon_fasta/$bact/*.fna.gz
python run/oric/Plectoneme/gordon_fasta/get_tpm_groups.py $bact
python run/oric/Plectoneme/run_gordon_plect.py high $bact &
python run/oric/Plectoneme/run_gordon_plect.py zero $bact
}

export -f bactXplect
# shuf -n 10 run/oric/Plectoneme/gordon_fasta/bact_list > run/oric/Plectoneme/gordon_fasta/tmp_bact_list

# parallel -a run/oric/Plectoneme/gordon_fasta/tmp_bact_list bactXplect 

parallel bactXplect  ::: 227318 2479767 39779 297314 2030927 1520817 214856 2049048

# find . -cmin -20 |rm -rf