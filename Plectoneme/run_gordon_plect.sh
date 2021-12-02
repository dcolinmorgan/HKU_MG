#!/bin/bash
#PBS -l nodes=1:ppn=1
#PBS -l mem=24g
#PBS -l walltime=96:00:00
#PBS -N run_gordon_plect
#PBS -q cgsd
#PBS -e .err
#PBS -o .out

module purge
module add miniconda3
module add parallel
source activate mypy3
# parallel bactXpat ::: $HOME/data/*bact/*/*

chmod +x run/oric/Plectoneme/run_gordon_plect.py

# python run/oric/Plectoneme/run_gordon_plect.py 1000 853 &
# python run/oric/Plectoneme/run_gordon_plect.py 500 853 &
python run/oric/Plectoneme/run_gordon_plect.py 0 853 &

# python run/oric/Plectoneme/run_gordon_plect.py 1000 821 &
# python run/oric/Plectoneme/run_gordon_plect.py 500 821 &
python run/oric/Plectoneme/run_gordon_plect.py 0 821 &

# python run/oric/Plectoneme/run_gordon_plect.py 1000 2479767
# python run/oric/Plectoneme/run_gordon_plect.py 500 2479767 &
python run/oric/Plectoneme/run_gordon_plect.py 0 2479767
