#!/bin/bash
#PBS -l nodes=1:ppn=1
#PBS -l mem=48g
#PBS -l walltime=96:00:00
#PBS -N relAb_ratio
#PBS -q cgsd
#PBS -e .err
#PBS -o .out

module purge
module add miniconda3
module add bowtie2
module add samtools
module add bcftools
module add parallel
#module add bamCoverage
source activate mypy3



function relAb_ratio {

# tmp=$HOME'/data/tmp/'
# bacteria=$(ls $tmp/*.bam)
# for bact in $bacteria
# do
bact=$1
# ratio=$1
python run/oric/relAb_ratio.py $bact #$ratio
# done

}

export -f relAb_ratio

parallel relAb_ratio ::: $HOME/data/meta*/*bedgraph  #../../groups/cgsd/gordonq/LauG_Metagenomics_CPOS-200710-CJX-3455a/filtered

# source activate miniconda3; source activate mypy3


# chmod +x run/oric/check_bact.py

# if [ $buffer != 0 ];then
#   python netZooPy/netZooPy/milipeed/benchmark/run_predScore.py -i $bench/red -o $bench/red/test/
  # find "data/MotifPipeline/camb_motif_"$buffer"_QCdelta/" -maxdepth 1 -type f -exec rm -rf {} \;
# else
# python run/oric/check_bact.py # -i $bench -o $bench/test/
