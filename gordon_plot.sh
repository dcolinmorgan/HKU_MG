#PBS -l nodes=1:ppn=48
#PBS -l mem=120gb
#PBS -l walltime=24:00:00
#PBS -N gordon_plot
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

# cd $PBS_O_WORKDIR
cd /groups/cgsd/gordonq/TSS_depth/nanopore_data/
#python gordon_plot.py $taxid


function gordon_plot {
taxid=$1
# taxid=cat tax
sample=$2
# sample=cat $samples
# for sample in `cat $idlist`
# do
python $HOME/run/oric/gordon_plot.py $taxid $sample



# done
}

export -f gordon_plot

parallel gordon_plot ::: 853 $(cat idlist_SRR)