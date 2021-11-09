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
tax_id=$1
# taxid=cat tax
pat_id=$2
mkdir -p $HOME/$tax_id/beds
mkdir -p $HOME/$tax_id/matrices
mkdir -p $HOME/$tax_id/heatmaps
# sample=cat $samples
# for sample in `cat $idlist`
# do
# python $HOME/run/oric/gordon_plot.py $taxid $sample
bamCoverage --bam $tax_id/bams/$pat_id.sorted.bam -p 48 -o $HOME/$tax_id/beds/$pat_id.bigwig --binSize 1

bamCoverage --bam $tax_id/bams/$pat_id.sorted.bam -p 48 -of bedgraph -o $HOME/$tax_id/beds/$pat_id.bedgraph --binSize 1


computeMatrix scale-regions -S $HOME/$tax_id/beds/$pat_id.bigwig -R $HOME/$tax_id/beds/$pat_id.bedgraph -a 1000 -b 1000 -o $HOME/$tax_id/matrices/$pat_id.matrixC.mat.gz

plotHeatmap -m $HOME/$tax_id/matrices/$pat_id.matrixC.mat.gz --kmeans 4 -o $HOME/$tax_id/heatmaps/$pat_id.png


# done
}

export -f gordon_plot

parallel gordon_plot ::: $(cat $HOME/idlist) ::: $(cat idlist_SRR)


