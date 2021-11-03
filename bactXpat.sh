#!/bin/bash
#PBS -l nodes=1:ppn=1
#PBS -l mem=24g
#PBS -l walltime=96:00:00
#PBS -N bactXpat
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



function bactXpat {
##ipython cannot access other permissioned directories
#for bact in $bacteria
#do
    # bactdir=$HOME'/data/hot_bact/*/*' #Staphylococcus_aureus*' #Escherichia_coli_O157*'
    # bacteria=$(ls $bactdir/*)

    outimg=$HOME'/img/heat/'
    outdat=$HOME'/data/meta/'
    tmpdir=$HOME'/data/tmp/'
    sampledirA='../../groups/cgsd/gordonq/LauG_Metagenomics_CPOS-200710-CJX-3455a/filtered/'
    sampledirB='../../groups/cgsd/hbyao/LauG_Metagenomics_CPOS-191021-DMS-326a/primary_seq'
    samples=$(ls $sampledirA*;ls $sampledirB*;)
    n=0
    bact=$1
    mkdir $outimg
    mkdir $outdat
    #pat=$2
    # printf $bact" "$pat" \n"
    if [ $((n%2)) -eq 0 ];then
    # ((n&2));then
        bac=$(basename "${bact%.fna}")
        eval "bowtie2-build "$bact" "$tmpdir$bac".fna"

        # baA=$(eval "basename "$bact" |cut -d _ -f1")
        # baB=$(eval "basename "$bact" |  cut -d _ -f2")
        grep "$baA $baB" data/doric10/OriC_db.txt | cut -f1,2,3,8 > data/doric10/tmp_oric.bed
        n=$((n+1))
    else
        n=$((n+1))
        continue
    fi
    # printf $n" is "$bact
    # printf $samples
    for sample in $samples
    do
    # sample=$2

        if [ -e $outdat$bac"_"$pat".bedgraph"];
        then
          continue
        fi

        pat=$(eval "basename "$sample" | cut -d _ -f1")
        printf "crossing "$bac" with "$pat"\n">>progress.txt

        eval "bowtie2 -x "$tmpdir$bac".fna -1 "$sampledir$pat"_clean_1.fastq -2 "$sampledir$pat"_clean_2.fastq -S "$tmpdir$bac"_bowtie.sam" ##align
        eval "samtools view -bS "$tmpdir$bac"_bowtie.sam > "$tmpdir$bac"_bowtie.bam"
        eval "samtools sort "$tmpdir$bac"_bowtie.bam -o "$tmpdir$bac"_"$pat"_bowtie.sorted.bam" ##sort and convert
        eval "samtools index "$tmpdir$bac"_"$pat"_bowtie.sorted.bam"
        eval "bamCoverage -b "$tmpdir$bac"_"$pat"_bowtie.sorted.bam -bs 1 -e -of bedgraph -o "$outdat$bac"_"$pat".bedgraph"
        eval "bamCoverage -b "$tmpdir$bac"_"$pat"_bowtie.sorted.bam -bs 1 -e -o "$outdat$bac"_"$pat".bigwig"
        # eval "cut -f2,3,4 "$pat".bedgraph > "$bac"_"$pat".bg"
        # eval "computeMatrix scale-regions -S "$outdat$bac"_"$pat".bigwig -R "$outdat$bac"_"$pat".bedgraph -a 1000 -b 1000 -o "$tmpdir$bac"_"$pat"matrix.mat.gz" ## -R data/tmp_oric.bed

        # eval "plotHeatmap -m "$tmpdir$bac"_"$pat"matrix.mat.gz --kmeans 4 -o "$outimg$bac"_"$pat".png"


    done
    ##set -- arrVarA
    ##eval "computeMatrix scale-regions -S $arrVarA -R $arrVarB -a 1000 -b 1000 -o matrix.mat.gz"

#done

}

export -f bactXpat

parallel bactXpat ::: $HOME/data/*bact/*/* #../../groups/cgsd/gordonq/LauG_Metagenomics_CPOS-200710-CJX-3455a/filtered/*

source activate miniconda3; source activate mypy3
chmod +x run/oric/bactXpat.py

# if [ $buffer != 0 ];then
#   python netZooPy/netZooPy/milipeed/benchmark/run_predScore.py -i $bench/red -o $bench/red/test/
  # find "data/MotifPipeline/camb_motif_"$buffer"_QCdelta/" -maxdepth 1 -type f -exec rm -rf {} \;
# else
python run/oric/bactXpat.py # -i $bench -o $bench/test/

##eval "computeMatrix scale-regions -S NC_016776R0414-THY.bigwig NC_016776R0415-CWH.bigwig NC_016776R0416-LMY.bigwig NC_016776R0404-CPL.bigwig -R NC_016776R0414-THY.bedgraph NC_016776R0415-CWH.bedgraph NC_016776R0416-LMY.bedgraph NC_016776R0404-CPL.bedgraph -a 1000 -b 1000 -o matrix.mat.gz"

##eval "plotHeatmap -m matrix.mat.gz --zMin -3 --zMax 3 --kmeans 4 -o testAAA.png"

