#!/bin/bash
#PBS -l nodes=1:ppn=1
#PBS -l mem=16g
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

bactdir='bact/'
bacteria=$(ls $bactdir*/*)



function bactXpat {
##ipython cannot access other permissioned directories
#for bact in $bacteria
#do
    sampledir='../../groups/cgsd/gordonq/LauG_Metagenomics_CPOS-200710-CJX-3455a/filtered/'
    samples=$(ls $sampledir*)
    n=0
    bact=$1
    pat=$2
    # printf $bact" "$pat" \n"
    if [ $((n%2)) -eq 0 ];then
    # ((n&2));then
        bac=$(eval "echo "$bact" | cut -d / -f3| cut -d . -f1")
        eval "bowtie2-build "$bact" "$bac".fna"
        baA=$(eval "echo "$bacteria" | cut -d / -f2| cut -d _ -f1")
        baB=$(eval "echo "$bacteria" | cut -d / -f2| cut -d _ -f2")
        grep "$baA $baB" data/OriC_db.txt | cut -f1,2,3,8 > data/tmp_oric.bed
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
        # printf $sample
        pat=$(eval "echo "$sample" | cut -d / -f8|cut -d _ -f1")
        printf "crossing "$bac" with "$pat"\n">>progress.txt

        eval "bowtie2 -x "$bac".fna -1 "$sampledir$pat"_clean_1.fastq -2 "$sampledir$pat"_clean_2.fastq -S "$bac"_bowtie.sam" ##align
        eval "samtools view -bS "$bac"_bowtie.sam > "$bac"_bowtie.bam"
        eval "samtools sort "$bac"_bowtie.bam -o "$bac"_"$pat"_bowtie.sorted.bam" ##sort and convert
        eval "samtools index "$bac"_"$pat"_bowtie.sorted.bam"
        eval "bamCoverage -b "$bac"_"$pat"_bowtie.sorted.bam -of bedgraph -o "$bac"_"$pat".bedgraph"
        eval "bamCoverage -b "$bac"_"$pat"_bowtie.sorted.bam -o "$bac"_"$pat".bigwig"
        eval "cut -f2,3,4 "$pat".bedgraph > "$bac"_"$pat".bg"
        eval "computeMatrix scale-regions -S "$bac"_"$pat".bigwig -R "$bac"_"$pat".bedgraph -a 1000 -b 1000 -o matrix.mat.gz" ## -R data/tmp_oric.bed
        printf "heatmapping "$baA" "$baB" with "$pat"\n" >> progress.txt
        eval "plotHeatmap -m matrix.mat.gz --kmeans 4 -o img/heat/"$bac"_"$pat".png"

        mv *.bigwig *.bedgraph data/meta
        # mv *.png img
        # rm *.sam *.mat.gz *.bam *.bt2 *.bam.bai
#

#        rm *.e* *.o*
    done
    ##set -- arrVarA
    ##eval "computeMatrix scale-regions -S $arrVarA -R $arrVarB -a 1000 -b 1000 -o matrix.mat.gz"

#done

}

export -f bactXpat

parallel bactXpat ::: bact/*/* #../../groups/cgsd/gordonq/LauG_Metagenomics_CPOS-200710-CJX-3455a/filtered/*

##eval "computeMatrix scale-regions -S NC_016776R0414-THY.bigwig NC_016776R0415-CWH.bigwig NC_016776R0416-LMY.bigwig NC_016776R0404-CPL.bigwig -R NC_016776R0414-THY.bedgraph NC_016776R0415-CWH.bedgraph NC_016776R0416-LMY.bedgraph NC_016776R0404-CPL.bedgraph -a 1000 -b 1000 -o matrix.mat.gz"

##eval "plotHeatmap -m matrix.mat.gz --zMin -3 --zMax 3 --kmeans 4 -o testAAA.png"

