#!/bin/sh

#####################################################################################################
# Author: Marek Cmero
#
# Script for generating 3, 4 and 5 cluster mixtures
# used in the SVclone paper. Note that patient 001 
# WGS bam files for bM and gM are required (Hong 2015).
#
# required software: samtools
#####################################################################################################

sample1=001bM
sample2=001gM

bam1=SM_${sample1}.recal.bam
bam2=SM_${sample2}.recal.bam

# tumour and sequence stats
cov1=51.5
cov2=58.9
pur1=0.49
pur2=0.46
pl1=2.26
pl2=2.22

# run vars
threads=8
mem=4

#####################################################################################################
# 3 cluster mixtures
#####################################################################################################

for i in $(seq 1 9); do
    p1=0.${i}
    p2=`bc -l <<< "1 - $p1"`
    p2=0${p2}

    base=${sample1}_p${p1}_${sample2}_p${p2}_merged
    merged_bam=${base}.bam

    echo "making 2 sample mix: ${sample1} prop:${p1}, ${sample2} p:${p2}"

    echo "subsampling $sample1"
    samtools view -@ $threads -s $p1 -b $bam1 > ${sample1}_p${p1}.bam
    samtools sort -@ $threads -m ${mem}G ${sample1}_p${p1}.bam -o ${sample1}_p${p1}_sorted.bam
    samtools index -@ $threads ${sample1}_p${p1}_sorted.bam

    echo "subsampling $sample2"
    samtools view -@ $threads -s $p2 -b $bam2 > ${sample2}_p${p2}.bam
    samtools sort -@ $threads -m ${mem}G ${sample2}_p${p2}.bam -o ${sample2}_p${p2}_sorted.bam
    samtools index -@ $threads ${sample2}_p${p2}_sorted.bam

    echo "merging..."
    samtools merge -n -@ $threads $merged_bam ${sample1}_p${p1}_sorted.bam ${sample2}_p${p2}_sorted.bam
    samtools sort -@ $threads -m ${mem}G $merged_bam -o ${base}_sorted.bam
    samtools index -@ $threads ${base}_sorted.bam
done

#####################################################################################################
# high-cluster mixtures preprocessing
#####################################################################################################

echo "preparing samples for high-cluster mixtures"
echo "splitting samples into odd and even chromosomes..."

samtools view -@ $threads -hb $bam1 1 3 5 7 9 11 13 15 17 19 21 > ${sample1}_odd_chroms.bam
samtools index -@ $threads ${sample1}_odd_chroms.bam

samtools view -@ $threads -hb $bam1 2 4 6 8 10 12 14 16 18 20 22 > ${sample1}_even_chroms.bam
samtools index -@ $threads ${sample1}_even_chroms.bam

samtools view -@ $threads -hb $bam2 1 3 5 7 9 11 13 15 17 19 21 > ${sample2}_odd_chroms.bam
samtools index -@ $threads ${sample2}_odd_chroms.bam

samtools view -@ $threads -hb $bam2 2 4 6 8 10 12 14 16 18 20 22 > ${sample2}_even_chroms.bam
samtools index -@ $threads ${sample2}_even_chroms.bam

echo "downsampling $sample2 chrom bams to normalise coverage for 5 cluster mixture..."

p1=`bc -l <<< "($cov1 / $pl1) * $pur1"`
p2=`bc -l <<< "($cov2 / $pl2) * $pur2"`
p=`bc -l <<< "$p1 / $p2"`

samtools view -@ $threads -s 4895${p} -b ${sample2}_odd_chroms.bam > ${sample2}_odd_chroms_downs.bam
samtools view -@ $threads -s 3987${p} -b ${sample2}_even_chroms.bam > ${sample2}_even_chroms_downs.bam

samtools index -@ $threads ${sample2}_odd_chroms_downs.bam
samtools index -@ $threads ${sample2}_even_chroms_downs.bam

#####################################################################################################
# 4-cluster mixture
#####################################################################################################

sample1=001bM
sample2=001gM
sample3=001bM
prop1=20
prop2=40
prop3=60
type1=odd
type2=odd
type3=even

base=p${prop1}_${prop2}_${prop3}_${sample1}_${sample2}_${sample3}_merged
merged_bam=${base}.bam
bam1=${sample1}_${type1}_chroms.bam
bam2=${sample2}_${type2}_chroms.bam
bam3=${sample3}_${type3}_chroms.bam

p1=`bc -l <<< "$prop1 / 100" | sed '/\./ s/\.\{0,1\}0\{1,\}$//'`
p2=`bc -l <<< "$prop2 / 100" | sed '/\./ s/\.\{0,1\}0\{1,\}$//'`
p3=`bc -l <<< "$prop3 / 100" | sed '/\./ s/\.\{0,1\}0\{1,\}$//'`

echo "making 4 sample mixture:"
echo "${sample1} $type1 chroms:${prop1}%"
echo "${sample2} $type2 chroms:${prop2}%"
echo "${sample3} $type3 chroms:${prop3}%"

echo "subsampling $p1 from $bam1..."
samtools view -@ $threads -s 3489${p1} -b $bam1 > ${sample1}_p${prop1}_${type1}_chroms.bam
samtools sort -@ $threads -m ${mem}G ${sample1}_p${prop1}_${type1}_chroms.bam \
    -o ${sample1}_p${prop1}_${type1}_chroms_sorted.bam
samtools index -@ $threads ${sample1}_p${prop1}_${type1}_chroms_sorted.bam
result1=${sample1}_p${prop1}_${type1}_chroms_sorted.bam

echo "subsampling $p2 from $bam2..."
samtools view -@ $threads -s 8973${p2} -b $bam2 > ${sample2}_p${prop2}_${type2}_chroms.bam
samtools sort -@ $threads -m ${mem}G ${sample2}_p${prop2}_${type2}_chroms.bam \
    -o ${sample2}_p${prop2}_${type2}_chroms_sorted.bam
samtools index -@ $threads ${sample2}_p${prop2}_${type2}_chroms_sorted.bam
result2=${sample2}_p${prop2}_${type2}_chroms_sorted.bam

echo "subsampling $p3 from $bam3..."
samtools view -@ $threads -s 7252${p3} -b $bam3 > ${sample3}_p${prop3}_${type3}_chroms.bam
samtools sort -@ $threads -m ${mem}G ${sample3}_p${prop3}_${type3}_chroms.bam \
    -o ${sample3}_p${prop3}_${type3}_chroms_sorted.bam
samtools index -@ $threads ${sample3}_p${prop3}_${type3}_chroms_sorted.bam
result3=${sample3}_p${prop3}_${type3}_chroms_sorted.bam

echo "merging results..."
samtools merge -n -@ $threads $merged_bam $result1 $result2 $result3
samtools sort -@ $threads -m ${mem}G $merged_bam -o ${base}_sorted.bam
samtools index -@ $threads ${base}_sorted.bam

#####################################################################################################
# 5-cluster mixture
#####################################################################################################

sample1=001bM
sample2=001bM
sample3=001gM
sample4=001gM
prop1=80
prop2=60
prop3=20
prop4=40

echo "making 5 sample mixture:"
echo "${sample1} odd chroms:${prop1}%"
echo "${sample2} even chroms:${prop2}%"
echo "${sample3} odd chroms: ${prop3}%"
echo "${sample4} even chroms: ${prop4}%"

base=p${prop1}_${prop2}_${prop3}_${prop4}_${sample1}_${sample2}_${sample3}_${sample4}.merged
merged_bam=${base}.bam
bam1=${sample1}_odd_chroms.bam
bam2=${sample2}_even_chroms.bam
bam3=${sample3}_odd_chroms_downs.bam
bam4=${sample4}_even_chroms_downs.bam

p1=`bc -l <<< "$prop1 / 100" | sed '/\./ s/\.\{0,1\}0\{1,\}$//'`
p2=`bc -l <<< "$prop2 / 100" | sed '/\./ s/\.\{0,1\}0\{1,\}$//'`
p3=`bc -l <<< "$prop3 / 100" | sed '/\./ s/\.\{0,1\}0\{1,\}$//'`
p4=`bc -l <<< "$prop4 / 100" | sed '/\./ s/\.\{0,1\}0\{1,\}$//'`

echo "subsampling $p1 from $bam1"
samtools view -@ $threads -s 7836${p1} -b $bam1 > ${sample1}_p${prop1}_odd_chroms.bam
samtools sort -@ $threads -m ${mem}G ${sample1}_p${prop1}_odd_chroms.bam \
    -o ${sample1}_p${prop1}_odd_chroms_sorted.bam
samtools index -@ $threads ${sample1}_p${prop1}_odd_chroms_sorted.bam
result1=${sample1}_p${prop1}_odd_chroms_sorted.bam

echo "subsampling $p2 from $bam2"
samtools view -@ $threads -s 0983${p2} -b $bam2 > ${sample2}_p${prop2}_even_chroms.bam
samtools sort -@ $threads -m ${mem}G ${sample2}_p${prop2}_even_chroms.bam \
    -o ${sample2}_p${prop2}_even_chroms_sorted.bam
samtools index -@ $threads ${sample2}_p${prop2}_even_chroms_sorted.bam
result2=${sample2}_p${prop2}_even_chroms_sorted.bam

echo "subsampling $p3 from $bam3"
samtools view -@ $threads -s 3746${p3} -b $bam3 > ${sample3}_p${prop3}_odd_chroms.bam
samtools sort -@ $threads -m ${mem}G ${sample3}_p${prop3}_odd_chroms.bam \
    -o ${sample3}_p${prop3}_odd_chroms_sorted.bam
samtools index -@ $threads ${sample3}_p${prop3}_odd_chroms_sorted.bam
result3=${sample3}_p${prop3}_odd_chroms_sorted.bam

echo "subsampling $p4 from $bam3"
samtools view -@ $threads -s 9353${p4} -b $bam4 > ${sample4}_p${prop4}_even_chroms.bam
samtools sort -@ $threads -m ${mem}G ${sample4}_p${prop4}_even_chroms.bam \
    -o ${sample4}_p${prop4}_even_chroms_sorted.bam
samtools index -@ $threads ${sample4}_p${prop4}_even_chroms_sorted.bam
result4=${sample4}_p${prop4}_even_chroms_sorted.bam

echo "merging results..."
samtools merge -n -@ $threads $merged_bam $result1 $result2 $result3 $result4
samtools sort -@ $threads -m ${mem}G $merged_bam -o ${base}_sorted.bam
samtools index -@ $threads ${base}_sorted.bam
