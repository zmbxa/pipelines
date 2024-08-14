#!/bin/bash
##-------------------------
# @Author: XiongXiong
# @Date: 2023/2/7
##-------------------------

[ $# != 2 ] && { echo "Usage: bash downsample.sh <path_to_DNA> <path_to_RNA>"
		 echo "<path_to_DNA> path to DNA bam file, usually in prefix_Assign2BC.bam format"
		 echo "<path_to_RNA> path to RNA bam file, usually in prefix_Assign2BC.bam format"
		 exit 1
		}

mkdir -p ./downsample
cd ./downsample

bam_dna=${1}
bam_rna=${2}

prefix_dna=`echo ${bam_dna} | sed 's/\/.*\///g; s/_Assign2BC//g; s/\.bam//g' `
prefix_rna=`echo ${bam_rna} | sed 's/\/.*\///g; s/_Assign2BC//g; s/\.bam//g' `

samtools sort -@ 40 -m 4G ${bam_dna} > tmp.dna.sorted.bam
samtools sort -@ 40 -m 4G ${bam_rna} > tmp.rna.sorted.bam

echo -e "Count\tDup\tLib\tName" > Downsample_Dup.txt
echo -e "0\t0\tDNA\t${prefix_dna}" >> Downsample_Dup.txt
echo -e "0\t0\tRNA\t${prefix_rna}" >> Downsample_Dup.txt

for prop in `seq -f "%.2f" 0.05 0.05 1 `; do
	picard DownsampleSam \
		I=tmp.dna.sorted.bam \
		O=tmp.downsampled.bam \
		P=${prop} \
		R=100 \
		ACCURACY=0.00001 \
		STRATEGY=ConstantMemory  

	bedtools bamtobed -i tmp.downsampled.bam | awk '{split($4,a,":"); inds=length(a)} {print $1, $2, a[inds], a[inds-3]":"a[inds-2]":"a[inds-1], $6, $4}' OFS='\t' > tmp.unsorted.txt
	sort --parallel 40 -k1,1 -k2,2n -k3,3 -k4,4 -k5,5 -k6,6 tmp.unsorted.txt | bedtools groupby -i - -g 1,2,3,4,5 -c 6 -o collapse | sed 's/,.*//g' | cut -f 6 > tmp.unsorted2.txt
	sort --parallel 40 -k1,1 tmp.unsorted2.txt > tmp.Read1.txt

	nSample=`cat tmp.unsorted.txt | wc -l `

	nLine=`cat tmp.Read1.txt | wc -l `
	nDup=`echo "scale=1; ${nSample}-${nLine}" | bc`
	pDup=`echo "scale=1; 100*${nDup}/${nSample}" | bc`
	echo -e "${nSample}\t${pDup}\tDNA\t${prefix_dna}" >> Downsample_Dup.txt	

	picard DownsampleSam \
		I=tmp.rna.sorted.bam \
		O=tmp.downsampled.bam \
		P=${prop} \
		R=100 \
		ACCURACY=0.00001 \
		STRATEGY=ConstantMemory  

	bedtools bamtobed -i tmp.downsampled.bam | awk '{split($4,a,":"); inds=length(a)} {print $1, $2, a[inds], a[inds-3]":"a[inds-2]":"a[inds-1], $6, $4}' OFS='\t' > tmp.unsorted.txt
	sort --parallel 40 -k1,1 -k2,2n -k3,3 -k4,4 -k5,5 -k6,6 tmp.unsorted.txt | bedtools groupby -i - -g 1,2,3,4,5 -c 6 -o collapse | sed 's/,.*//g' | cut -f 6 > tmp.unsorted2.txt
	sort --parallel 40 -k1,1 tmp.unsorted2.txt > tmp.Read1.txt

	nSample=`cat tmp.unsorted.txt | wc -l `

	nLine=`cat tmp.Read1.txt | wc -l `
	nDup=`echo "scale=1; ${nSample}-${nLine}" | bc`
	pDup=`echo "scale=1; 100*${nDup}/${nSample}" | bc`
	echo -e "${nSample}\t${pDup}\tRNA\t${prefix_rna}" >> Downsample_Dup.txt	

done
Rscript /storage/zhangyanxiaoLab/xiongxiong/scripts/paired_tag/analysis/line_Dup.R
rm tmp.dna.sorted.bam tmp.rna.sorted.bam tmp.downsampled.bam tmp.unsorted.txt tmp.unsorted2.txt tmp.Read1.txt
echo "Run Competed."
cd ..
