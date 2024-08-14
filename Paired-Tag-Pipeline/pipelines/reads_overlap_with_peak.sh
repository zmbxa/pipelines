#!/bin/bash
##-------------------------
# @Author: XiongXiong
# @Date: 2023/2/8
##-------------------------

show_usage="Usage: bash reads_overlap_with_peak.sh -b <input.bam> -n <name> -w <input.bw> -r <reference.bw> -p <thread>\n\n
			-b|--input-bam \t input bam file\n
			-n|--name \t description of data, in mm10_brain_H3K27ac format\n
			-w|--input-bw \t signal bigwig file\n
			-r|--reference \t reference bigwig file\n
			-p|--thread \t how many threads to use. Default:8\n
			-h|--help \t\t get help document\n"

GETOPT_ARGS=`getopt -o b:n:w:r:p:h -al input-bam:,name:,input-bw:,reference:,thread:,help -- "$@"`
eval set -- "$GETOPT_ARGS"

while [ -n "$1" ]
do
	case "$1" in
		-b|--input-bam) bamFile=$2; shift 2;;
		-n|--name) name=$2; shift 2;;
		-w|--input-bw) signalBW=$2; shift 2;;
		-r|--reference) refBW=$2; shift 2;;                
		-p|--thread) nthread=$2; shift 2;;
		-h|--help) echo -e $show_usage; exit 0;;
		--) break;;
		*) echo -e $1,$2,$show_usage; break;;
	esac
done

if [[ -z ${bamFile} || -z ${name} ]]; then
	echo -e $show_usage
	exit 0

elif [[ ! -d "/storage/zhangyanxiaoLab/xiongxiong/share/ChIP-seq/${name}" ]]; then
	echo "Sorry, no such peak exits!"
	exit 0

else
	
	if [[ -z ${nthread} ]]; then
		nthread=8
	fi
	
	if [[ -z ${signalBW} ]]; then
		
		samtools view -H ${bamFile} | awk '$1 == "@SQ" {OFS="\t";print $2,$3}' - | sed 's/.N://g' > tmp.size
		nFrag=`samtools view -@ ${nthread} ${bamFile} | wc -l `
		bedtools genomecov -ibam ${bamFile} -bg | awk '{printf "%s\t%.0f\t%.0f\t%.2f\n",$1,$2,$3,$4*1000000/"'$nFrag'"}' | awk '{$4/=1;print}' OFS='\t' > tmp.unsorted.bdg
		sort --parallel ${nthread} -k1,1 -k2,2n tmp.unsorted.bdg > tmp.bdg
		bedGraphToBigWig tmp.bdg tmp.size Experiment.bw
		signalBW='Experiment.bw'
		rm tmp.size tmp.unsorted.bdg tmp.bdg
	fi

	if [[ -z ${refBW} ]]; then
		refBW="/storage/zhangyanxiaoLab/xiongxiong/share/ChIP-seq/${name}/bed/${name}.bw"
	fi

	genome=` echo "${name}" | sed 's/_.*//g' `

	bedtools bamtobed -i ${bamFile} > tmp.bam.bed
	nOri=`cat tmp.bam.bed | wc -l `
	egrep -v '^#' /storage/zhangyanxiaoLab/xiongxiong/share/ChIP-seq/${name}/bed/${name}_peaks.xls | awk 'NR>2' > tmp.peaks.bed
	nExpr=`bedtools intersect -u -a tmp.bam.bed -b tmp.peaks.bed | wc -l `
	pExpr=`echo "scale=1; 100*${nExpr}/${nOri}" | bc `

	samtools view -H ${bamFile} | awk '$1 == "@SQ" {OFS="\t";print $2,$3}' - | sed 's/.N://g' > tmp.size
	bedtools shuffle -i tmp.bam.bed -g tmp.size > tmp.rand.bed
	nCtrl=`bedtools intersect -u -a tmp.rand.bed -b tmp.peaks.bed | wc -l `
	pCtrl=`echo "scale=1; 100*${nCtrl}/${nOri}" | bc `

	echo "${nExpr} (${pExpr}%) reads overlap with peak in experiment sample" | tee reads_overlap_peak.summary
	echo "${nCtrl} (${pCtrl}%) reads overlap with peak in random sample" | tee -a reads_overlap_peak.summary

	computeMatrix reference-point -S ${signalBW} ${refBW} --regionsFileName /storage/zhangyanxiaoLab/xiongxiong/Reference/refBed/${genome}_TSS.bed --skipZeros -a 2000 -b 2000 -bs 20 -p ${nthread} -o tmp.plotFile.gz
	plotProfile --yMin 0 -m tmp.plotFile.gz -out tss_enrichment.pdf --perGroup --refPointLabel TSS --dpi 500

	computeMatrix reference-point -S ${signalBW} ${refBW} --regionsFileName /storage/zhangyanxiaoLab/xiongxiong/share/ChIP-seq/${name}/bed/${name}_summits.bed --skipZeros -a 2000 -b 2000 -bs 20 -p ${nthread} -o tmp.plotFile.gz
	plotProfile --yMin 0 -m tmp.plotFile.gz -out peak_summit_enrichment.pdf --perGroup --refPointLabel Summit --dpi 500

	rm tmp.bam.bed tmp.peaks.bed tmp.size tmp.rand.bed tmp.plotFile.gz
fi
