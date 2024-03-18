#!/bin/bash
##-------------------------
# @Author: XiongXiong
# @Date: 2023/1/19
##-------------------------

show_usage="Usage: bash bam2mtx.sh -n <prefix> -g <reference-genome> -i <input.bam> -l <library-type> -o <output_directory> -p <thread> -s <False>\n
			e.g. bash bam2mtx.sh -n merged_DNA -g mm10 -i ./DNA.bam -l DNA\n\n
			-n|--prefix \t\t prefix of output\n
			-g|--reference-genome \t reference genome, must be one of hg38, ce11 and mm10\n
			-i|--input \t\t path to DNA bam\n
			-l|--library-type \t library type, must be one of DNA and RNA\n
			-o|--output \t\t output directory, Default ./matrix\n
			-p|--thread \t how many threads to use. Default:8\n
			-s|--sub \t True or False. Generate sub-library matrix, i.e., use bb:cc:dd barcode instead of aa:bb:cc:dd. Default:False\n
			-h|--help \t\t get help document\n"

GETOPT_ARGS=`getopt -o n:g:i:l:o:p:s:h -al prefix:,reference-genome:,input:,library-type:,output:,thread:,sub:,help -- "$@"`
eval set -- "$GETOPT_ARGS"

while [ -n "$1" ]
do
	case "$1" in
		-n|--prefix) prefix=$2; shift 2;;
		-g|--reference-genome) reference_genome=$2; shift 2;;
		-i|--input) input=$2; shift 2;;
		-l|--library-type) library=$2; shift 2;;
		-o|--output) output=$2; shift 2;;
		-p|--thread) nthread=$2; shift 2;;
		-s|--sub) subMtx=$2; shift 2;;
		-h|--help) echo -e $show_usage; exit 0;;
				--) break;;
		*) echo -e $1,$2,$show_usage; break;;
        esac
done

if [[ -z ${nthread} ]]; then
	nthread=8
fi

if [[ -z ${prefix} || -z ${reference_genome} || -z ${input} || -z ${library} ]]; then
	echo -e $show_usage
	exit 0
elif [[ ${library} != "DNA" && ${library} != "RNA" ]]; then
	echo -e "Unknown data type, should be one of DNA or RNA.\n"
	echo -e $show_usage
	exit 0
elif [[ ! -f ${input} ]]; then
	echo -e "Can not find bam file.\n"
	echo -e $show_usage
	exit 0
else
	if [[ -z ${output} ]]; then
		output=`pwd`"/matrix"
	fi

	if [[ -z ${subMtx} || ${subMtx} = "False" ]]; then
		bedtools bamtobed -i ${input} | awk '{print $1, $2, $3, substr($4,length($4)-21,11), $4}' OFS='\t' > tmp.UsefulReads.bed
	elif [[ ${subMtx} = "True" ]]; then
		bedtools bamtobed -i ${input} | awk '{print $1, $2, $3, substr($4,length($4)-18,8), $4}' OFS='\t' > tmp.UsefulReads.bed
	else
		echo -e "Invalid option:--sub.\n"
		echo -e $show_usage
		exit 0
	fi

	sort --parallel ${nthread} -k1,1 -k2,2n tmp.UsefulReads.bed > tmp.UsefulReads.Sorted.bed
	case ${library} in
		"DNA")
			awk '{print $1, $2, $3, $1"-"$2"-"$3, "bin"NR}' OFS='\t' /storage/zhangyanxiaoLab/niuyuxiao/annotations/refBed/binGenome/mm10_gfpCre_bin_5k.bed | bedtools intersect -sorted -F 0.5 -wa -wb -a - -b tmp.UsefulReads.Sorted.bed | cut -f 4,5,9,10 > tmp.Count.txt
			sort --parallel ${nthread} tmp.Count.txt | bedtools groupby -i - -g 1,2,3 -c 4 -o count > ${prefix}_CountMatrix.txt
			;;
		"RNA")
			featureCounts -T ${nthread} -a /storage/zhangyanxiaoLab/niuyuxiao/annotations/gtf/mm10_gfp_Cre.gtf -t gene -f -g gene_id --extraAttributes gene_name -M -o tmp.Counts ${input}
			awk 'NR>2' tmp.Counts | awk '{gsub(/\.[0-9]*/, "", $1); print $2, $3-1, $4, $1, $8, $5, $7}' OFS='\t' | awk '$5>0' > tmp.ExpressedGenes.bed
			sort --parallel ${nthread} -k1,1 -k2,2n tmp.ExpressedGenes.bed > ${prefix}_ExpressedGenes.bed
			bedtools intersect -sorted -F 0.5 -wa -wb -a ${prefix}_ExpressedGenes.bed -b tmp.UsefulReads.Sorted.bed | awk '{print $4, $7, $11, substr($12,length($12)-9,10)}' OFS='\t' > tmp.Count.txt
			sort -u --parallel ${nthread} tmp.Count.txt | bedtools groupby -i - -g 1,2,3 -c 4 -o count > ${prefix}_CountMatrix.txt
			;;
	esac

	python /storage/zhangyanxiaoLab/xiongxiong/scripts/paired_tag/scripts/bamToMatrix.py ${prefix}_CountMatrix.txt

	sort --parallel ${nthread} -u -k1,1n tmp.Barcode.txt | cut -f 2 > barcodes.tsv
	nBC=`cat barcodes.tsv | wc -l `

	sort --parallel ${nthread} -u -k1,1n tmp.Genes.txt | cut -f 2,3 > features.tsv
	nFea=`cat features.tsv | wc -l `

	nCount=`cat tmp.Matrix.txt | wc -l `

	echo -e "%%MatrixMarket matrix coordinate real general\n%\n${nFea} ${nBC} ${nCount}" | cat - tmp.Matrix.txt > matrix.mtx
	mkdir -p ${output}
	pigz -p ${nthread} *.tsv matrix.mtx
	mv *.tsv.gz matrix.mtx.gz ${output}
	rm -f tmp.UsefulReads.bed tmp.UsefulReads.Sorted.bed tmp.Counts tmp.Counts.summary tmp.ExpressedGenes.bed tmp.Count.txt tmp.Barcode.txt tmp.Genes.txt tmp.Matrix.txt
fi
