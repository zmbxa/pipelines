#!/bin/bash
##-------------------------
# @Author: XiongXiong
# @Date: 2023/1/19
##-------------------------

show_usage="Usage: bash pseudo_bulk.sh -i <input_bam> -t <cell-cluster> -o <output> -p <thread>\n
			e.g. bash pseudo_bulk.sh -i DNA.bam -t batch1_file_list.txt\n\n
			-i|--input-bam \t input bam file\n\n
			-t|--cell-cluster-file \t 2-column file in <brcode><TAB><cell-cluster> format\n\n
			-o|--output \t [optional]output directory\n\n
			-p|--thread \t how many threads to use. Default:8\n
			-h|--help \t\t get help document\n"

GETOPT_ARGS=`getopt -o i:t:p:h -al input:,cell-cluster-file:,thread:,help -- "$@"`
eval set -- "$GETOPT_ARGS"

while [ -n "$1" ]
do
	case "$1" in
		-i|--input) input_bam=$2; shift 2;;
		-t|--cell-cluster-file) cluster=$2; shift 2;;
		-o|--output) output=$2; shift 2;;
		-p|--thread) nthread=$2; shift 2;;
		-h|--help) echo -e $show_usage; exit 0;;
		--) break;;
		*) echo -e $1,$2,$show_usage; break;;
	esac
done

if [[ -z ${input_bam} || -z ${cluster} ]]; then
        echo -e $show_usage
        exit 0
else
	if [ -z ${output} ]; then
		output=`pwd`"/pseudo_bulk"
	fi
	
	if [[ -z ${nthread} ]]; then
		nthread=8
	fi

	mkdir -p ${output}/bam ${output}/bw
	samtools view -@ ${nthread} ${input_bam} | awk '{print substr($1,length($1)-21,11), $1}' OFS='\t' > tmp.reads
	sort --parallel ${nthread} tmp.reads > tmp.sorted.reads
	sed 's/ /_/g' ${cluster} | sort > tmp.cluster

	join -j 1 tmp.sorted.reads tmp.cluster | awk '{print $2, $3}' OFS='\t' > tmp.Cluster.Reads
	sort --parallel ${nthread} tmp.Cluster.Reads > tmp.Cluster.sorted.Reads

	samtools view -H ${input_bam} > tmp.header

	samtools view -@ ${nthread} ${input_bam} > tmp.all.reads
	sort --parallel ${nthread} tmp.all.reads > tmp.all.sorted.reads
	join -j 1 tmp.Cluster.sorted.Reads tmp.all.sorted.reads | sed 's/ /\t/g' > tmp.all.clusters.reads

	for i in `cut -f 2 tmp.cluster | sort -u `; do
		awk '$2=="'${i}'"' tmp.all.clusters.reads | cut -f 1,3- | cat tmp.header - | samtools view -@ ${nthread} -bh | samtools sort -m 4G -@ ${nthread} > ${output}/bam/DNA_Cluster${i}.bam
		nLine=`samtools view -@ ${nthread} ${output}/bam/DNA_Cluster${i}.bam | wc -l `
		nFrag=`echo "scale=1; $nLine/1" | bc`
		samtools view -H ${output}/bam/DNA_Cluster${i}.bam | awk '$1 == "@SQ" {OFS="\t";print $2,$3}' - | sed 's/.N://g' > tmp.size
		bedtools genomecov -ibam ${output}/bam/DNA_Cluster${i}.bam -bg | awk '{printf "%s\t%.0f\t%.0f\t%.2f\n",$1,$2,$3,$4*1000000/"'$nFrag'"}' | awk '{$4/=1;print}' OFS='\t' > tmp.bdg
		sort --parallel ${nthread} -k1,1 -k2,2n tmp.bdg > tmp.sorted.bdg
		bedGraphToBigWig tmp.sorted.bdg tmp.size ${output}/bw/DNA_Cluster${i}.bw
	done
	rm -f tmp.all.clusters.reads tmp.all.reads tmp.all.sorted.reads tmp.cluster tmp.Cluster.Reads tmp.Cluster.sorted.Reads tmp.header tmp.reads tmp.sorted.reads tmp.size tmp.bdg tmp.sorted.bdg
fi
