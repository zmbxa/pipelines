#!/bin/bash
##-------------------------
# @Author: XiongXiong
# @Date: 2023/2/2
##-------------------------

show_usage="Usage: bash filter_merge_bam.sh -i <input_file_list> -d <DNA_Cutoff> -r <RNA_Cutoff> -p <thread>\n
			e.g. bash filter_merge_bam.sh -i batch1_file_list.txt\n\n
			-i|--input \t input file list, paired bam files paths should be separated by comma with DNA first in one line;    e.g. /storage/zhangyanxiaoLab/xiongxiong/Data/projects/mouse_brain/paired-tag/batch1/DNA/HSQ046_bam/HSQ046_UsefulReads.bam,/storage/zhangyanxiaoLab/xiongxiong/Data/projects/mouse_brain/paired-tag/batch1/RNA/HSQ050_bam/HSQ050_UsefulReads.bam\n\n
			-d|--dnacutoff \t only cells with read number above <DNA_Cutoff> are kept\n\n
			-r|--rnacutoff \t only cells with read number above <RNA_Cutoff> are kept\n\n
			-p|--thread \t how many threads to use. Default:8\n
			-h|--help \t\t get help document\n"

GETOPT_ARGS=`getopt -o i:d:r:p:h -al input:,dnacutoff:,rnacutoff:,thread:,help -- "$@"`
eval set -- "$GETOPT_ARGS"

while [ -n "$1" ]
do
	case "$1" in
		-i|--input) input_file=$2; shift 2;;
		-d|--dnacutoff) dnacutoff=$2; shift 2;;
		-r|--rnacutoff) rnacutoff=$2; shift 2;;
		-p|--thread) nthread=$2; shift 2;;
		-h|--help) echo -e $show_usage; exit 0;;
		--) break;;
		*) echo -e $1,$2,$show_usage; break;;
	esac
done

if [[ -z ${input_file} ]]; then
        echo -e $show_usage
        exit 0
else
	if [ -z ${dnacutoff} ]; then
		dnacutoff=500
	fi

	if [ -z ${rnacutoff} ]; then
		rnacutoff=1000
	fi
	
	if [[ -z ${nthread} ]]; then
		nthread=8
	fi

	nPairs=`cat ${input_file} | wc -l `
	for ((pair=1; pair<=${nPairs}; pair++)); do
		dnaBam=`awk -F "," '{if(NR=="'${pair}'") print $1}' ${input_file} `
		rnaBam=`awk -F "," '{if(NR=="'${pair}'") print $2}' ${input_file} `
		samtools view -@ ${nthread} ${dnaBam} | awk '{{print substr($1, length($1)-18, 8), $1}}' OFS='\t' > tmp.names
		sort --parallel ${nthread} tmp.names | bedtools groupby -g 1 -c 2 -o count -i - > tmp.dna.Counts
		sort --parallel ${nthread} tmp.dna.Counts > tmp.dna.sorted.Counts

		samtools view -@ ${nthread} ${rnaBam} | awk '{{print substr($1, length($1)-18, 8), $1}}' OFS='\t' > tmp.names
		sort --parallel ${nthread} tmp.names | bedtools groupby -g 1 -c 2 -o count -i - > tmp.rna.Counts
		sort --parallel ${nthread} tmp.rna.Counts > tmp.rna.sorted.Counts

		join -j 1 tmp.dna.sorted.Counts tmp.rna.sorted.Counts | awk -v dnacutoff=${dnacutoff} -v rnacutoff=${rnacutoff} '{{if($2>=dnacutoff && $3>=rnacutoff) print $1}}' > tmp.counts

		mkdir -p ./filtered_bam/
		
		pairN=`printf "%02d" ${pair} `
		samtools view -H ${dnaBam} > tmp.header
		samtools view -@ ${nthread} ${dnaBam} | awk '{{print substr($1,length($1)-18,8), $0}}' OFS='\t' | grep -w -F -f tmp.counts - | awk -v pairN=${pairN} '{{print substr($2,1,length($2)-20)":"pairN":"substr($2,length($2)-18,19), $0}}' OFS='\t' | cut -f 1,4- | cat tmp.header - | samtools view -@ ${nthread} -bh | samtools sort -m 4G -@ ${nthread} > ./filtered_bam/sub_lib${pairN}_DNA.bam
		samtools view -H ${rnaBam} > tmp.header
		samtools view -@ ${nthread} ${rnaBam} | awk '{{print substr($1,length($1)-18,8), $0}}' OFS='\t' | grep -w -F -f tmp.counts - | awk -v pairN=${pairN} '{{print substr($2,1,length($2)-20)":"pairN":"substr($2,length($2)-18,19), $0}}' OFS='\t' | cut -f 1,4- | cat tmp.header - | samtools view -@ ${nthread} -bh | samtools sort -m 4G -@ ${nthread} > ./filtered_bam/sub_lib${pairN}_RNA.bam
	done

	samtools merge -f -h ${rnaBam} -@ ${nthread} ./filtered_bam/merged_RNA.bam ./filtered_bam/sub_lib[0-9]*_RNA.bam
	samtools merge -f -h ${dnaBam} -@ ${nthread} tmp.bam ./filtered_bam/sub_lib[0-9]*_DNA.bam
	bedtools bamtobed -i tmp.bam | cut -f 1,2,4,6 > tmp.reads
	sort --parallel ${nthread} -k1,1 -k2,2n -k4,4 -k3,3 tmp.reads | bedtools groupby -i - -g 1,2,4 -c 3 -o distinct | awk '{split($4,names,",")} {if(length(names)>10) print $4}' | sed 's/,/\n/g' | sort --parallel ${nthread} -u > tmp.RmNames
	samtools view -@ ${nthread} tmp.bam | grep -v -w -F -f tmp.RmNames - | cat tmp.header - | samtools view -@ ${nthread} -bh | samtools sort -m 4G -@ ${nthread} > ./filtered_bam/merged_DNA.bam
	rm tmp.bam tmp.counts tmp.dna.Counts tmp.dna.sorted.Counts tmp.header tmp.names tmp.reads tmp.RmNames tmp.rna.Counts tmp.rna.sorted.Counts
fi
