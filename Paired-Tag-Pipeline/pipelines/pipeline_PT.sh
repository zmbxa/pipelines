#!/bin/bash
##-------------------------
# @Author: XiongXiong, YanxiaoZhang and Yuxiao Niu
# @Date: 2024/5/30
# @Description: original pipeline from Xiong, remove dup part was modified by Yanxiao. Then add CreGFP into reference.
##-------------------------

script_dir=$(dirname $0)
show_usage="Usage: bash pipeline_PT.sh -g <genome> -l <library> -p <thread>\n\n
			-g|--reference-genome \t reference genome, must be one of hg38, mm10, mm10_gfp or ce11\n
			-l|--library-type \t library type, must be one of DNA and RNA\n
			-p|--thread \t how many threads to use. Default:8\n
			-h|--help \t\t get help document\n"

GETOPT_ARGS=`getopt -o g:l:p:d:h -al reference-genome:,library-type:,thread:,debug:,help -- "$@"`
eval set -- "$GETOPT_ARGS"

while [ -n "$1" ]
do
	case "$1" in
		-g|--reference-genome) reference_genome=$2; shift 2;;
		-l|--library-type) library=$2; shift 2;;
		-p|--thread) nthread=$2; shift 2;;
    -d|--debug) debug=$2; shift 2;;
		-h|--help) echo -e $show_usage; exit 0;;
		--) break;;
		*) echo -e $1,$2,$show_usage; break;;
	esac
done

if [[ -z ${reference_genome} || -z ${library} ]]; then
        echo -e $show_usage
        exit 0
else
	
	if [[ -z ${nthread} ]]; then
		nthread=8
	fi

  if [[ -z ${debug} ]]; then
    debug=-1
  fi


	rename _001.fastq.gz .fastq.gz ./fastq/*
  ls ./fastq/*_R1.fastq.gz | sed "s/_R1.fastq.gz//g; s/\.\/fastq\///g" | sort --parallel ${nthread} -u | xargs --max-procs=4 -n 1 -I {} sh -c "bash ${script_dir}/PT.main.sh {} $reference_genome $library $nthread $debug"

  echo "All samples completed."

fi
