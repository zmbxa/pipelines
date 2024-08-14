#!/bin/bash
### This script is for auto-download data from SRA using aira2, a list txt file of SRA accession or metadata is required. 
### Usage:
###	bash ~/pipelines/miniscripts/download_SRA.sh [path/to/SRR_Acc_List.txt | path/to/SraRunTable.txt]
### 

help() {
    sed -rn 's/^### ?//;T;p;' "$0"
}

if [[ $1 == 0 ]] || [[ "$1" == "-h" ]]|| [[ "$1" == "--help" ]]; then
  help
  exit 1
fi

file=$1

if [[ $(grep SRR $file) == '' ]];then
  echo Your input file has no useful accession! please check it again!
  help
  exit 1
fi

if [[ $(head -1 $file |awk -F ,  '{print NF}') == 1 ]];then
  echo Your input file is : $file , in accession list format.
  awk -F, '{ print "https://sra-pub-run-odp.s3.amazonaws.com/sra/"$1"/"$1 }' $file > aria2_download.txt
else
  echo Your input file is : $file , in SraRunTable format.
  awk -F, 'FNR!=1 { print "https://sra-pub-run-odp.s3.amazonaws.com/sra/"$1"/"$1 }' $file > aria2_download.txt
fi


if [[ $(ls aria2_download.txt) == 0 ]];then
  help
  exit 1
else
  echo continue aria2 download
fi

 /storage/zhangyanxiaoLab/niuyuxiao/anaconda3/bin/aria2c -c -x 6 -s 6 -j 4 -d sra/ -i aria2_download.txt    

 mkdir -p fastq/
 fasterq-dump -3 -e 4 -p -O fastq/ sra/*
 pigz -p 6 fastq/* 


