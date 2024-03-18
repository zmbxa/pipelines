### This script is for auto-merge 2 batches data, 2 paths of each batch and sublibrary start are required. 
### It will merge same id file of shallow and deep batch in ./fastq folder
### Usage:
###	bash ~/pipelines/Paired-Tag-Pipeline/mergefile.sh [path/to/shallowData/] [path/to/deepData/] [GY]
### 

#!/bin/bash
help() {
    sed -rn 's/^### ?//;T;p;' "$0"
}

if [[ -z "$1" || -z "$2" || -z "$3" ]] || [[ "$1" == "-h" ]]|| [[ "$1" == "--help" ]]; then
  help
  exit 1
fi

shallow=$1
deep=$2
pattern=$3



mkdir -p fastq
echo R1,R2,library,prefix,tissue,exp_prefix > sampleInfo.csv
for i in $(ls $shallow| grep "^${pattern}");do
# merge fq
echo $i  R1 merging...
cat ${shallow}/${i}/*R1*gz ${deep}/${i}/*R1*gz > fastq/${i}_R1.fastq.gz
echo  R2 merging...
zcat ${shallow}/${i}/*R2*gz ${deep}/${i}/*R2*gz |pigz -p 6 > fastq/${i}_R2.fastq.gz
path=$(echo $(ls $(pwd)/fastq/$i* )|sed "s/ /,/g")
if [ $(ls $shallow| grep "^${pattern}"|grep $i -n|cut -d ":" -f 1) -ge $(( $(ls $shallow| grep "^${pattern}"|wc -l) / 2 )) ];then 
  lib=RNA;exp=${pattern}${i##*${pattern}}
else 
  lib=DNA;exp=$(( ${i##*${pattern}} + $(( $(ls $shallow| grep "^${pattern}"|wc -l) / 2 )) ))
fi
echo $path,$lib,$i,brain_FC,$exp >>sampleInfo.csv
done

