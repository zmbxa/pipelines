### This scripts is for parallel HiCPro running step by step
usage(){
  echo "Usage: bash $0 -g <genome> -e <enzyme> -p <parallel sample numbers> -n <N_CPU for config bowtie2> -m <RAM for samtools sorting>"
  echo "Options:"
  echo "  -g,   Genome name in mm10, hg38 or hg19"
  echo "  -p,   Core numbers for snakemake, default is 4"
  echo "  -n,   number of CPUs using in bowtie2 mapping for each sample, default is 16"
  echo "  -e,   restrict enzyme used in the experiment, one of [DpnII, HaeIII or hindIII], default is DpnII"
  echo "  -m,   Mb, sorting mamory limit of each CPU during samtools sort, default is 2000M "
  exit 1
}

N_CORE=4
N_CPU=16
ENZYME=DpnII
ST_RAM=2000M
GENOME=mm10
while getopts ":g:p:n:e:m:" opt; do
  case ${opt} in
    g)GENOME=$OPTARG;;
    p)N_CORE=$OPTARG;;
    n)N_CPU=$OPTARG;;
    e)ENZYME=$OPTARG;;
    m)ST_RAM=$OPTARG;;
    *)usage;;
  esac
done

## modify config file
cp /storage/zhangyanxiaoLab/niuyuxiao/pipelines/HiC/config-hicpro.txt ${PWD}/
sed "s/CUSTOMGENOME/${GENOME}/g" config-hicpro.txt -i 
sed "s/CUSTOMEnzyme/${ENZYME}/g" config-hicpro.txt -i 
if [ $ENZYME != DpnII ];then
 case $ENZYME in
    hindIII )
      sed 's/GATCGATC/AAGCTAGCTT/g' -i config-hicpro.txt ;;
    HaeIII )
      sed 's/GATCGATC/GGCCGGCC/g' -i config-hicpro.txt ;;
    *)
      echo "Enzyme not recognized!"
      exit 1 ;;
  esac
fi

if [ $N_CPU != 16 ];then
  sed "s/^N_CP.*$/N_CPU = ${N_CPU}/g" config-hicpro.txt -i 
fi
if [ $ST_RAM != 2000M ];then
  sed "s/^SORT_RAM.*$/SORT_RAM = ${ST_RAM}/g" config-hicpro.txt -i 
fi

### make fastq tree directory
for sample in $(ls fastq/|cut -f 1 -d _|sort|uniq);do
mkdir fastq_hicpro/${sample}/ -p 
ln -s $(realpath fastq/${sample}*R1*gz) fastq_hicpro/${sample}/${sample}_R1.fastq.gz 
ln -s $(realpath fastq/${sample}*R2*gz)  fastq_hicpro/${sample}/${sample}_R2.fastq.gz
done

source /storage/zhangyanxiaoLab/niuyuxiao/anaconda3/etc/profile.d/conda.sh
conda activate snakemake

snakemake --cores $N_CORE -s ~/pipelines/HiC/HiCPro_step.smk --config GENOME=${GENOME}  ENZYME=${ENZYME} -p --use-conda --conda-frontend conda --dryrun
if [ $? -eq 0 ]; then
  echo done
snakemake --cores $N_CORE -s ~/pipelines/HiC/HiCPro_step.smk --config GENOME=${GENOME}  ENZYME=${ENZYME} -p --use-conda --conda-frontend conda
else
  echo "Please Check inputs or see runing log..."
  exit 1
fi

if [ $? -eq 0 ]; then
 python /storage/zhangyanxiaoLab/niuyuxiao/pipelines/reminder.py "Your HiC data mapping is successfully done! Congrats~" "Need a cup of latte?"
else
 echo "Please Check inputs or see runing log..."
 exit 1
fi






