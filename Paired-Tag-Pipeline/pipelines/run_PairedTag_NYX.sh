### This script is for paired-tag data mapping. 
## Original pipeline from Xiong. Remove dup part was modified by Yanxiao. Automatization by Yuxiao.

usage() {
  echo "Usage: bash $0 -g <genome> -c <cores_for_snakemake>"
  echo "Options:"
  echo "  -g,   Genome name in mm10, mm10_gfp, dm6, hg38, ce11"
  echo "  -c,   Core numbers for snakemake, default is 8"
  exit 1
}
N_CORE=8
while getopts ":g:c:s:" opt; do
  case ${opt} in
    g)GENOME=$OPTARG;;
    c)N_CORE=$OPTARG;;
    *)usage;;
  esac
done

if [ -z "$GENOME" ]; then
  usage
fi

echo "Mapping paired-tag data to $GENOME genome using $N_CORE cores."
begin_time=$(date +%F%n%T)
# make sample list
python /storage/zhangyanxiaoLab/niuyuxiao/pipelines/Paired-Tag-Pipeline/analysis/makeSampleInfo.py
# move files
for type in DNA RNA;do
mkdir -p $type/fastq
for sample in $(grep ${type} sample.csv|cut -f 1 -d ,);do 
  ln -s ${PWD}/fastq/${sample}*R1*.fastq.gz ${type}/fastq/${sample}_R1.fastq.gz
  ln -s ${PWD}/fastq/${sample}*R2*.fastq.gz ${type}/fastq/${sample}_R2.fastq.gz
done
done

source /storage/zhangyanxiaoLab/niuyuxiao/anaconda3/etc/profile.d/conda.sh
conda activate snakemake

snakemake --cores $N_CORE -s /storage/zhangyanxiaoLab/niuyuxiao/pipelines/Paired-Tag-Pipeline/pipelines/PT_mapping_modYX.smk --config GENOME=$GENOME --use-conda --conda-frontend conda --dryrun
if [ $? -eq 0 ]; then
  echo Fine! Start running
  snakemake --cores $N_CORE -s /storage/zhangyanxiaoLab/niuyuxiao/pipelines/Paired-Tag-Pipeline/pipelines/PT_mapping_modYX.smk --config GENOME=$GENOME --use-conda --conda-frontend conda 
else
  echo "Please Check inputs or see runing log..."
  exit 1
fi

if [ $? -eq 0 ]; then
 python /storage/zhangyanxiaoLab/niuyuxiao/pipelines/reminder_westlake.py "Your paired-tag mapping in ${begin_time} are successfully done! Congrats~" "See $PWD. Need a cup of Latte?"
else
 echo "Please Check mail setting..."
 exit 1
fi

