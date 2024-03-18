### This script is for paired-ended ChIP/CUT&Tag/ATAC data mapping. 
## trim adapter --> Bowtie2 alignment --> sambamba mark&remove duplicate  --> generate bigWig file

usage() {
  echo "Usage: bash $0 -g <genome> -c <cores_for_snakemake>"
  echo "Options:"
  echo "  -g,   Genome name in mm10, mm10_GFP, dm6 or hg38"
  echo "  -c,   Core numbers for snakemake, default is 4"
  exit 1
}
N_CORE=4
while getopts ":g:c:" opt; do
  case ${opt} in
    g)GENOME=$OPTARG;;
    c)N_CORE=$OPTARG;;
    *)usage;;
  esac
done
if [ -z "$GENOME" ]; then
  usage
fi

echo "Mapping ATAC/PE_ChIP data to $GENOME genome with $N_CORE cores."

source /storage/zhangyanxiaoLab/niuyuxiao/anaconda3/etc/profile.d/conda.sh
conda activate snakemake
snakemake --cores $N_CORE -s /storage/zhangyanxiaoLab/niuyuxiao/pipelines/ATAC_PEChIP_seq/ATACseq.smk --config GENOME=$GENOME --use-conda --conda-frontend conda --dryrun
if [ $? -eq 0 ]; then
 snakemake --cores $N_CORE -s /storage/zhangyanxiaoLab/niuyuxiao/pipelines/ATAC_PEChIP_seq/ATACseq.smk --config GENOME=$GENOME --use-conda --conda-frontend conda 
else
  echo "Please Check inputs or see runing log..."
 exit 1
fi

if [ $? -eq 0 ]; then
 python /storage/zhangyanxiaoLab/niuyuxiao/pipelines/reminder.py "Your bulk ATAC mapping are successfully done! Congrats~" "Need a cup of Latte?"
else
 echo "Please Check inputs or see runing log..."
 exit 1
fi