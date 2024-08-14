### This script is for paired-ended ChIP/CUT&Tag/ATAC data mapping. 
## trim adapter --> Bowtie2 alignment --> sambamba mark&remove duplicate  --> generate bigWig file

usage() {
  echo "Usage: bash $0 -g <genome> -c <cores_for_snakemake>"
  echo "Options:"
  echo "  -g,   Genome name in mm10, mm10_GFP, dm6, hg38,ce11,and lambda"
  echo "  -c,   Core numbers for snakemake, default is 4"
  echo "  -s,   Spike-in genome (optional)"
  exit 1
}
N_CORE=4
while getopts ":g:c:s:" opt; do
  case ${opt} in
    g)GENOME=$OPTARG;;
    c)N_CORE=$OPTARG;;
    s)SPIKE=$OPTARG;;
    *)usage;;
  esac
done

if [ -z "$GENOME" ]; then
  usage
fi

if [ -z "$SPIKE" ]; then
 echo "Mapping ATAC/PE_ChIP data to $GENOME genome using $N_CORE cores."
else
 echo "Mapping ATAC/PE_ChIP data to $GENOME genome with $SPIKE spike-in, using $N_CORE cores."
fi

#  

source /storage/zhangyanxiaoLab/niuyuxiao/anaconda3/etc/profile.d/conda.sh
conda activate snakemake

if [ -z "$SPIKE" ]; then
  snakemake --cores $N_CORE -s /storage/zhangyanxiaoLab/niuyuxiao/pipelines/ATAC_PEChIP_seq/ATACseq.smk --config GENOME=$GENOME --use-conda --conda-frontend conda --dryrun
  if [ $? -eq 0 ]; then
   snakemake --cores $N_CORE -s /storage/zhangyanxiaoLab/niuyuxiao/pipelines/ATAC_PEChIP_seq/ATACseq.smk --config GENOME=$GENOME --use-conda --conda-frontend conda 
  else
    echo "Please Check inputs or see runing log..."
  exit 1
  fi
else
 snakemake --cores $N_CORE -s /storage/zhangyanxiaoLab/niuyuxiao/pipelines/ATAC_PEChIP_seq/ATACseq_PEChIP_unmap2SpikeIn.smk --config GENOME=$GENOME spikeGENOME=$SPIKE --use-conda --conda-frontend conda --dryrun
 if [ $? -eq 0 ]; then
    snakemake --cores $N_CORE -s /storage/zhangyanxiaoLab/niuyuxiao/pipelines/ATAC_PEChIP_seq/ATACseq_PEChIP_unmap2SpikeIn.smk --config GENOME=$GENOME spikeGENOME=$SPIKE --use-conda --conda-frontend conda 
  else
    echo "Please Check inputs or see runing log..."
  exit 1
  fi
fi


if [ $? -eq 0 ]; then
 python /storage/zhangyanxiaoLab/niuyuxiao/pipelines/reminder_westlake.py "Your bulk ATAC/PE-ChIP mapping are successfully done! Congrats~" "Need a cup of Latte?"
else
 echo "Please Check mail setting..."
 exit 1
fi