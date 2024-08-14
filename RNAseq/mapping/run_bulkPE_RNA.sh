### This script is for paired-ended bulk RNA-seq data mapping. 
## fastp --> STAR  alignment --> featureCount  --> generate counts matrix file || fastp --> STAR  alignment --> TEcount -->generate TEcounts matrix

usage() {
  echo "Usage: bash $0 -g <genome> -c <cores_for_snakemake> -noTE"
  echo "Options:"
  echo "  -g,   Genome name in hg38, hg38_with_dm6, mm10, mm10_GFP, mm10_with_dm6(in1), dm6 or Zokor3"
  echo "  -c,   Core numbers for snakemake, default is 4"
  echo "  -n, --noTE,   Do not mapping to TE reference"
  exit 1
}

NO_TE=false
N_CORE=4
while getopts ":g:c:noTE:n" opt; do
  case ${opt} in
    g)GENOME=$OPTARG;;
    c)N_CORE=$OPTARG;;
    n)NO_TE=true;;
    noTE)NO_TE=true;;
    *)usage;;
  esac
done

if [ -z "$GENOME" ]; then
  usage
fi

source /storage/zhangyanxiaoLab/niuyuxiao/anaconda3/etc/profile.d/conda.sh
conda activate snakemake

if $NO_TE; then
  echo "Mapping bulk RNAseq data to $GENOME genome with $N_CORE cores without TE"
  SMK_FILE=/storage/zhangyanxiaoLab/niuyuxiao/pipelines/RNAseq/mapping/bulk_PE_RNA.smk
else
  echo "Mapping bulk RNAseq data to $GENOME genome with $N_CORE cores with TE"
  SMK_FILE=/storage/zhangyanxiaoLab/niuyuxiao/pipelines/RNAseq/mapping/bulk_PE_RNA_TE.smk
fi

snakemake --cores $N_CORE -s ${SMK_FILE} --config GENOME=$GENOME --use-conda --conda-frontend conda --dryrun
if [ $? -eq 0 ]; then
 snakemake --cores $N_CORE -s ${SMK_FILE} --config GENOME=$GENOME --use-conda --conda-frontend conda 
else
  echo "Please Check inputs or see runing log..."
  exit 1
fi

if [ $? -eq 0 ]; then
 python /storage/zhangyanxiaoLab/niuyuxiao/pipelines/reminder.py "Your bulk RNAseq mapping are successfully done! Congrats~" "Need a cup of Latte?"
else
 echo "Please Check inputs or see runing log..."
 exit 1
fi



