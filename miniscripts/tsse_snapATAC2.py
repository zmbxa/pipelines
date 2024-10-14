import subprocess
# activating conda environment
conda_env = "/storage/zhangyanxiaoLab/niuyuxiao/anaconda3/envs/snapATAC2"
#import os 
#os.system("source /storage/zhangyanxiaoLab/niuyuxiao/anaconda3/etc/profile.d/conda.sh")
#os.system('conda activate ' + conda_env)
#activate_cmd = f"conda activate {conda_env}"
#subprocess.run(activate_cmd, shell=True,executable='/bin/bash')

import snapatac2 as snap
import argparse
import glob 
import os  
import anndata as ad
import gzip
import pandas as pd

parser = argparse.ArgumentParser(description="This script is for calculating TSS-enrichment score from bam files using snapATAC2.")
parser.add_argument("-b","--bam",nargs="+",help="REQUIRED, bam files that you want to assess, like ./passorted_bam.bam",dest="bam_files")
parser.add_argument("-s","--species",help="specify species for chrom_size to use, mm10 for default",default="mm10",dest="species")
parser.add_argument("-f","--fragments-dir",help="directory contain fragments file, default=NONE, this script will generate fragments file under ./tsse/frags",default="./tsse/frags",dest="frag_dir")
parser.add_argument("-O","--output_dir",help="output directory, ./tsse if not specified",default="./tsse",dest="output_dir")
parser.add_argument("-o","--output_tsse",help="output tsse score file name, tsse.csv if not specified",default="tsse.csv",dest="tsseOut")
parser.add_argument("-t","--thread",help="number od thread,default is 10",default=10)
args = parser.parse_args()

if not args.bam_files:
  print("BAM file is required!")
  parser.print_usage()
  exit(1)

if not hasattr(snap.genome,args.species):
  print("This species is not available in snapATAC2!")
  parser.print_usage()
  exit(1)


if not os.path.exists(args.output_dir):
  os.makedirs(args.output_dir)
  
if not os.path.exists(args.frag_dir):
  os.makedirs(args.frag_dir)
 
# make a df to save tsse score, indexes are sample_id prefixes
result_df = pd.DataFrame(columns=['index', 'tsse'])

### begin analysis, for each inputed bam file
for f in args.bam_files:
  name=f.rsplit("/",1)[-1].split(".bam")[0]
  fragFile=args.frag_dir+"/"+name+".frag.tsv.gzip"
  print("Generating fragments file: "+fragFile)
  # generate gzip fragment to FragFile, for bulk CUT&Tag regex
  snap.pp.make_fragment_file(f,fragFile,compression='gzip',barcode_regex= r"^(.*?):\d+.*$")
  
  # import data for TSSe calculation
  data = snap.pp.import_data(fragment_file=fragFile,chrom_sizes=getattr(snap.genome,args.species),sorted_by_barcode=False)
  
  # calc tsse, value will store to the 'data' object
  print("Calculating tsse...")
  snap.metrics.tsse(data, getattr(snap.genome,args.species))
  
  result_df.loc[len(result_df)] = [name,data.obs['tsse'].values[0]]  
  

result_df.to_csv(args.output_dir+"/"+args.tsseOut, index=False, mode='a', header=True)
  
  
  
  
