# import standard python libraries
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os, subprocess
# Import python package for working with cooler files and tools for analysis
import cooler
import cooltools.lib.plotting
from packaging import version
import cooltools
import bioframe
import argparse

### get args
parser = argparse.ArgumentParser(description="This script is for calculate PC1 using cooltools")
parser.add_argument("-f","--file",help="REQUIRED, cooler files that you want to calculate on, like ./HiRes_pseudobulk.mcool",dest="file")
parser.add_argument("-r","--resolution",help="REQUIRED, cooler file resolution you want to use, default is 100000",dest="res",default=100000)
args = parser.parse_args()

resolution=args.res
cooler_file=args.file

mm10_genome = bioframe.load_fasta('/storage/zhangyanxiaoLab/share/fasta/mm10.fa');

base_name = os.path.basename(cooler_file).split('_c')[0]
print(base_name)
# Load the cooler file
clr = cooler.Cooler(f'{cooler_file}::resolutions/{resolution}')
# Get the bins
bins = clr.bins()[:]
# Compute GC content
gc_cov = bioframe.frac_gc(bins[['chrom', 'start', 'end']], mm10_genome)
# Define the view DataFrame
view_df = pd.DataFrame({'chrom': clr.chromnames,'start': 0,
                        'end': clr.chromsizes.values,'name': clr.chromnames})
# Compute eigenvectors
cis_eigs = cooltools.eigs_cis(clr,gc_cov,view_df=view_df,n_eigs=3,)

# Extract the eigenvector track
eigenvector_track = cis_eigs[1]
  
# Save the eigenvector track to a text file
  
output_file = f'{base_name}_eigenvector_{resolution}.tsv'
print(f'Saved eigenvector track to {output_file}')
eigenvector_track.to_csv(output_file, index=False, sep='\t')
