#!/bin/python

# sample_barcode.txt
'''
barcode sample
01      H3K27ac
02      H3K27ac
03      H3K27ac
04      H3K27ac
05      H3K27ac
06      H3K27ac
07      H3K9me3
08      H3K9me3
09      H3K9me3
10      H3K9me3
11      H3K9me3
12      H3K9me3
'''

import argparse
import pandas as pd
from multiprocessing import Pool

parser = argparse.ArgumentParser(description='filter bam based on QNAMES')
parser.add_argument('--bam', type=str, dest="bamf", help='bam file prefix')
parser.add_argument('--statH', type=str, dest="statH", help='input statH matrix')
parser.add_argument('-o', '--outPrefix', type=str, dest="outPrefix", help='output prefix')
parser.add_argument('-p', '--cores',type=int,dest="cpu",help='num of CPUs',default=30)
args = parser.parse_args()

import numpy as np
import pysam
from time import perf_counter as pc

NCPU = args.cpu


def run():
  start_time = pc()
  """ init input files """
  bamf = args.bamf
  statHf = args.statH
  outPrefix = args.outPrefix
  print("filter out bam files")
  generate_bams(bamf, statHf, outPrefix)
  end_time = pc()
  print('Used (secs): ', end_time - start_time)

def generate_bams(bamf, statHf, prefix):
  o_stat_H = pd.read_csv(statHf, sep='\t')
  print(np.unique(o_stat_H['sample']))
  p = Pool(NCPU)
  for idx in np.unique(o_stat_H['sample']):
        p.apply_async(generate_bam_worker, (bamf, o_stat_H, idx, prefix))
  p.close()
  p.join()

def generate_bam_worker(bamf, o_stat_H, sample, prefix):
  print(f"hello {sample}")
  name = bamf 
  bamF = pysam.AlignmentFile(name)
  qnames = list(o_stat_H[o_stat_H['sample']==sample]['barcode'].astype(str))
  nbarcode = len(qnames[0].split(':'))
  qnames_set = set(qnames)
  bam_fname = prefix + "." + sample + ".bam"
  print("For metaCell =", sample, "The filtered bam is writing to:", bam_fname)
  obam = pysam.AlignmentFile(bam_fname, "wb", template=bamF)
  for b in bamF.fetch(until_eof=True):
    if ':'.join(b.query_name.split(':')[-1-nbarcode:-1]) in qnames_set:
        obam.write(b)
  obam.close()
  bamF.close()
  print("metaCell =", sample," witting finished.")



if __name__ == "__main__":
  """filter bam based on QNAMES"""
  run()

