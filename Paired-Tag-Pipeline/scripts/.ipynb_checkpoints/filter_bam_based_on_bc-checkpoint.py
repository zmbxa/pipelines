#!/bin/python

import argparse
import numpy as np
import pandas as pd
import pysam
from time import perf_counter as pc

parser = argparse.ArgumentParser(description='filter bam based on QNAMES')
parser.add_argument('--bam', type=str, dest="bamf", help='bam file prefix')
parser.add_argument('--statH', type=str, dest="statH", help='input statH matrix')
args = parser.parse_args()

def run():
    
    bamf = args.bamf
    statHf = args.statH

    qnames_set = set(pd.read_csv(statHf, sep='\t')["BC"])
    bamF = pysam.AlignmentFile(bamf)
    bam_fname = bamf.split('.bam')[0] + "_filtered_bc.bam"  
    obam = pysam.AlignmentFile(bam_fname, "wb", template=bamF)

    for b in bamF.fetch(until_eof=True):
        if ":".join(b.query_name.split(':')[-4:-1]) in qnames_set:
            obam.write(b)
    obam.close()
    bamF.close()
    

if __name__ == "__main__":
    """filter bam based on barcode"""
    run()
