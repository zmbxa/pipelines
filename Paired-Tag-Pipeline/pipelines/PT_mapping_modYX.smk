#snakefile 

import pandas as pd
import os
import re

meta = pd.read_csv("sample.csv")

print(meta)
SAMPLES=list(meta['sampleID'])
DNAsamples=list(meta.query('Type=="DNA"')['sampleID'])
RNAsamples=list(meta.query('Type=="RNA"')['sampleID'])

GENOME = config['GENOME']
#NCORE = config[NCORE]

BAM2MTX='/storage/zhangyanxiaoLab/niuyuxiao/pipelines/Paired-Tag-Pipeline/pipelines/bam2mtx_gfp_AAVCRE_ERT.sh' if GENOME=='mm10_gfp'  else "/storage/zhangyanxiaoLab/niuyuxiao/pipelines/Paired-Tag-Pipeline/pipelines/bam2mtx.sh"

rule require_all:
  input:
#    "merge/bam_list.txt",
    expand("DNA/{sample}/bam/{sample}_UsefulReads.bam",sample=DNAsamples),
    expand("DNA/{sample}/process/{sample}_CB.Counts",sample=DNAsamples),
    expand("RNA/{sample}/bam/{sample}_UsefulReads.bam",sample=RNAsamples),
    expand("RNA/{sample}/process/{sample}_CB.Counts",sample=RNAsamples),
    "merge/filtered_bam/merged_DNA.bam",
    "merge/filtered_bam/merged_RNA.bam",
    "merge/DNA_matrix/matrix.mtx.gz",
    "merge/RNA_matrix/matrix.mtx.gz",
    "03_Sub-lib_Read_Number.pdf",
    "04_Barcode_Cumulative_Prop.pdf",
    "reportSummary.csv",
    "MapQuality_Bar.pdf",
    "RNA_barcodesDistribution.pdf",
    "DNA_barcodesDistribution.pdf",

rule mapping_DNA:
  input:
    expand("DNA/fastq/{sample}_R1.fastq.gz",sample=DNAsamples),
  output:
    expand("DNA/{sample}/bam/{sample}_UsefulReads.bam",sample=DNAsamples),
    expand("DNA/{sample}/process/{sample}_CB.Counts",sample=DNAsamples),
  conda:
    "paired_tag"
  params:
    genome=GENOME,
  shell:
    '''
      cd DNA
      bash /storage/zhangyanxiaoLab/niuyuxiao/pipelines/Paired-Tag-Pipeline/pipelines/pipeline_PT.sh -l DNA -g {params.genome}
      cd ..
    '''

rule mapping_RNA:
  input:
    expand("RNA/fastq/{sample}_R1.fastq.gz",sample=RNAsamples),
  output:
    expand("RNA/{sample}/bam/{sample}_UsefulReads.bam",sample=RNAsamples),
    expand("RNA/{sample}/process/{sample}_CB.Counts",sample=RNAsamples),
  conda:
    "paired_tag"
  params:
    genome=GENOME,
  shell:
    '''
      cd RNA
      bash /storage/zhangyanxiaoLab/niuyuxiao/pipelines/Paired-Tag-Pipeline/pipelines/pipeline_PT.sh -l RNA -g {params.genome}
      cd ..
    '''

rule list_make:
    output:
        "bam_list.txt"
    run:
        path=os.getcwd()
        for row in range(len(DNAsamples)):
          f=open("bam_list.txt","a")
          f.write(path+"/DNA/"+DNAsamples[row]+"/bam/"+DNAsamples[row]+"_UsefulReads.bam"+","+path+"/RNA/"+RNAsamples[row]+"/bam/"+RNAsamples[row]+"_UsefulReads.bam")
          f.write('\n')
        f.close()

rule merge:
  input:
    expand("DNA/{sample}/bam/{sample}_UsefulReads.bam",sample=DNAsamples),
    expand("RNA/{sample}/bam/{sample}_UsefulReads.bam",sample=RNAsamples),
    bamlist="bam_list.txt"
  output:
    "merge/filtered_bam/merged_DNA.bam",
    "merge/filtered_bam/merged_RNA.bam",
  conda:
    "paired_tag"
  shell:
    '''
    mkdir -p merge
    mv {input.bamlist} merge/bam_list.txt
    cd merge
    bash /storage/zhangyanxiaoLab/niuyuxiao/pipelines/Paired-Tag-Pipeline/pipelines/filter_merge_bam.sh -i bam_list.txt  -d 300 -r 500
    cd ..
    '''

rule bam2mtx:
  input:
    "merge/filtered_bam/merged_{type}.bam",
#    "merge/filtered_bam/merged_RNA.bam"
  output:
    "merge/{type}_matrix/matrix.mtx.gz",
#    RNA="merge/RNA_matrix/matrix.mtx.gz"
  conda:
    "paired_tag"
  params:
    bam2mtx=BAM2MTX
  shell:
    '''
    cd merge
    bash {params.bam2mtx} -n merged_{wildcards.type} -g mm10 -i {input} -l {wildcards.type} -o ./{wildcards.type}_matrix
    cd ..
    '''

rule QC_plots_xiong:
  input:
    expand("DNA/{sample}/process/{sample}_CB.Counts",sample=DNAsamples),
    expand("RNA/{sample}/process/{sample}_CB.Counts",sample=RNAsamples),
  output:
    "03_Sub-lib_Read_Number.pdf",
    "04_Barcode_Cumulative_Prop.pdf"
  run:
    command1="Rscript /storage/zhangyanxiaoLab/niuyuxiao/pipelines/Paired-Tag-Pipeline/analysis/dot_reads_numbers_total.R 300 500 "
    path=os.getcwd()
    for row in range(len(DNAsamples)):
      command1=command1+"DNA/"+DNAsamples[row]+"/process/"+DNAsamples[row]+"_CB.Counts,RNA/"+RNAsamples[row]+"/process/"+RNAsamples[row]+"_CB.Counts "
    command2="Rscript /storage/zhangyanxiaoLab/niuyuxiao/pipelines/Paired-Tag-Pipeline/analysis/line_barcodeProp_total.R "+(' '.join(["DNA/"+sample+"/process/"+sample+"_CB.Counts" for sample in DNAsamples]))+" "+(' '.join(["RNA/"+sample+"/process/"+sample+"_CB.Counts" for sample in RNAsamples]))
    os.system(command1)
    os.system(command2)

rule mapping_summary:
  input:
    expand("DNA/{sample}/process/{sample}_CB.Counts",sample=DNAsamples)
  output:
    "reportSummary.csv" 
  shell:
    '''
    bash /storage/zhangyanxiaoLab/niuyuxiao/pipelines/Paired-Tag-Pipeline/pipelines/ReportSummary.sh .
    '''

rule QC_plots_yuxiao:
  input:
    "reportSummary.csv"
  output:
    "MapQuality_Bar.pdf",
    "RNA_barcodesDistribution.pdf",
    "DNA_barcodesDistribution.pdf",
  shell:
    "Rscript /storage/zhangyanxiaoLab/niuyuxiao/pipelines/Paired-Tag-Pipeline/analysis/QC_bar_BCbox.R $PWD"
