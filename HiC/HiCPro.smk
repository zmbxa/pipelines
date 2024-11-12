### This snakemake file can process HiC data data ###
### The workdir need to be changed as data containing folder (include fastq/) and all .fq file in fastq/[name]/[name]_R1/R2.fastq.gz ###
### !!! Change SAMPLES to a Tuple contain sample names (for snakemake final output)!!! ###
### workflow is similar to HiC-Pro but will work in parallel ###
### Usage:  snakemake --cores 4 -s HiCPro.smk --use-conda --conda-frontend conda

import pandas as pd
import os
import re

f = [file for file in os.listdir("fastq") if os.path.isfile(os.path.join("fastq", file))]
m =[ re.search("(.*?)(|_1|_2|_r1|_r2|_R1|_R2)(_001)?.fastq.*",x) for x in f]
name = [x.group(1) for x in m if x is not None]
#
FASTQ_DICT = dict()
for idx in range(len(f)):
  if name[idx] in FASTQ_DICT:
    FASTQ_DICT[name[idx]].append("fastq/"+f[idx])
  else:
    FASTQ_DICT[name[idx]] = ["fastq/"+f[idx]]

for sample in FASTQ_DICT:
  FASTQ_DICT[sample].sort()
SAMPLES = FASTQ_DICT.keys()
GENOME = config["GENOME"]
CHROM_SIZE="/storage/zhangyanxiaoLab/niuyuxiao/tools/HiC-Pro-3.1.0/annotation/chrom_"+GENOME+".sizes"

rule all:
  input:
    allVP=expand("HiCPro_out/hic_results/data/{sample}/{sample}.allValidPairs",sample=SAMPLES),
    mpairStat=expand("HiCPro_out/hic_results/stats/{sample}/{sample}_allValidPairs.mergestat",sample=SAMPLES),
    avp_stat=expand("HiCPro_out/hic_results/stats/{sample}/{sample}_allValidPairs.mergestat",sample=SAMPLES),
    hicfile=expand("hicFile_juicer/{sample}_srt.hic",sample=SAMPLES),
    qc="HiCPro_out/all_sample_qc.txt",

rule HiCpro:
  input:
    inputFolder="fastq_hicpro/{sample}/{sample}/",
    config="config-hicpro.txt",
  output:
    allVP="HiCPro_out/hic_results/data/{sample}/{sample}.allValidPairs",
    mpairStat="HiCPro_out/hic_results/stats/{sample}/{sample}.mpairstat",
    avp_stat="HiCPro_out/hic_results/stats/{sample}/{sample}_allValidPairs.mergestat",
  shell:
    '''
    /storage/zhangyanxiaoLab/niuyuxiao/tools/HiC-Pro-3.1.0/bin/HiCPro4smk -i {input.inputFolder} -o HiCPro_out -c {input.config}
    '''

rule avp2hic:
  input:
    allVP="HiCPro_out/hic_results/data/{sample}/{sample}.allValidPairs",
  params:
    chromSize=CHROM_SIZE,
  output:
    hicfile="hicFile_juicer/{sample}_srt.hic",
    juicebox=temp("hicFile_juicer/{sample}_allValidPairs.pre_juicebox_sorted"),
  shell:
    '''
    echo "Generating Juicebox input files for ${wildcards.sample}..."
    awk '{{$4=$4!="+"; $7=$7!="+"}} $2<=$5{{print $1, $4, $2, $3, 0, $7, $5, $6, 1, $11, $12 }}$5<$2{{ print $1, $7, $5, $6, 0, $4, $2, $3, 1, $12, $11 }}' {input.allVP} | sort -T ./ -k3,3d  -k7,7d -S 250G  > {output.juicebox}
    
    echo "Running Juicebox for {wildcards.sample}..."
    java -XX:ParallelGCThreads=20 -Xmx50g -jar ~/tools/juicer/juicer/CPU/common/juicer_tools_1.22.01.jar pre {output.juicebox} {output.hicfile} {params.chromSize}
    '''

rule QC:
  input:
    mpairStat=expand("HiCPro_out/hic_results/stats/{sample}/{sample}_allValidPairs.mergestat",sample=SAMPLES),
    avp_stat=expand("HiCPro_out/hic_results/stats/{sample}/{sample}_allValidPairs.mergestat",sample=SAMPLES),
  output:
    qc="HiCPro_out/all_sample_qc.txt"
  shell:
    '''
    cd HiCPro_out/
    python /storage/zhangyanxiaoLab/niuyuxiao/pipelines/HiC/QC_HiCPro.py  -s {SAMPLES}
    cd ..
    python ~/pipelines/reminder.py "Youcan check your hic qc now!" "$(cat {output.qc})"
    '''



