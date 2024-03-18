### This snakemake file can process paired-end bulk RNAseq data (mouse mm10) ###
### The workdir need to be changed as data containing folder (include fastq/) and all .fq file in fastq/[name]/[name]_1/2.fq.gz ###
### !!! Change SAMPLES to a Tuple contain sample names (for snakemake final output)!!! ###
### workflow:  *** fastp ---> STAR ---> featureCounts ---> TECounts ***  output fetureCount and TECount matrix  ###
### Usage:  snakemake --cores 4 -s bulk_PE_RNA.smk --use-conda --conda-frontend conda

import pandas as pd
import os
import re

f = os.listdir("fastq")
m =[ re.search("(.*?)(|_1|_2|_r1|_r2|_R1|_R2)(_001)?.fastq.*",x) for x in f]
name = [x.group(1) for x in m]
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

rule all:
    input:
        "rna/bulk_rna_counts.txt",
        "te/te_counts.txt",
#        "test.txt"

rule fastp:
    input:
        r1="fastq/{sample}/{sample}_1.fastq.gz",
        r2="fastq/{sample}/{sample}_2.fastq.gz",
    output:
        r1_clean="fastq/{sample}/{sample}_R1.clean.fastq.gz",
        r2_clean="fastq/{sample}/{sample}_R2.clean.fastq.gz",
        report_json="fastq/{sample}/{sample}.fastp.json",
        report_html="fastq/{sample}/{sample}.fastp.html"
    shell:
        "fastp -i {input.r1} -I {input.r2} -o {output.r1_clean} -O {output.r2_clean} -j {output.report_json} -h {output.report_html}"


rule star_align:
    input:
        r1_clean="fastq/{sample}/{sample}_R1.clean.fastq.gz",
        r2_clean="fastq/{sample}/{sample}_R2.clean.fastq.gz"
    output:
        bam="rna/{sample}/{sample}.Aligned.sortedByCoord.out.bam"
    shell:
        "STAR --runThreadN 8 "
        "--genomeDir /storage/zhangyanxiaoLab/share/STAR_index/mm10 "
        "--readFilesIn {input.r1_clean} {input.r2_clean} "
        "--readFilesCommand zcat "
        "--outFileNamePrefix rna/{wildcards.sample}/{wildcards.sample}. "
        "--outSAMtype BAM SortedByCoordinate "
        "--quantMode TranscriptomeSAM GeneCounts"

rule featureCounts:
    input:
        bam="rna/{sample}/{sample}.Aligned.sortedByCoord.out.bam",
        gtf="/storage/zhangyanxiaoLab/xiongxiong/Reference/Annotation/mm10.gencode.annotation.gtf"
    output:
        counts="rna/{sample}/{sample}_counts.txt"
    shell:
        "featureCounts -T 40 -a {input.gtf} -t exon -g gene_id --extraAttributes gene_name -o {output.counts} {input.bam} "

rule combine_counts:
    input:
        expand("rna/{sample}/{sample}_counts.txt", sample=SAMPLES)
    output:
        "rna/bulk_rna_counts.txt"
    shell:
        '''
        cut -f 1,7 {input[0]} |grep -v '^#' > {output}
        for file in {input}; do
            cut -f 8 $file |grep -v '^#' | paste {output} - > tmp && mv tmp {output}
        done
        '''

rule star_align_te:
    input:
        r1_clean="fastq/{sample}/{sample}_R1.clean.fastq.gz",
        r2_clean="fastq/{sample}/{sample}_R2.clean.fastq.gz"
    output:
        bam="te/{sample}/{sample}.Aligned.sortedByCoord.out.bam"
    shell:
        "STAR --runThreadN 8 "
        "--genomeDir /storage/zhangyanxiaoLab/share/STAR_index/mm10 "
        "--readFilesIn {input.r1_clean} {input.r2_clean} "
        "--readFilesCommand zcat "
        "--outFileNamePrefix te/{wildcards.sample}/{wildcards.sample}. "
        "--outSAMtype BAM SortedByCoordinate "
        "--quantMode TranscriptomeSAM GeneCounts "
        "--winAnchorMultimapNmax 100 "
        "--outFilterMultimapNmax 100"
   
   
rule tetranscripts:
    input:
        "te/{sample}/{sample}.Aligned.sortedByCoord.out.bam"
    output:
        "te/{sample}/{sample}.cntTable"
    conda:
        "TEtranscripts"
    shell:
        '''
        TEcount --sortByPos --format BAM --mode multi  --project {wildcards.sample} --outdir te/{wildcards.sample} \
        --GTF /storage/zhangyanxiaoLab/xiongxiong/Reference/Annotation/mm10.gencode.annotation.gtf \
        --TE /storage/zhangyanxiaoLab/qihongjian/datasets/mm10_rmsk_TE.gtf \
        -b {input}  
        '''
        
rule combine_te_counts:
    input:
        expand("te/{sample}/{sample}.cntTable", sample=SAMPLES)
    output:
        "te/te_counts.txt"
    shell:
        '''
        cut -f 1 {input[0]} > {output}
        for file in {input}; do
            cut -f 2 $file | paste {output} - > tmpfile && mv tmpfile {output}
        done
        '''



