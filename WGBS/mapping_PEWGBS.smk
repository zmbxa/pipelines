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
print("Found samples:")
print(SAMPLES)


rule all:
  input:
    expand("bam/{sample}_bismark_bt2_pe.bam",sample=SAMPLES),
    expand("bam/{sample}_bt2_pe.deduplicated.bam",sample=SAMPLES),
    expand("qc/{sample}_bismark_bt2_PE_report.html",sample=SAMPLES),
    expand("bam/{sample}_bt2_pe.deduplication_report.txt",sample=SAMPLES),
    expand("meth_ext/{sample}_bt2_pe.deduplicated.bismark.cov.gz",sample=SAMPLES),
    expand("meth_ext/{sample}_bt2_pe.deduplicated.CpG_report.txt",sample=SAMPLES),
    expand("methylKit/{sample}_bt2_pe.deduplicated.CpG_report.txt_CG_methylkit.txt",sample=SAMPLES)

rule fastqc:
  input:
    r1="fastq/{sample}_1.fastq.gz",
    r2="fastq/{sample)_2.fastq.gz"
  output:
    "qc/{sample}_1_fastqc.html",
    "qc/{sample}_2_fastqc.html"
  shell:
    "fastqc -o qc {input.r1} {input.r2} -t 4 "

rule Trim_adapter:
  input:
    R1="fastq/{sample}_1.fastq.gz",
    R2="fastq/{sample}_2.fastq.gz",
  output:
    temp("fastq_trim/{sample}_1_val_1.fq.gz"),
    temp("fastq_trim/{sample}_2_val_2.fq.gz"),
    "qc/{sample}_1_val_1_fastqc.html",
    "qc/{sample}_2_val_2_fastqc.html"
#    log="logs/{sample}.log"
  conda:"cutadapt"
  shell:
    '''
    trim_galore  -q 25 --phred33 --length 20 --stringency 3 --paired -o fastq_trim  {input.R1} {input.R2} --core 4 --gzip --path_to_cutadapt /storage/zhangyanxiaoLab/niuyuxiao/anaconda3/envs/cutadapt/bin/cutadapt --fastqc_args "-o qc -t 4 "
    '''

rule Bismark_mapping:
  input:
    R1="fastq_trim/{sample}_1_val_1.fq.gz",
    R2="fastq_trim/{sample}_2_val_2.fq.gz",
  output:
    bam="bismark_out/{sample}_bismark_bt2_pe.bam",
    mapRepo="reports/{sample}_mapping_report.txt",
    #bismarkQC="reports/{sample}_mapping_report.html"
  shell:
    '''
    bismark --genome /storage/zhangyanxiaoLab/niuyuxiao/annotations/Bismark_ref/hg38 --parallel 8 \
     --gzip --un --ambiguous -q -L 20 -D 20 -R 3 --score_min L,0,-0.4 -1 {input.R1} -2 {input.R2} -o bismark_out 
    mv bismark_out/{wildcards.sample}_1_val_1_bismark_bt2_PE_report.txt reports/{wildcards.sample}_mapping_report.txt
    rename _1_val_1_bismark_bt2 _bismark_bt2 bismark_out/{wildcards.sample}* 
    #bismark2report --alignment_report {output.mapRepo} --dir reports/ 
    '''

rule rmDup:
  input:
    bam="bismark_out/{sample}_bismark_bt2_pe.bam",
    #map_repo="reports/{sample}_bismark_bt2_PE_report.txt",
  output:
    bam="bismark_out/{sample}.deduplicated.bam",
    repo="reports/{sample}.deduplication_report.txt",
    #bismarkQC="qc/{sample}_bt2_pe.deduplication_report.html"
  shell:
    '''
    deduplicate_bismark -p --bam {input.bam} --output_dir bismark_out -o {wildcards.sample}.bam
    #bismark2report --alignment_report {input.map_repo} --dedup_report {output.repo} -o {wildcards.sample}_bt2_pe.deduplication_report.html --dir qc
    mv bismark_out/{wildcards.sample}.deduplication_report.txt {output.repo}
    '''

rule methy_ext:
  input:
    ref="/storage/zhangyanxiaoLab/niuyuxiao/annotations/Bismark_ref/hg38",
    bam="bam/{sample}_bt2_pe.deduplicated.bam",
  output:
    bdg="meth_ext/{sample}_bt2_pe.deduplicated.bedGraph.gz",
    cov="meth_ext/{sample}_bt2_pe.deduplicated.bismark.cov.gz",
    cpgRepo="meth_ext/{sample}_bt2_pe.deduplicated.CpG_report.txt"
  shell:
    '''
    bismark_methylation_extractor -p --no_overlap --parallel 30 --bedGraph --buffer_size 90G --cytosine_report --output meth_ext --genome_folder {input.ref} {input.bam}
    '''
rule bismark2methylkit:
  input:
    cpgRepo="meth_ext/{sample}_bt2_pe.deduplicated.CpG_report.txt"
  output:
    meKit="methylKit/{sample}_bt2_pe.deduplicated.CpG_report.txt_CG_methylkit.txt"
  shell:
    '''
    cd meth_ext
    perl ~/pipelines/WGBS/bismark2methylkit.pl $(basename {input.cpgRepo})
    mv $(basename {input.cpgRepo})_CG_methylkit.txt ../methylKit/
    cd ..
    '''


rule reminder:
    input:
        "merge/DNA_matrix/matrix.mtx.gz",
        "merge/RNA_matrix/matrix.mtx.gz",
    output:
        "completeEmail.log"
    shell:
        '''
        python ~/pipelines/reminder.py "Your WGBS mapping pipeline are successfully done! Congrats~" "Need a cup of Latte?" > {output}
        '''
        














