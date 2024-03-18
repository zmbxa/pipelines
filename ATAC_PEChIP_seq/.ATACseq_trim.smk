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
    expand("bam/{sample}.filt.nodup.srt.bam",sample=SAMPLES),
    expand("bigWig/{sample}.filt.nodup.srt.bw",sample=SAMPLES),

rule Trim_adapter:
  input:
    R1="fastq/{sample}_1.fastq.gz",
    R2="fastq/{sample}_2.fastq.gz",
  output:
    temp("fastq_trim/{sample}_val_1.fq.gz"),
    temp("fastq_trim/{sample}_val_2.fq.gz"),
    log="logs/{sample}.log"
  conda:"cutadapt"
  shell:
    '''
    trim_galore  --paired {input.R1} {input.R2} --core 4 --output_dir fastq_trim --basename {wildcards.sample} 2> {output.log}
    '''

rule bowtie2_align:
  input:
    R1="fastq_trim/{sample}_val_1.fq.gz",
    R2="fastq_trim/{sample}_val_2.fq.gz",
  output:
    temp("bam/{sample}.bam")
  shell:
    '''
    BOWTIE2_REF=/storage/zhangyanxiaoLab/share/bowtie2_index/mm10

    bowtie2 -x $BOWTIE2_REF \
    -1 {input.R1}  -2 {input.R2} \
    -X 2000 -p 5  | samtools view -@ 8 -bh > {output}
    '''

rule bam_rmDup:
  input:
    "bam/{sample}.bam"
  output:
    nodup=temp("bam/{sample}.nodup.bam"),
    nodupbai=temp("bam/{sample}.nodup.bam.bai"),
    sorted=temp("bam/{sample}.sorted.bam"),
    bai=temp("bam/{sample}.sorted.bam.bai")
  shell:
    '''
    sambamba sort --tmpdir=$PWD -t 10 -m 4G {input} {output.sorted}
    sambamba markdup -p -r -t 6 --overflow-list-size 600000 {output.sorted} {output.nodup}
    '''
rule bam_filt:
  input:
    bam = "bam/{sample}.nodup.bam",
  output:
    bam = "bam/{sample}.filt.nodup.srt.bam",
    bai = "bam/{sample}.filt.nodup.srt.bam.bai",
    stat = "qc/{sample}.dup.flagstat.qc"
  shell:
    "samtools view -F 1804 -b {input.bam} > {output.bam};"
    "sambamba index -t 4 {output.bam};"
    "samtools flagstat {output.bam} > {output.stat};"

rule bam2bw:
  input:
    "bam/{sample}.filt.nodup.srt.bam"
  output:
    "bigWig/{sample}.filt.nodup.srt.bw"
  shell:
    '''
    source /storage/zhangyanxiaoLab/share/Pipelines/environments/python3env/bin/activate
    bamCoverage -b {input} -o {output} --outFileFormat bigwig -bs 50 --numberOfProcessors 6 --normalizeUsing RPKM
    
    '''









