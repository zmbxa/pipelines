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
    expand("logs/{sample}.adapterTrim.log",sample=SAMPLES),
    expand("logs/{sample}.bowtie2.log",sample=SAMPLES),
    expand("qc/{sample}.raw.flagstat.qc",sample=SAMPLES),
    expand("qc/{sample}.markdup.qc",sample=SAMPLES),
    expand("qc/{sample}.dup.flagstat.qc",sample=SAMPLES),
    "all_sample.qc.txt"

rule Trim_adapter:
  input:
    R1="fastq/{sample}_1.fastq.gz",
    R2="fastq/{sample}_2.fastq.gz",
  output:
    temp("fastq_trim/{sample}_val_1.fq.gz"),
    temp("fastq_trim/{sample}_val_2.fq.gz"),
    log="logs/{sample}.adapterTrim.log"
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
    bam=temp("bam/{sample}.bam"),
    log="logs/{sample}.bowtie2.log",
    raw_qc="qc/{sample}.raw.flagstat.qc",
  shell:
    '''
    BOWTIE2_REF=/storage/zhangyanxiaoLab/share/bowtie2_index/hg38

    bowtie2 -x $BOWTIE2_REF \
    -1 {input.R1}  -2 {input.R2} \
    -X 2000 -p 5  2>{output.log} | samtools view -@ 8 -bh > {output.bam}
    samtools flagstat {output.bam} > {output.raw_qc}
    '''

rule bam_rmDup:
  input:
    "bam/{sample}.bam"
  output:
    nodup=temp("bam/{sample}.nodup.bam"),
    nodupbai=temp("bam/{sample}.nodup.bam.bai"),
    sorted=temp("bam/{sample}.sorted.bam"),
    bai=temp("bam/{sample}.sorted.bam.bai"),
    log="qc/{sample}.markdup.qc"
  shell:
    '''
    sambamba sort --tmpdir=$PWD -t 10 -m 4G {input} {output.sorted}
    sambamba markdup -r -t 6 --overflow-list-size 600000 {output.sorted} {output.nodup} 2>&1 |tee  {output.log}
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
    
rule lastqc:
  input:
    raw=expand("qc/{sample}.raw.flagstat.qc",sample=SAMPLES),
    dup=expand("qc/{sample}.markdup.qc",sample=SAMPLES),
    trim=expand("logs/{sample}.adapterTrim.log",sample=SAMPLES),
  output:
    "{sample}.qc.txt"
  run:
      out = open(output[0],'w')
      out.write("\t".join(["Sample","Total","Trim","Mapped","Filtered","Uniq","Adapter%","Map%","Dup%"])+"\n")
      for idx in range(len(input.raw)):
        samples = re.match(r"qc\/(.*).raw.flagstat.qc",input.raw[idx]).groups()[0]
        raw_file = open(input.raw[idx], 'r')
        for line in raw_file:
          words = line.strip().split(' ')
          if len(words)>5:
            if words[4] == "total":
              trim = str(int(int(words[0])/2))
            elif words[3] == "properly":
              mapped = str(int(int(words[0])/2))
        raw_file.close()
        dup_file = open(input.dup[idx],'r')
        for line in dup_file:
          words = line.strip().split(' ')
          if words[0] == "sorted":
            filt = words[1]
          elif words[0] == "found":
            duplicates = words[1]
            nodup = str(int(filt)-int(duplicates))
        trim_file = open(input.trim[idx],'r')
        for line in trim_file:
          words = line.strip().split()
          if words[:3] == ['Total', 'reads', 'processed:']:
            total= words[3].replace(',','')
          if words[:3] == ['Reads', 'with', 'adapters:']:
            adapter_p = "%.2f"%(float(words[4][1:][:-2])/100)
         
        map_p = "%.2f"%(float(mapped)/float(trim))
        dup_p = "%.2f"%(float(duplicates)/float(mapped))
        out.write("\t".join([samples,total,trim,mapped,filt,nodup,adapter_p,map_p,dup_p])+"\n")










