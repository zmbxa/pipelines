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

BOWTIE_DICT = {"hg38":"/storage/zhangyanxiaoLab/share/bowtie2_index/hg38",
            "mm10":"/storage/zhangyanxiaoLab/share/bowtie2_index/mm10",
            "mm10_GFP":"/storage/zhangyanxiaoLab/share/bowtie2_index/GFP_AAVCre_index/mm10_gfp_AAVcre",
            "dm6":"/storage/zhangyanxiaoLab/share/bowtie2_index/dm6"
            }

bt_index = BOWTIE_DICT[GENOME]

rule all:
  input:
    expand("bam/{sample}.filt.nodup.srt.bam",sample=SAMPLES),
    expand("bigWig/{sample}.filt.nodup.srt.bw",sample=SAMPLES),
    "all_sample.qc.txt",


rule Trim_adapter:
  input:
    R1="fastq/{sample}_1.fastq.gz",
    R2="fastq/{sample}_2.fastq.gz",
  output:
    temp("fastq_trim/{sample}_val_1.fq.gz"),
    temp("fastq_trim/{sample}_val_2.fq.gz"),
    log="logs/{sample}.trim.log"
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
    raw_qc="qc/{sample}.raw.flagstat.qc"
  params:
    bowtie2_index=bt_index
  shell:
    '''
    BOWTIE2_REF={params.bowtie2_index}

    bowtie2 -x $BOWTIE2_REF \
    -1 {input.R1}  -2 {input.R2} \
    -X 2000 -p 5  | samtools view -@ 8 -bh > {output.bam}
    sambamba flagstat -t 6 {output.bam} |&tee {output.raw_qc}
     '''

rule bam_rmDup:
  input:
    "bam/{sample}.bam"
  output:
    nodup=temp("bam/{sample}.nodup.bam"),
    nodupbai=temp("bam/{sample}.nodup.bam.bai"),
    sorted=temp("bam/{sample}.sorted.bam"),
    log = "logs/{sample}.markDup.log",
  shell:
    '''
    sambamba sort --tmpdir=$PWD -t 10 -m 4G {input} {output.sorted}
    sambamba markdup -p -r -t 6 --overflow-list-size 600000 {output.sorted} {output.nodup} |&tee {output.log}
    '''
rule bam_filt:
  input:
    bam = "bam/{sample}.nodup.bam",
  output:
    bam = "bam/{sample}.filt.nodup.srt.bam",
    bai = "bam/{sample}.filt.nodup.srt.bam.bai",
    dup = "qc/{sample}.dup.flagstat.qc"
  shell:
    "samtools view -F 1804 -b {input.bam} > {output.bam};"
    "sambamba index -t 4 {output.bam};"
    "sambamba flagstat -t 6 {output.bam} |&tee {output.dup};"

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

rule last_qc:
  input:
    raw = expand("qc/{sample}.raw.flagstat.qc",sample=SAMPLES),
    dup = expand("logs/{sample}.markDup.log",sample=SAMPLES),
    trim = expand("logs/{sample}.trim.log",sample=SAMPLES)
  output:
    "{sample}.qc.txt"
  threads:1
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
      dup_file = open(input.dup[idx], 'r')
      for line in dup_file:
        if "sorted" in line:
          filt = line.split()[1]
        elif "found" in line:
          duplicates = int(line.split()[1]) / 2
      dup_file.close()
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
      os.system('python /storage/zhangyanxiaoLab/niuyuxiao/pipelines/reminder.py "Your bulk ATAC mapping are successfully done! Congrats~" "Need a cup of Latte?"')








