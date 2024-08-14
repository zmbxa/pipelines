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
spikeGENOME = config["spikeGENOME"]

BOWTIE_DICT = {"hg38":"/storage/zhangyanxiaoLab/share/bowtie2_index/hg38",
            "mm10":"/storage/zhangyanxiaoLab/share/bowtie2_index/mm10",
            "mm10_with_lambda":"/storage/zhangyanxiaoLab/niuyuxiao/annotations/bowtie2_index/mm10_lambda/mm10_lambda",
            "mm10_GFP":"/storage/zhangyanxiaoLab/share/bowtie2_index/GFP_AAVCre_index/mm10_gfp_AAVcre",
            "dm6":"/storage/zhangyanxiaoLab/share/bowtie2_index/dm6",
            "ce11":"/storage/zhangyanxiaoLab/share/bowtie2_index/ce11",
            "lambda":"/storage/zhangyanxiaoLab/niuyuxiao/annotations/bowtie2_index/lambda/lambda"
            }


#BOWTIE_DICT = {
#            "hg38_with_dm6":"/storage/zhangyanxiaoLab/niuyuxiao/annotations/bowtie2_index/hg38_with_dm6/hg38_with_dm6",
#            "mm10_with_lambda":"/storage/zhangyanxiaoLab/niuyuxiao/annotations/bowtie2_index/mm10_lambda/mm10_lambda",
#            "dm6":"/storage/zhangyanxiaoLab/share/bowtie2_index/dm6"
#            }

bt_index_main = BOWTIE_DICT[GENOME]
bt_index_spike = BOWTIE_DICT[spikeGENOME]

rule all:
  input:
    expand("bam/{sample}.filt.srt.bam",sample=SAMPLES),
    expand("bigWig/{sample}.scaled.filt.srt.bw",sample=SAMPLES),
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

rule bowtie2_align_main:
  input:
    R1="fastq_trim/{sample}_val_1.fq.gz",
    R2="fastq_trim/{sample}_val_2.fq.gz",
  output:
    bam=temp("bam/{sample}.bam"),
    raw_qc="qc/{sample}.raw.flagstat.qc"
  params:
    bowtie2_index=bt_index_main
  shell:
    '''
    BOWTIE2_REF={params.bowtie2_index}

    bowtie2 -x $BOWTIE2_REF \
    -1 {input.R1}  -2 {input.R2} \
    -X 2000 -p 5  | samtools view -@ 8 -bh > {output.bam}
    sambamba flagstat -t 6 {output.bam} |&tee {output.raw_qc}
     '''
     
rule unmap_reads_extract:
  input:
    bam=temp("bam/{sample}.bam"),
  output:
    um_bam=temp("unmapped/bam/{sample}_unmap.bam"),
    umsrt_bam="unmapped/bam/{sample}_unmap.srt.bam",
    um_fq1="unmapped/fastq/{sample}_unmap_1.fastq",
    um_fq2="unmapped/fastq/{sample}_unmap_2.fastq",
  #param:
  shell:
    '''
    sambamba view -t 6 --num-filter="13"  {input.bam} -f bam -h -o {output.um_bam}
    sambamba sort -o {output.umsrt_bam} -N -t 6 {output.um_bam}
    bedtools bamtofastq -i {output.umsrt_bam} -fq {output.um_fq1} -fq2 {output.um_fq2}

    '''
  
     
rule bowtie2_align_spikeIn:
  input:
    R1="unmapped/fastq/{sample}_unmap_1.fastq",
    R2="unmapped/fastq/{sample}_unmap_2.fastq",
  output:
    bam=temp("bam/{sample}_spikein.bam"),
    raw_qc="qc/{sample}_spikein.raw.flagstat.qc"
  params:
    bowtie2_index=bt_index_spike
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
    sorted="bam/{sample}.sorted.bam",
    log = "logs/{sample}.markDup.log",
  shell:
    '''
    sambamba sort --tmpdir=$PWD -t 10 -m 4G {input} {output.sorted}
    sambamba markdup -p -r -t 6 --overflow-list-size 600000 {output.sorted} {output.nodup} |&tee {output.log}
    '''
rule bam_filt:
  input:
    bam = "bam/{sample}.sorted.bam",
  output:
    bam = "bam/{sample}.filt.srt.bam",
    bai = "bam/{sample}.filt.srt.bam.bai",
    dup = "qc/{sample}.dup.flagstat.qc"
  shell:
    "samtools view -F 1804 -b {input.bam} > {output.bam};"
    "sambamba index -t 4 {output.bam};"
    "sambamba flagstat -t 6 {output.bam} |&tee {output.dup};"

rule tsse:
  input:
    expand("bam/{sample}.filt.nodup.srt.bam",sample=SAMPLES),
  output:
    "tsse/tsse.csv"
  conda:"snapATAC2"
  shell:
    '''
    python /storage/zhangyanxiaoLab/niuyuxiao/pipelines/miniscripts/tsse_snapATAC2.py -b bam/*.filt.nodup.srt.bam
    '''

#rule tsse:
#  input:
#    expand("bam/{sample}.filt.nodup.srt.bam",sample=SAMPLES),
#  output:
#    "tsse/tsse.csv"
#  conda:"snapATAC2"
#  shell:
#    '''
#    python /storage/zhangyanxiaoLab/niuyuxiao/pipelines/miniscripts/tsse_snapATAC2.py -b bam/*.filt.nodup.srt.bam    
#    '''

rule last_qc:
  input:
    raw = expand("qc/{sample}.raw.flagstat.qc",sample=SAMPLES),
    spike_raw = expand("qc/{sample}_spikein.raw.flagstat.qc",sample=SAMPLES),
    dup = expand("logs/{sample}.markDup.log",sample=SAMPLES),
    trim = expand("logs/{sample}.trim.log",sample=SAMPLES),
  output:
    "{sample}.qc.txt"
  threads:1
  run:
    out = open(output[0],'w')
    out.write("\t".join(["Sample","Total","Trim","Mapped","Filtered","Uniq","MappedToSpikeIn","ScaleFactor","Total_map%","Adapter%","Map%","Dup%","SpikeIn%"])+"\n")
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
      spike_raw_file = open(input.spike_raw[idx], 'r')
      for line in spike_raw_file:
        words = line.strip().split(' ')
        if len(words)>5:
          if words[3] == "properly":
            mapped2spike = str(int(int(words[0])/2))
            scaleFactor = float(10000)/float(mapped2spike)
      spike_raw_file.close()
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
      map_p = "%.4f"%(float(mapped)/float(trim))
      total_map_p = "%.4f"%((float(mapped)+float(mapped2spike))/float(trim))
      spikein_ratio = str(float(mapped2spike)/(float(mapped)+float(mapped2spike)))
      dup_p = "%.2f"%(float(duplicates)/float(mapped))
      out.write("\t".join([samples,total,trim,mapped,filt,nodup,mapped2spike,str(scaleFactor),total_map_p,adapter_p,map_p,dup_p,spikein_ratio])+"\n")

rule bam2bw:
  input:
    bam = "bam/{sample}.filt.srt.bam",
    qc = "all_sample.qc.txt"
  output:
    "bigWig/{sample}.scaled.filt.srt.bw"
  shell:
    '''
    sf=$(grep {wildcards.sample} {input.qc} |cut -f 8)
    source /storage/zhangyanxiaoLab/share/Pipelines/environments/python3env/bin/activate
    bamCoverage --scaleFactor ${{sf}} -b {input.bam} -o {output} --outFileFormat bigwig -bs 50 --numberOfProcessors 6 --normalizeUsing RPKM
    
    '''






