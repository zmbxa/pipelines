### This snakemake file can process paired-end bulk RNAseq data (Zokor) ###
### The workdir need to be changed as data containing folder (include fastq/) and all .fq file in fastq/[name]_1/2.fastq.gz ###
### !!! Change SAMPLES to a Tuple contain sample names (for snakemake final output)!!! ###
### workflow:  *** fastp ---> STAR ---> featureCounts ---> TECounts ***  output fetureCount and TECount matrix  ###
### Usage:  snakemake --cores 4 -s ~/pipelines/RNAseq/mapping/bulk_PE_RNA.reminder.smk --conda-frontend conda

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
        expand("rna/{sample}/{sample}.Log.final.out",sample=SAMPLES),
        expand("rna/{sample}/{sample}.nodup.srt.bam",sample=SAMPLES),
        expand("qc/{sample}/{sample}.dup.qc",sample=SAMPLES),
        "all_sample.qc.txt",
        #"te/te_counts.txt",

rule fastp:
    input:
        r1="fastq/{sample}_1.fastq.gz",
        r2="fastq/{sample}_2.fastq.gz",
    output:
        r1_clean=temp("fastq/{sample}/{sample}_R1.clean.fastq.gz"),
        r2_clean=temp("fastq/{sample}/{sample}_R2.clean.fastq.gz"),
        report_json="fastq/{sample}/{sample}.fastp.json",
        report_html="fastq/{sample}/{sample}.fastp.html"
    shell:
        "fastp -i {input.r1} -I {input.r2} -o {output.r1_clean} -O {output.r2_clean} -j {output.report_json} -h {output.report_html}"


rule star_align:
    input:
        r1_clean="fastq/{sample}/{sample}_R1.clean.fastq.gz",
        r2_clean="fastq/{sample}/{sample}_R2.clean.fastq.gz"
    output:
        bam=temp("rna/{sample}/{sample}.Aligned.sortedByCoord.out.bam"),
        raw_qc="rna/{sample}/{sample}.raw.flagstat.qc",
        star_out = "rna/{sample}/{sample}.Log.final.out",
        temp=temp("rna/{sample}/{sample}.Aligned.toTranscriptome.out.bam")
    shell:
        '''
        STAR --runThreadN 8 \
        --genomeDir  /storage/zhangyanxiaoLab/niuyuxiao/annotations/STAR_index/Zokor2 \
        --readFilesIn {input.r1_clean} {input.r2_clean} \
        --readFilesCommand zcat \
        --outFileNamePrefix rna/{wildcards.sample}/{wildcards.sample}. \
        --outSAMtype BAM SortedByCoordinate \
        --quantMode TranscriptomeSAM GeneCounts
        sambamba flagstat -t 6 {output.bam} |&tee {output.raw_qc}
        '''

rule bam_rmdup:
    input:
        bam = "rna/{sample}/{sample}.Aligned.sortedByCoord.out.bam",
    output:
        bam="rna/{sample}/{sample}.nodup.srt.bam",
        qc="qc/{sample}/{sample}.dup.qc",
        log="qc/{sample}/{sample}.rmdup.log"
    shell:
        '''
        sambamba markdup -t 6 -r {input.bam} {output.bam}  |&tee {output.log}
        sambamba flagstat -t 6 {output.bam} |&tee {output.qc}
        '''

rule featureCounts:
    input:
        bam="rna/{sample}/{sample}.nodup.srt.bam",
        gtf="/storage/zhangyanxiaoLab/niuyuxiao/annotations/gtf/Zokor2.gtf"
    output:
        counts="rna/{sample}/{sample}_counts.txt"
    shell:
        "featureCounts -T 40 -a {input.gtf} -t exon -g gene_id --extraAttributes gene_name -o {output.counts} {input.bam} "

rule lastqc:
    input:
        raw = expand("rna/{sample}/{sample}.Log.final.out",sample=SAMPLES),
        dup = expand("qc/{sample}/{sample}.rmdup.log",sample=SAMPLES)
    output:
        "{sample}.qc.txt"
    threads: 1
    run:
      out = open(output[0],'w')
      out.write("\t".join(["Sample","Total","Mapped","Filtered","Uniq","Map%","Dup%"])+"\n")
      for idx in range(len(input.raw)):
        samples = re.match("rna\/(.*)\/Log.final.out",input.raw[idx]).groups()[0]
        raw_file = open(input.raw[idx], 'r')
        for line in raw_file:
          words = line.strip().split('|')
          print(words)
          if words[0] == "Number of input reads ":
              total = words[1]
          elif words[0] == "Uniquely mapped reads number ": 
              uniq_mapped = words[1]
          elif words[0] == "Number of reads mapped to multiple loci ":
              multi_mapped = words[1]
        mapped = str(int(uniq_mapped)+int(multi_mapped))
        raw_file.close()
        dup_file = open(input.dup[idx], 'r')
        for line in dup_file:
          if "sorted" in line:
            filt = line.split()[1]
          elif "found" in line:
            duplicates = line.split()[1] // 2
        dup_file.close()
        nodup = str(int(filt)-int(duplicates))
        map_p = "%.2f"%(float(mapped)/float(total))
        dup_p = "%.2f"%(float(duplicates)/float(filt))
        out.write("\t".join([samples,total,mapped,filt,nodup,map_p,dup_p])+"\n")

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
        python ~/pipelines/reminder.py "Your bulk RNA mapping are successfully done! Congrats~" "Need a cup of Latte?"
        '''


### filt QC extract



