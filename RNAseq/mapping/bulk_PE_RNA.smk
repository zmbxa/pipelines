### This snakemake file can process paired-end bulk RNAseq data  ###
### The workdir need to be changed as data containing folder (include fastq/) and all .fq file in fastq/[name]_1/2.fastq.gz ###
### Please SPECIFY species in command line use --config GENOME=mm10/hg38/...
### fastp --> STAR --> remove dup (sambamba) --> generate bigwig* --> featureCount -->lastQC
### Usage:  snakemake --cores 4 -s ~/pipelines/RNAseq/mapping/bulk_PE_RNA.smk --config GENOME=<mm10/hg38/Zokor2/...> --conda-frontend conda

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

GTF_DICT = {"hg38":"/storage/zhangyanxiaoLab/xiongxiong/Reference/Annotation/hg38.gencode.annotation.gtf",
            "hg38_with_dm6":"/storage/zhangyanxiaoLab/niuyuxiao/annotations/gtf/hg38_dm6in1.gtf",
            "mm10":"/storage/zhangyanxiaoLab/xiongxiong/Reference/Annotation/mm10.gencode.annotation.gtf",
            "mm10_GFP":"/storage/zhangyanxiaoLab/niuyuxiao/annotations/gtf/mm10_gfp_AAVCre.gtf",
            "dm6":"/storage/zhangyanxiaoLab/xiongxiong/Reference/Annotation/dm6.gencode.annotation.gtf",
            "Zokor2":"/storage/zhangyanxiaoLab/niuyuxiao/annotations/gtf/Zokor2.gtf",
            "Zokor3":"/storage/zhangyanxiaoLab/niuyuxiao/annotations/gtf/EBaileyi.gtf",
            "Sta_aureus":"/storage/zhangyanxiaoLab/niuyuxiao/annotations/gtf/Staphylococcus_aureus_gca_001018655.ASM101865v2.60.gtf"
            }
STAR_DICT = {"hg38":"/storage/zhangyanxiaoLab/share/STAR_index/hg38",
            "hg38_with_dm6":"/storage/zhangyanxiaoLab/niuyuxiao/annotations/STAR_index/hg38_with_dm6",
            "dm6":"/storage/zhangyanxiaoLab/share/STAR_index/dm6",
            "mm10":"/storage/zhangyanxiaoLab/share/STAR_index/mm10",
            "mm10_GFP":"/storage/zhangyanxiaoLab/niuyuxiao/annotations/STAR_index/mm10_gfp_AAVCre",
            "Zokor2":"/storage/zhangyanxiaoLab/niuyuxiao/annotations/STAR_index/Zokor2",
            "Zokor3":"/storage/zhangyanxiaoLab/niuyuxiao/annotations/STAR_index/EBaileyi",
            "Sta_aureus":"/storage/zhangyanxiaoLab/niuyuxiao/annotations/STAR_index/Sta_aureus"
            }

GTF = GTF_DICT[GENOME]
STAR_IND = STAR_DICT[GENOME]

rule all:
    input:
        "rna/bulk_rna_counts.txt",
        expand("rna/{sample}/{sample}.Log.final.out",sample=SAMPLES),
        expand("rna/{sample}/{sample}.Aligned.sortedByCoord.out.bam",sample=SAMPLES),
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
        bam="rna/{sample}/{sample}.Aligned.sortedByCoord.out.bam",
        raw_qc="rna/{sample}/{sample}.raw.flagstat.qc",
        star_out = "rna/{sample}/{sample}.Log.final.out",
        temp="rna/{sample}/{sample}.Aligned.toTranscriptome.out.bam"
    params:
        star_index=STAR_IND
    shell:
        '''
        STAR --runThreadN 8 \
        --genomeDir  {params.star_index} \
        --readFilesIn {input.r1_clean} {input.r2_clean} \
        --readFilesCommand zcat \
        --outFileNamePrefix rna/{wildcards.sample}/{wildcards.sample}. \
        --outSAMtype BAM SortedByCoordinate \
        --quantMode TranscriptomeSAM GeneCounts
        sambamba flagstat -t 6 {output.bam} |&tee {output.raw_qc}
        '''
rule featureCounts:
    input:
        bam = "rna/{sample}/{sample}.Aligned.sortedByCoord.out.bam",
    output:
        counts="rna/{sample}/{sample}_counts.txt"
    params:
        gtf=GTF
    shell:
        "featureCounts -T 40 -a {params.gtf} -t exon -g gene_id --extraAttributes gene_name -o {output.counts} {input.bam} "

rule lastqc:
    input:
        raw = expand("rna/{sample}/{sample}.Log.final.out",sample=SAMPLES),
    output:
        "{sample}.qc.txt"
    threads: 1
    run:
      out = open(output[0],'w')
      out.write("\t".join(["Sample","Total","Mapped","Map%"])+"\n")
      for idx in range(len(input.raw)):
        samples = re.match("rna\/(.*)\/(.*)Log.final.out",input.raw[idx]).groups()[0]
        raw_file = open(input.raw[idx], 'r')
        for line in raw_file:
          words = line.strip().split('|')
          #print(words)
          if words[0] == "Number of input reads ":
              total = words[1]
          elif words[0] == "Uniquely mapped reads number ": 
              uniq_mapped = words[1]
          elif words[0] == "Number of reads mapped to multiple loci ":
              multi_mapped = words[1]
        mapped = str(int(uniq_mapped)+int(multi_mapped))
        raw_file.close()
        map_p = "%.2f"%(float(mapped)/float(total))
        out.write("\t".join([samples,total,mapped,map_p])+"\n")

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
        python ~/pipelines/reminder.py "Your bulk RNA mapping are successfully done! Congrats~" "$(cat all_sample.qc.txt)"
        '''



