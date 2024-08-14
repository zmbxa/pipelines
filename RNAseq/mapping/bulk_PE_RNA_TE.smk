### This snakemake file can process paired-end bulk RNAseq data (mouse mm10) ###
### The workdir need to be changed as data containing folder (include fastq/) and all .fq file in fastq/[name]/[name]_1/2.fq.gz ###
### !!! Change SAMPLES to a Tuple contain sample names (for snakemake final output)!!! ###
### workflow:  *** fastp ---> STAR ---> featureCounts ---> TECounts ***  output fetureCount and TECount matrix  ###
### Usage:  snakemake --cores 4 -s bulk_PE_RNA.smk --use-conda --conda-frontend conda

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
            "hg38_with_dm6":"/storage/zhangyanxiaoLab/niuyuxiao/annotations/gtf/hg38_dm6in.gtf",
            #"hg38_with_dm6":"/storage/zhangyanxiaoLab/niuyuxiao/annotations/gtf/hg38_dm6in1.gtf",
            "mm10":"/storage/zhangyanxiaoLab/xiongxiong/Reference/Annotation/mm10.gencode.annotation.gtf",
            "mm10_with_dm6in1":"/storage/zhangyanxiaoLab/niuyuxiao/annotations/gtf/mm10_dm6in1.gtf", 
            "mm10_with_dm6":"/storage/zhangyanxiaoLab/niuyuxiao/annotations/gtf/mm10_dm6in.gtf",
            "mm10_GFP":"/storage/zhangyanxiaoLab/niuyuxiao/annotations/gtf/mm10_gfp_AAVCre.gtf",
            "dm6":"/storage/zhangyanxiaoLab/xiongxiong/Reference/Annotation/dm6.gencode.annotation.gtf",
            "Zokor2":"/storage/zhangyanxiaoLab/niuyuxiao/annotations/gtf/Zokor2.gtf",
            "Zokor3":"/storage/zhangyanxiaoLab/niuyuxiao/annotations/gtf/EBaileyi.gtf"
            }
STAR_DICT = {"hg38":"/storage/zhangyanxiaoLab/share/STAR_index/hg38",
            "hg38_with_dm6":"/storage/zhangyanxiaoLab/niuyuxiao/annotations/STAR_index/hg38_with_dm6",
            #"hg38_with_dm6":"/storage/zhangyanxiaoLab/niuyuxiao/annotations/STAR_index/hg38_with_dm6in1",
            "dm6":"/storage/zhangyanxiaoLab/share/STAR_index/dm6",
            "mm10":"/storage/zhangyanxiaoLab/share/STAR_index/mm10",
            "mm10_with_dm6in1":"/storage/zhangyanxiaoLab/niuyuxiao/annotations/STAR_index/mm10_with_dm6",
            "mm10_with_dm6":"/storage/zhangyanxiaoLab/niuyuxiao/annotations/STAR_index/mm10_with_dm6",
            "mm10_GFP":"/storage/zhangyanxiaoLab/niuyuxiao/annotations/STAR_index/mm10_gfp_AAVCre",
            "Zokor2":"/storage/zhangyanxiaoLab/niuyuxiao/annotations/STAR_index/Zokor2",
            "Zokor3":"/storage/zhangyanxiaoLab/niuyuxiao/annotations/STAR_index/EBaileyi"
            }
TE_GTFDICT = {"hg38":"/storage/zhangyanxiaoLab/niuyuxiao/annotations/hg38_rmsk_TE.gtf",
            "hg38_with_dm6":"/storage/zhangyanxiaoLab/niuyuxiao/annotations/hg38_rmsk_TE.gtf",
            "dm6":"/storage/zhangyanxiaoLab/niuyuxiao/annotations/gtf/dm6_rmsk.gtf",
            "mm10":"/storage/zhangyanxiaoLab/niuyuxiao/annotations/gtf/mm10_rmsk_TE.gtf",
            "mm10_with_dm6in1":"/storage/zhangyanxiaoLab/niuyuxiao/annotations/gtf/mm10_rmsk_TE.gtf",
            "mm10_with_dm6":"/storage/zhangyanxiaoLab/niuyuxiao/annotations/gtf/mm10_rmsk_TE.gtf",
            "mm10_GFP":"/storage/zhangyanxiaoLab/niuyuxiao/annotations/gtf/mm10_rmsk_TE.gtf",
            "Zokor2":"/storage/zhangyanxiaoLab/niuyuxiao/annotations/gtf/EBaileyi_v3_rmsk.TE.gtf",
            "Zokor3":"/storage/zhangyanxiaoLab/niuyuxiao/projects/non_Model_TE_Zokor/analysis/repeatMasker/fromZokor/EBaileyi.rmsk3.gtf"
            }


GTF = GTF_DICT[GENOME]
STAR_IND = STAR_DICT[GENOME]
TEGTF = TE_GTFDICT[GENOME]

rule all:
    input:
        "rna/bulk_rna_counts.txt",
        expand("rna/{sample}/{sample}.Aligned.sortedByCoord.out.bam",sample=SAMPLES),
        "all_sample.qc.txt",
        "te/te_counts.txt",
#        "test.txt"

rule fastp:
    input:
        r1="fastq/{sample}_1.fastq.gz",
        r2="fastq/{sample}_2.fastq.gz",
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
        bam="rna/{sample}/{sample}.Aligned.sortedByCoord.out.bam",
        raw_qc="rna/{sample}/{sample}.raw.flagstat.qc",
        star_out = "rna/{sample}/{sample}.Log.final.out",
        temp=temp("rna/{sample}/{sample}.Aligned.toTranscriptome.out.bam")
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
        bam="rna/{sample}/{sample}.Aligned.sortedByCoord.out.bam",
    params:
        gtf=GTF
    output:
        counts="rna/{sample}/{sample}_counts.txt"
    shell:
        "featureCounts -T 40 -a {params.gtf} -t exon -g gene_id --extraAttributes gene_name -o {output.counts} {input.bam} "

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
        bam="te/{sample}/{sample}.Aligned.sortedByCoord.out.bam",
        te_qc="te/{sample}/{sample}_TEstar.flagstat.qc",
        STAR_out_te="te/{sample}/{sample}.Log.final.out"
    params:
        star_index=STAR_IND
    shell:
        '''
        STAR --runThreadN 8 \
        --genomeDir {params.star_index} \
        --readFilesIn {input.r1_clean} {input.r2_clean} \
        --readFilesCommand zcat \
        --outFileNamePrefix te/{wildcards.sample}/{wildcards.sample}. \
        --outSAMtype BAM SortedByCoordinate \
        --quantMode TranscriptomeSAM GeneCounts \
        --winAnchorMultimapNmax 100 \
        --outFilterMultimapNmax 100
        sambamba flagstat -t 6 {output.bam} |&tee {output.te_qc}
       '''
   
rule tetranscripts:
    input:
        "te/{sample}/{sample}.Aligned.sortedByCoord.out.bam"
    output:
        "te/{sample}/{sample}.cntTable"
    params:
        geneGtf=GTF,
        TEgtf=TEGTF
    conda:
        "TEtranscripts"
    shell:
        '''
        TEcount --sortByPos --format BAM --mode multi  --project {wildcards.sample} --outdir te/{wildcards.sample} \
        --GTF {params.geneGtf} \
        --TE {params.TEgtf} \
        -b {input}  
        '''

rule lastqc:
    input:
        raw = expand("rna/{sample}/{sample}.Log.final.out",sample=SAMPLES),
        te = expand("te/{sample}/{sample}.Log.final.out",sample=SAMPLES)
    output:
        "{sample}.qc.txt"
    threads: 1
    run:
      out = open(output[0],'w')
      out.write("\t".join(["Sample","Total","Mapped","Map%","Total_withTE","Mapped_withTE","withTE_Map%"])+"\n")
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
        TE_file = open(input.te[idx], 'r')
        for line in TE_file:
          words = line.strip().split('|')
          #print(words)
          if words[0] == "Number of input reads ":
              TEtotal = words[1]
          elif words[0] == "Uniquely mapped reads number ": 
              TEuniq_mapped = words[1]
          elif words[0] == "Number of reads mapped to multiple loci ":
              TEmulti_mapped = words[1]
        mappedTE = str(int(TEuniq_mapped)+int(TEmulti_mapped))
        TE_file.close()
        map_p = "%.2f"%(float(mapped)/float(total))
        TEmap_p = "%.2f"%(float(mappedTE)/float(TEtotal))
        out.write("\t".join([samples,total,mapped,map_p,TEtotal,mappedTE,TEmap_p])+"\n")
        
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
        python ~/pipelines/reminder.py "Your bulk RNA mapping are successfully done! Congrats~" "Need a cup of Latte?"
        '''



