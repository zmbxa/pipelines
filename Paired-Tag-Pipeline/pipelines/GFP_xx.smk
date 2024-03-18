#snakefile 

import pandas as pd
import os

## setdb1 AAV-injected brain FC;Barcode: 01-06 control, 07-12 Cre
samples = pd.read_csv("/storage/zhangyanxiaoLab/niuyuxiao/projects/aging_SETDB1/data/paired_tag/setdb1_GFP/230528_GY_setdb1KO_GFP/sampleInfo.csv")
workdir:"/storage/zhangyanxiaoLab/niuyuxiao/projects/aging_SETDB1/data/paired_tag/setdb1_GFP/230528_GY_setdb1KO_GFP"

print(samples)

rule require_all:
    input:
        "bam_list.txt",
        "DNA_matrix/matrix.mtx.gz",
        "RNA_matrix/matrix.mtx.gz",
        expand("{sample}/{sample}/process/{sample}_CB.Counts",sample=samples["prefix"].unique()),
        expand("{sample}/{sample}/bam/{sample}_UsefulReads.bam",sample=samples["prefix"].unique())

        #expand("{sample}/matrix/matrix.mtx.gz", sample=samples["prefix"].unique()),
        #expand("qc/{dna}_{rna}_reads.pdf", zip, dna=samples[samples["library"]=="DNA"]["prefix"].unique(), rna=samples[samples["library"]=="RNA"]["prefix"].unique())

rule paired_tag_bam:
    input:
        R1 = lambda wildcards: samples[samples["prefix"]==wildcards.prefix]["R1"].values[0],
        R2 = lambda wildcards: samples[samples["prefix"]==wildcards.prefix]["R2"].values[0]
    output:
        "{prefix}/{prefix}/process/{prefix}_CB.Counts",
        "{prefix}/{prefix}/bam/{prefix}_UsefulReads.bam"
    params:
        library = lambda wildcards: samples[samples["prefix"]==wildcards.prefix]["library"].values[0]
    conda:
        "paired_tag"
    shell:
        """
        mkdir -p {wildcards.prefix}
        cd {wildcards.prefix}
        rm -rf *
        mkdir -p fastq
        ln -s {input.R1} fastq/{wildcards.prefix}_R1.fastq.gz
        ln -s {input.R2} fastq/{wildcards.prefix}_R2.fastq.gz
        bash ~/pipelines/Paired-Tag-Pipeline/pipelines/pipeline_PT_gfp.sh -g mm10 -l {params.library} -p 4
        """

rule list_make:
    output:
        "bam_list.txt"
    run:
        for row in samples.head(4).itertuples():
          f=open("bam_list.txt","a")
          f.write("/storage/zhangyanxiaoLab/niuyuxiao/projects/aging_SETDB1/data/paired_tag/setdb1_GFP/230528_GY_setdb1KO_GFP/"+getattr(row,'prefix')+"/"+getattr(row,'prefix')+"/bam/"+getattr(row,'prefix')+"_UsefulReads.bam"+",/storage/zhangyanxiaoLab/niuyuxiao/projects/aging_SETDB1/data/paired_tag/setdb1_GFP/230528_GY_setdb1KO_GFP/"+getattr(row,'exp_prefix')+"/"+getattr(row,'exp_prefix')+"/bam/"+getattr(row,'exp_prefix')+"_UsefulReads.bam")
          f.write('\n')
        f.close()

rule filt_merge:
    input:
        "bam_list.txt"
    output:
        "filtered_bam/merged_DNA.bam",
        "filtered_bam/merged_RNA.bam"
    shell:
        "bash /storage/zhangyanxiaoLab/xiongxiong/scripts/paired_tag/pipelines/filter_merge_bam.sh -i {input} -d 300 -r 500 "

rule bam2mtx:
    input:
        DNAbam="filtered_bam/merged_DNA.bam",
        RNAbam="filtered_bam/merged_RNA.bam"
    output:
        DNAmtx="DNA_matrix/matrix.mtx.gz",
        RNAmtx="RNA_matrix/matrix.mtx.gz"
    shell:
        '''
        bash /storage/zhangyanxiaoLab/xiongxiong/scripts/paired_tag/pipelines/bam2mtx.sh -n merged_DNA -g mm10 -i {input.DNAbam} -l DNA -o ./DNA_matrix
        bash /storage/zhangyanxiaoLab/xiongxiong/scripts/paired_tag/pipelines/bam2mtx.sh -n merged_RNA -g mm10 -i {input.RNAbam} -l DNA -o ./RNA_matrix
        '''


        


