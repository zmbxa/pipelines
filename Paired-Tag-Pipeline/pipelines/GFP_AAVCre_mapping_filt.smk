#snakefile 

import pandas as pd
import os

## setdb1 AAV-injected brain FC;Barcode: 01-06 control, 07-12 Cre
samples = pd.read_csv("sampleInfo.csv")

print(samples)

rule require_all:
    input:
        #"merge/bam_list.txt",
        expand("{sample}/{sample}/process/{sample}_CB.Counts",sample=samples["prefix"].unique()),
        expand("{sample}/{sample}/bam/{sample}_UsefulReads.bam",sample=samples["prefix"].unique()),
        "merge/filtered_bam/merged_DNA.bam",
        "merge/filtered_bam/merged_RNA.bam",
        "merge/DNA_matrix/matrix.mtx.gz",
        "merge/RNA_matrix/matrix.mtx.gz",
        "03_Sub-lib_Read_Number.pdf",
        "04_Barcode_Cumulative_Prop.pdf",
        "completeEmail.log"

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
        bash ~/pipelines/Paired-Tag-Pipeline/pipelines/pipeline_PT_gfpAAVCre.sh -g mm10 -l {params.library} -p 4
        """

rule list_make:
    output:
        "bam_list.txt"
    run:
        path=os.getcwd()
        for row in samples.head(round(len(samples)/2)).itertuples():
          f=open("bam_list.txt","a")
          f.write(path+"/"+getattr(row,'prefix')+"/"+getattr(row,'prefix')+"/bam/"+getattr(row,'prefix')+"_UsefulReads.bam"+","+path+"/"+getattr(row,'exp_prefix')+"/"+getattr(row,'exp_prefix')+"/bam/"+getattr(row,'exp_prefix')+"_UsefulReads.bam")
          f.write('\n')
        f.close()

rule merge:
    input:
        expand("{sample}/{sample}/bam/{sample}_UsefulReads.bam",sample=samples["prefix"].unique()),
        bamlist="bam_list.txt"
    output:
        "merge/filtered_bam/merged_DNA.bam",
        "merge/filtered_bam/merged_RNA.bam"
    conda:
        "paired_tag"
    shell:
        '''
        mkdir -p merge
        mv {input.bamlist} merge/bam_list.txt
        cd merge
        bash /storage/zhangyanxiaoLab/xiongxiong/scripts/paired_tag/pipelines/filter_merge_bam.sh -i bam_list.txt  -d 300 -r 500
        cd ..
        '''
rule bam2mtx:
    input:
        "merge/filtered_bam/merged_DNA.bam",
        "merge/filtered_bam/merged_RNA.bam"
    output:
        DNA="merge/DNA_matrix/matrix.mtx.gz",
        RNA="merge/RNA_matrix/matrix.mtx.gz"
    conda:
        "paired_tag"
    shell:
        '''
        cd merge
        bash /storage/zhangyanxiaoLab/niuyuxiao/pipelines/Paired-Tag-Pipeline/pipelines/bam2mtx_gfpAAVCre.sh -n merged_DNA -g mm10 -i filtered_bam/merged_DNA.bam -l DNA -o ./DNA_matrix
        bash /storage/zhangyanxiaoLab/niuyuxiao/pipelines/Paired-Tag-Pipeline/pipelines/bam2mtx_gfpAAVCre.sh -n merged_RNA -g mm10 -i filtered_bam/merged_RNA.bam -l RNA -o ./RNA_matrix
        cd ..
        '''

rule QC_plot:
    input:
        expand("{sample}/{sample}/process/{sample}_CB.Counts",sample=samples["prefix"].unique())
    output:
        "03_Sub-lib_Read_Number.pdf",
        "04_Barcode_Cumulative_Prop.pdf"
    run:
        commond1="Rscript /storage/zhangyanxiaoLab/xiongxiong/scripts/paired_tag/analysis/dot_reads_numbers.R 300 500 "
        path=os.getcwd()
        for row in samples.head(round(len(samples)/2)).itertuples():
          commond1=commond1+path+"/"+getattr(row,'prefix')+"/"+getattr(row,'prefix')+"/process/"+getattr(row,'prefix')+"_CB.Counts,"+path+"/"+getattr(row,'exp_prefix')+"/"+getattr(row,'exp_prefix')+"/process/"+getattr(row,'exp_prefix')+"_CB.Counts "
        commond2="Rscript /storage/zhangyanxiaoLab/xiongxiong/scripts/paired_tag/analysis/line_barcodeProp.R "
        path=os.getcwd()
        for row in samples.itertuples():
          commond2=commond2+path+"/"+getattr(row,'prefix')+"/"+getattr(row,'prefix')+"/process/"+getattr(row,'prefix')+"_CB.Counts "
        os.system(commond1)
        os.system(commond2)

rule reminder:
    input:
        "merge/DNA_matrix/matrix.mtx.gz",
        "merge/RNA_matrix/matrix.mtx.gz",
    output:
        "completeEmail.log"
    shell:
        '''
        python ~/pipelines/reminder.py "Your paired-tag pipeline are successfully done! Congrats~" "Need a cup of Latte?" > {output}
        '''
        
        


