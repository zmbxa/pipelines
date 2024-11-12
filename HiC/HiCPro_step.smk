### This snakemake file can process HiC data data ###
### The workdir need to be changed as data containing folder (include fastq/) and all .fq file in fastq/[name]/[name]_R1/R2.fastq.gz ###
### !!! Change SAMPLES to a Tuple contain sample names (for snakemake final output)!!! ###
### workflow is similar to HiC-Pro but will work in parallel ###
### Usage:  snakemake --cores 4 -s HiCPro.smk --use-conda --conda-frontend conda

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
ALL_SAMPLES = FASTQ_DICT.keys()
GENOME = config["GENOME"]
ENZYME = config["ENZYME"]
if ENZYME == "DpnII":
    CUTSITE = "GATC"
elif ENZYME == "hindIII":
    CUTSITE = "AAGCTAGCTT"
elif ENZYME == "HaeIII":
    CUTSITE = "GGCCGGCC"
else:
    raise ValueError(f"Unsupported enzyme: {ENZYME}")

FRAGFILE = "/storage/zhangyanxiaoLab/niuyuxiao/tools/HiC-Pro-3.1.0/annotation/"+GENOME+"_"+ENZYME+".bed"
CHROM_SIZE="/storage/zhangyanxiaoLab/niuyuxiao/tools/HiC-Pro-3.1.0/annotation/chrom_"+GENOME+".sizes"
BOWTIE_DICT = {"hg38":"/storage/zhangyanxiaoLab/share/bowtie2_index/hg38",
            "mm10":"/storage/zhangyanxiaoLab/share/bowtie2_index/mm10",
            "mm10_with_lambda":"/storage/zhangyanxiaoLab/niuyuxiao/annotations/bowtie2_index/mm10_lambda/mm10_lambda",
            "mm10_GFP":"/storage/zhangyanxiaoLab/share/bowtie2_index/GFP_AAVCre_index/mm10_gfp_AAVcre",
            "dm6":"/storage/zhangyanxiaoLab/share/bowtie2_index/dm6",
            "ce11":"/storage/zhangyanxiaoLab/share/bowtie2_index/ce11",
            "lambda":"/storage/zhangyanxiaoLab/niuyuxiao/annotations/bowtie2_index/lambda/lambda"
            }
bt_index = BOWTIE_DICT[GENOME]



### Samples ###
# ALL_SAMPLES = ["WT", "KO"]

### Define the size of the genomic bins used to generate the HiC contact matrix ###
ALL_BIN_SIZE = [2000,5000,10000,20000,40000,150000,500000,1000000]

### Define the rules ###
ALL_RAW_MAPS = expand("HiCPro_out/hic_results/matrix/{sample}/raw/{bsize}/{sample}_{bsize}{suffix}", sample = ALL_SAMPLES, bsize = ALL_BIN_SIZE, suffix = ["_abs.bed", ".matrix"] * len(ALL_BIN_SIZE))
ALL_ICE_MATRIX = expand("HiCPro_out/hic_results/matrix/{sample}/iced/{bsize}/{sample}_{bsize}_iced.matrix", sample = ALL_SAMPLES, bsize = ALL_BIN_SIZE)


### Rule all ###
rule all:
	input:
		ALL_ICE_MATRIX


### List of rules to be divided into rules .smk format file ###

rule fastp:
	input:
		fastq1="fastq_hicpro/{sample}/{sample}_R1.fastq.gz",
		fastq2="fastq_hicpro/{sample}/{sample}_R2.fastq.gz"
	output:
		trimmed1="fastp/{sample}_R1.fastp.fastq.gz",
		trimmed2="fastp/{sample}_R2.fastp.fastq.gz"
	params:
		parameters="-w 5 -t 1 -A -Q -L"
	log:
		fastpLog="logs/{sample}.fastpLog"
	shell:
		"fastp -i {input.fastq1} -I {input.fastq2} \
		-o {output.trimmed1} -O {output.trimmed2} \
		{params.parameters} \
		2> {log.fastpLog}"


rule alignStep1:
	input:
		aln1="fastp/{sample}_{pos}.fastp.fastq.gz"
	output:
		mapped="mapped_reads/{sample}_{pos}.bam",
		unmapped="unmapped_reads/{sample}_{pos}.unmap.fastq"
	params:
		index=bt_index,
		log="los/{sample}_{pos}.bowtie2.log",
		sm="SM:{sample}_{pos}"
	log:
		alignStep1_log="logs/{sample}_{pos}.alignStep1.log"

	shell:
		"bowtie2 --very-sensitive \
		-L 30 \
		--score-min \
		L,-0.6,-0.2 \
		--end-to-end \
		--reorder \
		--un {output.unmapped} \
		--rg-id BMG --rg \
		{params.sm} -p 3 \
		-x {params.index} \
		-U {input.aln1} \
		2>> {log.alignStep1_log} | \
		samtools view -F 4 -bS - > {output.mapped}"

rule cutsite_trimming:
	input:
		tocut1="unmapped_reads/{sample}_{pos}.unmap.fastq"
	output:
		cut1="unmapped_reads/{sample}_{pos}.unmap.trimmed.fastq"
	params:
		cutSite = CUTSITE
	log:
		cutLog="logs/{sample}_{pos}.unmapped_trimming.log"
	shell:
		"/storage/zhangyanxiaoLab/niuyuxiao/tools/HiC-Pro-3.1.0/scripts/cutsite_trimming \
		--fastq {input.tocut1} \
		--cutsite {params.cutSite} \
		--out {output.cut1} > \
		{log.cutLog} \
		2>&1"

rule alignStep2:
	input:
		aln1="unmapped_reads/{sample}_{pos}.unmap.trimmed.fastq"
	output:
		mapped="mapped_reads/{sample}_{pos}.bwt2glob.unmap_bwt2loc.bam"
	params:
		index=bt_index,
		sm="SM:{sample}_{pos}_genome.bwt2glob.unmap"
	log:
		alignStep2_log="logs/{sample}_{pos}.alignStep2.log"
	shell:
		"bowtie2 --very-sensitive \
		-L 20 \
		--score-min \
		L,-0.6,-0.2 \
		--end-to-end \
		--reorder \
		--rg-id BML --rg \
		{params.sm} -p 8 \
		-x {params.index} \
		-U {input.aln1} \
		2>> {log.alignStep2_log} | \
		samtools view -bS - > {output.mapped}"

rule mappingCombine_merge:
	input:
		bam1="mapped_reads/{sample}_{pos}.bam",
		bam2="mapped_reads/{sample}_{pos}.bwt2glob.unmap_bwt2loc.bam"
	output:
		mergedBam="mapped_reads/bwt2/{sample}/{sample}_{pos}.bwt2merged.bam"
	params:
		thread=10,
		parameters="-n -f"
	shell:
		"samtools merge \
		-@ {params.thread} \
		{params.parameters} \
		{output.mergedBam} \
		{input.bam1} \
		{input.bam2}"

rule mappingCombine_sort:
	input:
		toSort="mapped_reads/bwt2/{sample}/{sample}_{pos}.bwt2merged.bam"
	output:
		sorted="mapped_reads/bwt2/{sample}/{sample}_{pos}.bwt2merged.sorted.bam"
	params:
		thread=10,
		memory=12,
		parameters="-n -T",
		temporary="tmp/{sample}_{pos}.genome"
	shell:
		"samtools sort \
		-@ {params.thread} \
		-m {params.memory}G \
		{params.parameters} {params.temporary} \
		-o {output.sorted} \
		{input.toSort}"

rule mergeSAM:
	input:
		toMerge_f="mapped_reads/bwt2/{sample}/{sample}_R1.bwt2merged.sorted.bam",
		toMerge_r="mapped_reads/bwt2/{sample}/{sample}_R2.bwt2merged.sorted.bam"
	output:
		mergedBam="mapped_reads/bwt2/{sample}/{sample}.bwt2pairs.bam"
	params:
		parameters="-q 10 -t -v"
	log:
		mergesam="{sample}.mergeSam.log"
	shell:
		"python /storage/zhangyanxiaoLab/niuyuxiao/tools/HiC-Pro-3.1.0/scripts/mergeSAM.py \
		{params.parameters} \
		-f {input.toMerge_f} \
		-r {input.toMerge_r} \
		-o {output.mergedBam} \
		2>{log.mergesam}"

rule mapped_2_hic_fragments:
	input:
		toMap="mapped_reads/bwt2/{sample}/{sample}.bwt2pairs.bam"
	output:
		mapped="HiCPro_out/hic_results/data/{sample}/{sample}.bwt2pairs.{suffix}"
	params:
		digestedGenome=FRAGFILE,
		parameters="-v -a",
		outdir="HiCPro_out/hic_results/data/{sample}"
	shell:
		"""
		python /storage/zhangyanxiaoLab/niuyuxiao/tools/HiC-Pro-3.1.0/scripts/mapped_2hic_fragments.py \
		{params.parameters} \
		-f {params.digestedGenome} \
		-r {input.toMap} \
		-o {params.outdir}
		
		"""

rule merge_valid_interactions:
	input:
		toMerge="HiCPro_out/hic_results/data/{sample}/{sample}.bwt2pairs.validPairs"
	output:
		Merged="HiCPro_out/hic_results/data/{sample}/{sample}.allValidPairs"
	shell:
		"""
		LANG=en; sort -T tmp -S 50% -k2,2V -k3,3n -k5,5V -k6,6n \
		-m {input.toMerge} | \
		awk -F"\t" 'BEGIN{{c1=0;c2=0;s1=0;s2=0}}(c1!=$2 || c2!=$5 || s1!=$3 || s2!=$6){{print;c1=$2;c2=$5;s1=$3;s2=$6}}' > \
		{output.Merged}

		"""

rule build_raw_maps:
	input:
		toBuild="HiCPro_out/hic_results/data/{sample}/{sample}.allValidPairs"
	output:
		Built="HiCPro_out/hic_results/matrix/{sample}/raw/{bsize}/{sample}_{bsize}{suffix}"
	params:
		chrSize=CHROM_SIZE,
		ifile="/dev/stdin",
		m_format="upper",
		bsize_par=ALL_BIN_SIZE,
		Built_par="HiCPro_out/hic_results/matrix/{sample}/raw/{bsize}/{sample}_{bsize}"
	shell:
		"""
		cat {input.toBuild} | \
		/storage/zhangyanxiaoLab/niuyuxiao/tools/HiC-Pro-3.1.0/scripts/build_matrix \
		--matrix-format {params.m_format} \
		--binsize {params.bsize_par} \
		--chrsizes {params.chrSize} \
		--ifile {params.ifile} \
		--oprefix {params.Built_par}
		"""

rule ice_normalisation:
	input:
		toNorm="HiCPro_out/hic_results/matrix/{sample}/raw/{bsize}/{sample}_{bsize}.matrix"
	output:
		Norm="HiCPro_out/hic_results/matrix/{sample}/iced/{bsize}/{sample}_{bsize}_iced.matrix"
	params:
		filter_Low = 0.02,
		filter_High = 0,
		max_Iter = 100,
		eps = 0.1,
		bias = 1,
		verbose = 1
	log:
		ice_logs = "logs/{sample}_{bsize}_ice.log"
	shell:
		"""
		/storage/zhangyanxiaoLab/niuyuxiao/anaconda3/envs/HiCPro/bin/ice --results_filename {output.Norm} \
		--filter_low_counts_perc {params.filter_Low} \
		--filter_high_counts_perc {params.filter_High} \
		--max_iter {params.max_Iter} \
		--eps {params.eps} \
		--remove-all-zeros-loci \
		--output-bias {params.bias} \
		--verbose {params.verbose} \
		{input.toNorm} \
		>{log.ice_logs}

		"""
rule avp2hic:
  input:
    allVP="HiCPro_out/hic_results/data/{sample}/{sample}.allValidPairs",
  params:
    chromSize=CHROM_SIZE,
  output:
    hicfile="hicFile_juicer/{sample}_srt.hic",
    juicebox=temp("hicFile_juicer/{sample}_allValidPairs.pre_juicebox_sorted"),
  shell:
    '''
    echo "Generating Juicebox input files for ${wildcards.sample}..."
    awk '{{$4=$4!="+"; $7=$7!="+"}} $2<=$5{{print $1, $4, $2, $3, 0, $7, $5, $6, 1, $11, $12 }}$5<$2{{ print $1, $7, $5, $6, 0, $4, $2, $3, 1, $12, $11 }}' {input.allVP} | sort -T ./ -k3,3d  -k7,7d -S 250G  > {output.juicebox}
    
    echo "Running Juicebox for {wildcards.sample}..."
    java -XX:ParallelGCThreads=20 -Xmx50g -jar ~/tools/juicer/juicer/CPU/common/juicer_tools_1.22.01.jar pre {output.juicebox} {output.hicfile} {params.chromSize}
    '''

rule QC:
  input:
    mpairStat=expand("HiCPro_out/hic_results/stats/{sample}/{sample}_allValidPairs.mergestat",sample=ALL_SAMPLES),
    avp_stat=expand("HiCPro_out/hic_results/stats/{sample}/{sample}_allValidPairs.mergestat",sample=ALL_SAMPLES),
  output:
    qc="HiCPro_out/all_sample_qc.txt"
  shell:
    '''
    cd HiCPro_out/
    python /storage/zhangyanxiaoLab/niuyuxiao/pipelines/HiC/QC_HiCPro.py  -s {ALL_SAMPLES}
    cd ..
    python ~/pipelines/reminder.py "Youcan check your hic qc now!" "$(cat {output.qc})"
    '''

	