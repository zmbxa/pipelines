#!/bin/bash
##-------------------------
# @Author: XiongXiong
# @Date: 2023/1/20
##-------------------------


bash /storage/zhangyanxiaoLab/xiongxiong/scripts/paired_tag/pipelines/pipeline_PT.sh -g mm10 -l DNA
bash /storage/zhangyanxiaoLab/xiongxiong/scripts/paired_tag/pipelines/pipeline_PT.sh -g mm10 -l RNA

Rscript /storage/zhangyanxiaoLab/xiongxiong/scripts/paired_tag/analysis/dot_reads_numbers.R 300 500 \
	/storage/zhangyanxiaoLab/xiongxiong/Data/Mouse/liver/DNA/GY016_process/GY016_CB.Counts,/storage/zhangyanxiaoLab/xiongxiong/Data/Mouse/liver/RNA/GY020_process/GY020_CB.Counts \
	/storage/zhangyanxiaoLab/xiongxiong/Data/Mouse/liver/DNA/GY018_process/GY018_CB.Counts,/storage/zhangyanxiaoLab/xiongxiong/Data/Mouse/liver/RNA/GY022_process/GY022_CB.Counts

Rscript /storage/zhangyanxiaoLab/xiongxiong/scripts/paired_tag/analysis/line_barcodeProp.R \
	/storage/zhangyanxiaoLab/xiongxiong/Data/Mouse/liver/DNA/GY016_process/GY016_CB.Counts \
	/storage/zhangyanxiaoLab/xiongxiong/Data/Mouse/liver/DNA/GY018_process/GY018_CB.Counts \
	/storage/zhangyanxiaoLab/xiongxiong/Data/Mouse/liver/RNA/GY020_process/GY020_CB.Counts \
	/storage/zhangyanxiaoLab/xiongxiong/Data/Mouse/liver/RNA/GY022_process/GY022_CB.Counts



bash /storage/zhangyanxiaoLab/xiongxiong/scripts/paired_tag/pipelines/filter_merge_bam.sh -i bam_list.txt -d 300 -r 500
bam_list.txt:
	/storage/zhangyanxiaoLab/xiongxiong/Data/projects/mouse_brain/paired-tag/batch1/DNA/HSQ043_bam/HSQ043_UsefulReads.bam,/storage/zhangyanxiaoLab/xiongxiong/Data/projects/mouse_brain/paired-tag/batch1/RNA/HSQ047_bam/HSQ047_UsefulReads.bam \
	/storage/zhangyanxiaoLab/xiongxiong/Data/projects/mouse_brain/paired-tag/batch1/DNA/HSQ044_bam/HSQ044_UsefulReads.bam,/storage/zhangyanxiaoLab/xiongxiong/Data/projects/mouse_brain/paired-tag/batch1/RNA/HSQ048_bam/HSQ048_UsefulReads.bam \


bash /storage/zhangyanxiaoLab/xiongxiong/scripts/paired_tag/pipelines/bam2mtx.sh -n merged_DNA -g mm10 -i /storage/zhangyanxiaoLab/xiongxiong/Data/Mouse/liver/merge/filtered_bam/merged_DNA.bam -l DNA -o ./DNA_martrix
bash /storage/zhangyanxiaoLab/xiongxiong/scripts/paired_tag/pipelines/bam2mtx.sh -n merged_RNA -g mm10 -i /storage/zhangyanxiaoLab/xiongxiong/Data/Mouse/liver/merge/filtered_bam/merged_RNA.bam -l RNA -o ./RNA_martrix

Rscript /storage/zhangyanxiaoLab/xiongxiong/scripts/paired_tag/analysis/seurat.R \
	/storage/zhangyanxiaoLab/xiongxiong/Data/Mouse/liver/merge/RNA_martrix/

bash /storage/zhangyanxiaoLab/xiongxiong/scripts/paired_tag/pipelines/pseudo_bulk.sh -i /storage/zhangyanxiaoLab/xiongxiong/Data/Mouse/liver/merge/filtered_bam/merged_DNA.bam -t Cell_Clusters.txt
batch1_file_list.txt: <brcode><TAB><cell_cluster>.
