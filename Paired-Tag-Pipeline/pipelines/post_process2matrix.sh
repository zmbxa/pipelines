ls dna/*/bam/*Useful*bam > tmp.dna_list.txt
ls rna/*/bam/*Useful*bam > tmp.rna_list.txt
paste tmp.dna_list.txt tmp.rna_list.txt -d, > bam_list.txt 

bash /storage/zhangyanxiaoLab/zhangyanxiao/github/Paired-Tag-Pipeline/pipelines/filter_merge_bam.sh -i bam_list.txt -d 0 -r 100
bash /storage/zhangyanxiaoLab/zhangyanxiao/github/Paired-Tag-Pipeline/pipelines/bam2mtx.sh -n merged_DNA -g mm10 -i filtered_bam/merged_RNA.bam -l RNA -o ./RNA_matrix_0_100



