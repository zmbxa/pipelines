#!/bin/bash

echo "barcode sample" > barcode.txt
echo "01 H3K9me3" >> barcode.txt
echo "02 H3K9me3" >> barcode.txt
echo "03 H3K9me3" >> barcode.txt
echo "04 H3K27ac" >> barcode.txt
echo "05 H3K27ac" >> barcode.txt
echo "06 H3K27ac" >> barcode.txt
echo "07 H3K9me3" >> barcode.txt
echo "08 H3K9me3" >> barcode.txt
echo "09 H3K9me3" >> barcode.txt
echo "10 H3K27ac" >> barcode.txt
echo "11 H3K27ac" >> barcode.txt
echo "12 H3K27ac" >> barcode.txt


echo "DNA_H3K27ac_merge.bam" > merge_list_dna_H3K27ac_bam.txt
echo "DNA_H3K9me3_merge.bam" > merge_list_dna_H3K9me3_bam.txt
echo "RNA_merge.bam" > merge_list_rna_bam.txt

dna_rna=("HSQ102_HSQ106" "HSQ103_HSQ107" "HSQ105_HSQ109")

for pair in "${dna_rna[@]}"
do
    dna=$(echo "$pair" | cut -d "_" -f 1)
    rna=$(echo "$pair" | cut -d "_" -f 2)
    
    # split bam by antibody
    dna_file="../../deep/$dna/$dna/bam/${dna}_UsefulReads.bam"
    samtools index $dna_file
    python ~/github/Paired-Tag-Pipeline/scripts/extract_sample_barcode_from_bam.py --bam $dna_file --statH barcode.txt  --outPrefix $dna
    
    rna_file="../../deep/$rna/$rna/bam/${rna}_UsefulReads.bam"
    samtools index $rna_file
    python ~/github/Paired-Tag-Pipeline/scripts/extract_sample_barcode_from_bam.py --bam $rna_file --statH barcode.txt  --outPrefix $rna
    
    # filter for H3K27ac dna: 500 rna: 500 
    bc_counts="/storage/zhangyanxiaoLab/qihongjian/projects/paired_tag/data/xx_pipeline/deep/qc/${dna}_${rna}_bc_counts.txt"
    awk 'NR == 1 || ($2 >= 500 && $3 >= 500)' $bc_counts > filtered_bc.txt
    
    samtools index $dna.metacell_H3K27ac.bam
    python ~/github/Paired-Tag-Pipeline_modified/scripts/filter_bam_based_on_bc.py --bam $dna.metacell_H3K27ac.bam --statH filtered_bc.txt
    echo "${dna}.metacell_H3K27ac_filtered_bc.bam ${dna}" >> merge_list_dna_H3K27ac_bam.txt

    samtools index $rna.metacell_H3K27ac.bam
    python ~/github/Paired-Tag-Pipeline_modified/scripts/filter_bam_based_on_bc.py --bam $rna.metacell_H3K27ac.bam --statH filtered_bc.txt
    
    
    # filter for H3K9me3 dna 200 rna 500
    awk 'NR == 1 || ($2 >= 200 && $3 >= 500)' $bc_counts > filtered_bc.txt
    samtools index $dna.metacell_H3K9me3.bam
    python ~/github/Paired-Tag-Pipeline_modified/scripts/filter_bam_based_on_bc.py --bam $dna.metacell_H3K9me3.bam --statH filtered_bc.txt
    echo "${dna}.metacell_H3K9me3_filtered_bc.bam ${dna}" >> merge_list_dna_H3K9me3_bam.txt
    
    samtools index $rna.metacell_H3K9me3.bam
    python ~/github/Paired-Tag-Pipeline_modified/scripts/filter_bam_based_on_bc.py --bam $rna.metacell_H3K9me3.bam --statH filtered_bc.txt
    
    # merge RNA
    samtools merge -f ${rna}_filtered.bam ${rna}.*filtered_bc.bam
    echo "${rna}_filtered.bam ${rna}" >> merge_list_rna_bam.txt
    
    echo $pair
done


echo "start to merge and add lib"
perl ~/github/Paired-Tag/perlscripts/mergeBam.pl merge_list_rna_bam.txt
perl ~/github/Paired-Tag/perlscripts/mergeBam.pl merge_list_dna_H3K27ac_bam.txt
perl ~/github/Paired-Tag/perlscripts/mergeBam.pl merge_list_dna_H3K9me3_bam.txt

wait

echo "start to get matrix" 

bash ~/github/Paired-Tag-Pipeline_modified/pipelines/bam2mtx.sh -n DNA_H3K9me3_merge -g mm10 -i DNA_H3K9me3_merge_sorted.bam -l DNA -o DNA_H3K9me3_matrix
bash ~/github/Paired-Tag-Pipeline_modified/pipelines/bam2mtx.sh -n DNA_H3K27ac_merge -g mm10 -i DNA_H3K27ac_merge_sorted.bam -l DNA -o DNA_H3K27ac_matrix
bash ~/github/Paired-Tag-Pipeline_modified/pipelines/bam2mtx.sh -n RNA_merge -g mm10 -i RNA_merge_sorted.bam -l RNA -o RNA_matrix

rm HSQ*
rm *merge.bam
rm *CountMatrix.txt
rm *ExpressedGenes.bed

# # filter high pileups. 
# python ~/github/Paired-Tag/remove_pileup/count_pileups.py PT_merge_sorted.bam PT_merge_sorted.tally.txt
# python ~/github/Paired-Tag/remove_pileup/remove_pileups.py PT_merge_sorted.bam PT_merge_sorted.tally.txt PT_merge_filtered.10.bam 10
# # bamcoverage 
# bamCoverage --bam PT_merge_sorted.bam --outFileFormat bigwig --outFileName hsq1216_prefilter.bw --binSize 25 --numberOfProcessors 8 --normalizeUsing RPKM
# bamCoverage --bam PT_merge_filtered.10.bam --outFileFormat bigwig --outFileName hsq1216_filtered.bw --binSize 25 --numberOfProcessors 8 --normalizeUsing RPKM


