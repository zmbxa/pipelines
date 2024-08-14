#cd $1
echo $pwd
bams=$(ls */bam/*_UsefulReads.bam)

for bam in $bams ; 
do 
samtools index $bam
bamCoverage --bam $bam --outFileFormat bigwig --outFileName $bam.bw --binSize 25 --numberOfProcessors 8 --normalizeUsing RPKM
done

#tss_bed=/storage/zhangyanxiaoLab/zhangyanxiao/annotations/mm10/mm10.gencode.vM25.annotation.gene.tss1k.bed
tss_bed=/storage/zhangyanxiaoLab/xiongxiong/Reference/refBed/mm10_TSS.bed
#bw=$(ls */bw/*.bw)
bw=$(ls */bam/*.bw)
computeMatrix reference-point -p 20 -S ${bw} \
    --regionsFileName $tss_bed \
    --skipZeros \
    -a 2000 \
    -b 2000 \
    -bs 20 \
    -o output.tss_matrix.gz \
    --maxThreshold 5000

plotProfile -m output.tss_matrix.gz \
    -out tss_plot.pdf \
    --perGroup \
    --refPointLabel TSS \
    --dpi 500 \
    --averageType mean \
    --plotTitle "Test" 
