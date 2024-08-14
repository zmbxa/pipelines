#!/bin/bash
##-------------------------
# @Author: XiongXiong
# @Date: 2023/7/8
##-------------------------


prefix=$1
reference_genome=$2
library=$3
nthread=$4
debug=$5

if [ $reference_genome == "mm10_gfp" ];then 
# echo mm10_gfp
bt_ref=/storage/zhangyanxiaoLab/niuyuxiao/annotations/bowtie2_index/GFP_AAVCRE_ERT_index/mm10_gfp_AAVCRE_ERT
star_ref=/storage/zhangyanxiaoLab/niuyuxiao/annotations/STAR_index/mm10_gfp_AAVCre_ERT
else
bt_ref=/storage/zhangyanxiaoLab/share/bowtie2_index/${reference_genome}
star_ref=/storage/zhangyanxiaoLab/share/STAR_index/${reference_genome}
fi




    # test if a report is completed.
    if [ -f ${prefix}/${prefix}.Report.txt ]
    then
    if (( $(grep Completed  ${prefix}/${prefix}.Report.txt|wc -l) == 1))
    then echo "already done!"
	 exit
    fi
    fi


    echo $prefix $reference_genome $library $nthread


		mkdir ./${prefix}
		cd ./${prefix}

		zcat ../fastq/${prefix}_R2.fastq.gz | awk '{if(NR%4==1) print $1; else print $0}' > tmp.R2.fq
		awk '{if(NR%4==1) print $1}' tmp.R2.fq > tmp.R2.ReadNames
		awk 'NR%4==2' tmp.R2.fq > tmp.R2.seqs
		paste tmp.R2.ReadNames tmp.R2.seqs > tmp.R2.Reads

		zcat ../fastq/${prefix}_R1.fastq.gz | awk '{if(NR%4==1) print $1; else print $0}' > tmp.R1.fq

		echo -e "# reference genome: ${reference_genome}\n# library: ${library}\n"  | tee ${prefix}.Report.txt
		echo -e "============= ${prefix} =============" | tee -a ${prefix}.Report.txt
		nRaw=`cat tmp.R2.ReadNames | wc -l`
		pRaw=`echo "scale=1; 100*${nRaw}/${nRaw}" | bc`
		echo "Raw: ${nRaw} (${pRaw}%)" | tee -a ${prefix}.Report.txt

		case ${library} in
			"DNA")
				awk '{if($2~/T[A,T,C,G]{4}CCTGCA/) print $1, "Contamination"}' OFS='\t' tmp.R2.Reads > tmp.ReadType.txt
				cut -f 1 tmp.ReadType.txt > tmp.conta.names
				awk '{print substr($1,18,150)}' tmp.R2.seqs | awk -F "GTGGCC|ATCCAC|GCGGCC" '{print "N"$1,$2,$3,$4,$5,$6}' > tmp.seqs
				awk '{print "GTGGCC"substr($2,1,length($2)-7)"ATCCAC"substr($3,1,length($3)-4)}' tmp.seqs | paste tmp.R2.ReadNames - | grep -v -w -F -f tmp.conta.names - | awk '{if(length($2)<66 && length($2)>54 && substr($2,length($2),1)=="A") print ">"substr($1,2,length($1))"\n"substr($2,1,length($2)-1)}' > ${prefix}_Bait.fa;;
			"RNA")
				awk '{if($2~/A[A,T,C,G]{4}GCGGCC/) print $1, "Contamination"}' OFS='\t' tmp.R2.Reads > tmp.ReadType.txt
				cut -f 1 tmp.ReadType.txt > tmp.conta.names
				awk '{print substr($1,18,150)}' tmp.R2.seqs | awk -F "GTGGCC|ATCCAC|CCTGCA" '{print "N"$1,$2,$3,$4,$5,$6}' > tmp.seqs
				awk '{print "GTGGCC"substr($2,1,length($2)-7)"ATCCAC"substr($3,1,length($3)-4)}' tmp.seqs | paste tmp.R2.ReadNames - | grep -v -w -F -f tmp.conta.names - | awk '{if(length($2)<66 && length($2)>54 && substr($2,length($2),1)=="T") print ">"substr($1,2,length($1))"\n"substr($2,1,length($2)-1)}' > ${prefix}_Bait.fa;;
			*)
			 echo "Unknown data type"; break;;
		esac

		bowtie2 -x /storage/zhangyanxiaoLab/xiongxiong/index/bowtie2/paired_tag_barcode_bait/Bait -p ${nthread} --local --no-unal --norc -f ${prefix}_Bait.fa | samtools view -@ ${nthread} -bhS - > ${prefix}_Bait.bam
		samtools view -@ ${nthread} ${prefix}_Bait.bam | awk '{print "@"$1}' > tmp.FL.txt
		paste tmp.R2.ReadNames tmp.R2.seqs tmp.seqs | grep -w -F -f tmp.FL.txt - | awk '{print ">"substr($1,2,length($1))"\n"substr($2,11,7)""substr($4,length($4)-6,7)""substr($5,length($5)-3,4)}' OFS='\t' > ${prefix}_Barcode.fa
		bowtie /storage/zhangyanxiaoLab/xiongxiong/index/bowtie/paired_tag_cellular_barcode/PT_Barcode -p ${nthread} --no-unal -v 1 -m 1 --norc -f ${prefix}_Barcode.fa -S tmp.BC.sam

		nAssignCB=`samtools view -@ ${nthread} tmp.BC.sam | wc -l `
		nFaCB=`cat ${prefix}_Barcode.fa | wc -l `
		pAssignCB=`echo "scale=0; 200*${nAssignCB}/${nFaCB}" | bc`
		
		if [ ${pAssignCB} -lt 90 ]; then
			echo "Trying new version of barcode 2, as less than 90% barcode successfully assigned to cellular barcode..."
			bowtie /storage/zhangyanxiaoLab/xiongxiong/index/bowtie/paired_tag_cellular_barcode/new_BC2/PT_Barcode -p ${nthread} --no-unal -v 1 -m 1 --norc -f ${prefix}_Barcode.fa -S tmp.BC.sam
			sed -i '3i# barcode: new' ${prefix}.Report.txt
		else
			sed -i '3i# barcode: old' ${prefix}.Report.txt
		fi
		
		samtools view -@ ${nthread} -bhS tmp.BC.sam > ${prefix}_Barcode.bam		
		
		awk '{print $1, "FullyLigated"}' OFS='\t' tmp.FL.txt >> tmp.ReadType.txt
		cut -f 1 tmp.ReadType.txt | grep -v -w -F -f - tmp.R2.Reads | awk '{if($2~/^G{15,}/) print $1, "polyG"}' OFS='\t' >> tmp.ReadType.txt
		cut -f 1 tmp.ReadType.txt | grep -v -w -F -f - tmp.R2.ReadNames | awk '{print $1, "Random"}' OFS='\t' >> tmp.ReadType.txt

		awk '{if($2=="FullyLigated") print $1}' tmp.ReadType.txt | grep -A 3 -w -F -f - tmp.R1.fq | awk '$0 != "--"' | awk '{if(NR%4==3) print "+"; else print $0}' > tmp.FullyLigated.R1.fq
		awk '{if($2=="FullyLigated") print $1}' tmp.ReadType.txt | grep -A 3 -w -F -f - tmp.R2.fq | awk '$0 != "--"' | awk '{if(NR%4==3) print "+"; else print $0}' > tmp.FullyLigated.R2.fq

		nFL=`awk '$2=="FullyLigated"' tmp.ReadType.txt | wc -l `
		pFL=`echo "scale=1; 100*${nFL}/${nRaw}" | bc`
		echo "Fully ligated: ${nFL} (${pFL}%)" | tee -a ${prefix}.Report.txt

		# Reads without genomic sequence insert
		case ${library} in
			"DNA")
				nUnload=`awk 'substr($1,1,35)~/CTGTCTCTTATACAC/' tmp.FullyLigated.R1.fq | wc -l `
				nUnload2=`awk 'substr($1,1,35)~/CTCAAGCACGTGGAT/' tmp.FullyLigated.R1.fq | wc -l `;;
			"RNA")
				nUnload=`awk '$0~/G{5,}[A,T,C,G]{,30}A{10,}/' tmp.FullyLigated.R1.fq | wc -l `
				nUnload2=`awk '$0~/AAGCAGTG/ || $0~/GTATCAA/' tmp.FullyLigated.R1.fq | wc -l `;;
			*)
			 echo "Unknown data type"; break;;
		esac

		pUnload=`echo "scale=1; 100*${nUnload}/${nFL}" | bc`
		echo -e "+++ Unload: ${nUnload} (${pUnload}%)" | tee -a ${prefix}.Report.txt

		pUnload2=`echo "scale=1; 100*${nUnload2}/${nFL}" | bc`
		echo -e "+++ Unload(no bc1): ${nUnload2} (${pUnload2}%)" | tee -a ${prefix}.Report.txt

		nConta=`awk '$2=="Contamination"' tmp.ReadType.txt | wc -l `
		pConta=`echo "scale=1; 100*${nConta}/${nRaw}" | bc`
		echo "Cross Contamination: ${nConta} (${pConta}%)" | tee -a ${prefix}.Report.txt

		nRandom=`awk '$2=="Random"' tmp.ReadType.txt | wc -l `
		pRandom=`echo "scale=1; 100*${nRandom}/${nRaw}" | bc`
		echo "Random: ${nRandom} (${pRandom}%)" | tee -a ${prefix}.Report.txt

		nPolyG=`awk '$2=="polyG"' tmp.ReadType.txt | wc -l `
		pPolyG=`echo "scale=1; 100*${nPolyG}/${nRaw}" | bc`
		echo "polyG: ${nPolyG} (${pPolyG}%)" | tee -a ${prefix}.Report.txt

		nTn5ME=`awk '{print $1, substr($2,11,140)}' OFS='\t' tmp.R2.Reads | grep "AGATGTGTATAAGAGACAG" | wc -l `
		pTn5ME=`echo "scale=1; 100*${nTn5ME}/${nRaw}" | bc`
		echo -e "\nTn5ME: ${nTn5ME} (${pTn5ME}%)" | tee -a ${prefix}.Report.txt

		case ${library} in
			"DNA")
				nBC1=`awk '{print $1, substr($2,11,140)}' OFS='\t' tmp.R2.Reads | grep "GGCCAGAGCATTCGA" | wc -l `;;
			"RNA")
				nBC1=`awk '{print $1, substr($2,11,140)}' OFS='\t' tmp.R2.Reads | grep "GGCCAGAGCATTCGT" | wc -l `;;
			*)
			 echo "Unknown data type"; break;;
		esac

		pBC1=`echo "scale=1; 100*${nBC1}/${nRaw}" | bc`
		echo "Barcode1: ${nBC1} (${pBC1}%)" | tee -a ${prefix}.Report.txt

		nBC2=`awk '{print $1, substr($2,11,140)}' OFS='\t' tmp.R2.Reads | grep "GTGCGAACTCAGACC" | wc -l `
		pBC2=`echo "scale=1; 100*${nBC2}/${nRaw}" | bc`
		echo "Barcode2: ${nBC2} (${pBC2}%)" | tee -a ${prefix}.Report.txt

		nBC3=`awk '{print $1, substr($2,11,140)}' OFS='\t' tmp.R2.Reads | grep "GTGGCCGATGT" | wc -l `
		pBC3=`echo "scale=1; 100*${nBC3}/${nRaw}" | bc`
		echo "Barcode3: ${nBC3} (${pBC3}%)" | tee -a ${prefix}.Report.txt

		#----------------------------------------
		# extract UMI, adapter trimming
		fastp --umi --umi_loc read2 --umi_len 10 --json tmp.FullyLigated.json --html tmp.FullyLigated.QC.html --correction -q 20 -e 20 -l 30 --thread 16 --in1 tmp.FullyLigated.R1.fq --out1 tmp.FullyLigated.QC.R1.fq.gz --in2 tmp.FullyLigated.R2.fq --out2 tmp.FullyLigated.QC.R2.fq.gz 2> /dev/null
		trim_galore --no_report_file --cores 6 --gzip tmp.FullyLigated.QC.R1.fq.gz 2> /dev/null

		case ${library} in
			"DNA")
				trim_galore --no_report_file --cores 6 -a AAACATCGGCCAC --gzip tmp.FullyLigated.QC.R1_trimmed.fq.gz 2> /dev/null
				trim_galore --no_report_file --cores 6 -a CAAGCACGTGGAT --gzip tmp.FullyLigated.QC.R1_trimmed_trimmed.fq.gz 2> /dev/null
				trim_galore --no_report_file --cores 6 -a CATCTGCGGCCGC --gzip tmp.FullyLigated.QC.R1_trimmed_trimmed_trimmed.fq.gz 2> /dev/null
				trim_galore --no_report_file --cores 6 -a CTGTCTCTTATAC --gzip tmp.FullyLigated.QC.R1_trimmed_trimmed_trimmed_trimmed.fq.gz 2> /dev/null ## trim barcode2 directly ligated with barcode3 + barcode2 + barcode1
				mv tmp.FullyLigated.QC.R1_trimmed_trimmed_trimmed_trimmed_trimmed.fq.gz ${prefix}_FullyLigated_trimmed_R1.fastq.gz

				QC=`zcat ${prefix}_FullyLigated_trimmed_R1.fastq.gz | wc -l `
				nQC=`echo "${QC}/4" | bc`
				pQC=`echo "scale=1; 100*${nQC}/${nFL}" | bc`
				echo -e "\nQuality control: ${nQC} (${pQC}%)" | tee -a ${prefix}.Report.txt

				bowtie2 -p ${nthread} -x ${bt_ref} -U ${prefix}_FullyLigated_trimmed_R1.fastq.gz | samtools view -bh -@ ${nthread} > tmp.${prefix}_FullyLigated.bam
				samtools view -@ ${nthread} -bh -F 260 tmp.${prefix}_FullyLigated.bam > ${prefix}_FullyLigated_mapped.bam
				samtools view -@ ${nthread} -bh -f 4 tmp.${prefix}_FullyLigated.bam > ${prefix}_FullyLigated_unmapped.bam
				;;
			"RNA")
				awk '{if(NR%4==1) print $1}' tmp.R1.fq > tmp.R1.ReadNames
				awk 'NR%4==2' tmp.R1.fq > tmp.R1.seqs
				paste tmp.R1.ReadNames tmp.R1.seqs > tmp.R1.Reads
				awk '{if($2!~/G{5,}[A,T,C,G]{,30}A{10,}/ && $2!~/AAGCAGTG/ && $2!~/GTATCAA/) print $1}' tmp.R1.Reads > tmp.loaded.names
				zcat tmp.FullyLigated.QC.R1_trimmed.fq.gz | grep -A 3 -w -F -f tmp.loaded.names - | awk '$0 != "--"' > tmp.FullyLigated.QC.R1_trimmed.fq
				cutadapt --cores 6 -g GGGGGGGGG tmp.FullyLigated.QC.R1_trimmed.fq > tmp.FullyLigated.QC.R1_trimmed_trimmed.fq 2> /dev/null
				trim_galore --no_report_file --cores 6 --gzip -a ACGAATGCTCTGGC tmp.FullyLigated.QC.R1_trimmed_trimmed.fq 2> /dev/null ## trim barcode2 directly ligated with barcode3 + barcode2 + barcode1
				trim_galore --no_report_file --cores 6 -a CCTGCAGGNNNNACGAATGCTCTGGCCTNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN tmp.FullyLigated.QC.R1_trimmed_trimmed_trimmed.fq.gz 2> /dev/null ## trim N6 primer
				trim_galore --no_report_file --cores 6 -a AAAAAAAAAAAAAAAACCTGCAGGNNNNACGAATGCTCTGGCCTNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN tmp.FullyLigated.QC.R1_trimmed_trimmed_trimmed_trimmed.fq.gz 2> /dev/null ### trim oligo-dT primer
				trim_galore --no_report_file --cores 6 -a AAAAAAAAAAAAAA tmp.FullyLigated.QC.R1_trimmed_trimmed_trimmed_trimmed_trimmed.fq.gz 2> /dev/null ## trim polyA tail
				mv tmp.FullyLigated.QC.R1_trimmed_trimmed_trimmed_trimmed_trimmed_trimmed.fq.gz ${prefix}_FullyLigated_trimmed_R1.fastq.gz
				QC=`zcat ${prefix}_FullyLigated_trimmed_R1.fastq.gz | wc -l `
				nQC=`echo "${QC}/4" | bc`
				pQC=`echo "scale=1; 100*${nQC}/${nFL}" | bc`
				echo -e "\nQuality control: ${nQC} (${pQC}%)" | tee -a ${prefix}.Report.txt

				STAR --outSAMunmapped Within --runThreadN ${nthread} --genomeDir ${star_ref} --readFilesIn ${prefix}_FullyLigated_trimmed_R1.fastq.gz --readFilesCommand zcat --outFileNamePrefix tmp.FullyLigated\. --outSAMtype BAM Unsorted
				samtools view -@ ${nthread} -bh -F 260 tmp.FullyLigated.Aligned.out.bam > ${prefix}_FullyLigated_mapped.bam
				samtools view -@ ${nthread} -bh -f 4 tmp.FullyLigated.Aligned.out.bam > ${prefix}_FullyLigated_unmapped.bam
				;;
			*)
			 echo "Unknown data type"; break;;
		esac

		nMapping=`samtools view -@ ${nthread} ${prefix}_FullyLigated_mapped.bam | wc -l `
		pMapping=`echo "scale=1; 100*${nMapping}/${nQC}" | bc`
		echo -e "Mapping: ${nMapping} (${pMapping}%)" | tee -a ${prefix}.Report.txt

#		samtools view -h -@ ${nthread} ${prefix}_FullyLigated_mapped.bam | awk '/^@/ || $3!~/[_,L]/' | samtools view -bh -@ ${nthread} | samtools sort -m 1G -@ ${nthread} > tmp.FullyLigated.Dup.bam
    samtools sort -m 1G -@ ${nthread} ${prefix}_FullyLigated_mapped.bam > tmp.FullyLigated.Dup.bam
    nContig=`samtools view -@ ${nthread} tmp.FullyLigated.Dup.bam | awk '$3=="chrM"' | wc -l `

#     nContig=`samtools view -@ ${nthread} tmp.FullyLigated.Dup.bam | awk '$3~/[_]/ || $3=="chrM"' | wc -l `

#    nContig=`samtools view -@ ${nthread} ${prefix}_FullyLigated_mapped.bam | awk '$3~/[_,L,M]/' | wc -l `
		pContig=`echo "scale=1; 100*${nContig}/${nMapping}" | bc`
		echo -e "Contigs | chrM: ${nContig} (${pContig}%)" | tee -a ${prefix}.Report.txt
 #   nContig=`samtools view -@ ${nthread} ${prefix}_FullyLigated_mapped.bam | awk '$3~/[_,L]/' | wc -l `
		nOthers=`echo "scale=1; ${nMapping}-${nContig}" | bc`

		samtools view -H tmp.FullyLigated.Dup.bam > tmp.header

		samtools view -@ ${nthread} ${prefix}_Barcode.bam | awk '{if($3!="*") print $1, $3}' OFS='\t' > tmp.unsorted.txt
		sort --parallel ${nthread} -k1,1 tmp.unsorted.txt > tmp.BC.txt
		samtools view -@ ${nthread} tmp.FullyLigated.Dup.bam | awk '{print substr($1, 1, length($1)-11), $0}' OFS='\t' > tmp.unsorted.txt
		sort --parallel ${nthread} -k1,1 tmp.unsorted.txt | join -j 1 tmp.BC.txt - | awk '{umi=substr($3,length($3)-9,10)} {print $1":"$2":"umi, $0}' | cut -d ' ' -f 1,5- | sed 's/ /\t/g' | cat tmp.header - | samtools view -@ ${nthread} -bh > ${prefix}_Assign2BC.bam

		nAssigned=`samtools view -@ ${nthread} ${prefix}_Assign2BC.bam | wc -l `
		pAssigned=`echo "scale=1; 100*${nAssigned}/${nOthers}" | bc`
		echo -e "Assigned: ${nAssigned} (${pAssigned}%)" | tee -a ${prefix}.Report.txt
		bedtools bamtobed -i ${prefix}_Assign2BC.bam | awk '{split($4,a,":"); inds=length(a)} {print $1, $2, a[inds], a[inds-3]":"a[inds-2]":"a[inds-1], $6, $4}' OFS='\t' > tmp.unsorted.txt

		case ${library} in
			"DNA")
				#sort --parallel ${nthread} -k1,1 -k2,2n -k3,3 -k4,4 -k5,5 -k6,6 tmp.unsorted.txt | bedtools groupby -i - -g 1,2,3,4 -c 6 -o collapse | sed 's/,.*//g' | cut -f 5 > tmp.unsorted2.txt
        sort --parallel ${nthread} -k4,4 -k3,3 -k1,1 -k2,2n tmp.unsorted.txt | bedtools groupby -i - -g 1,2,3,4 -c 6 -o collapse | sed 's/,.*//g' | cut -f 5 > tmp.unsorted2.txt
				;;
			"RNA")
				#sort --parallel ${nthread} -k1,1 -k2,2n -k3,3 -k4,4 -k5,5 -k6,6 tmp.unsorted.txt | bedtools groupby -i - -g 1,3,4 -c 6 -o collapse | sed 's/,.*//g' | cut -f 4 > tmp.unsorted2.txt
        sort --parallel ${nthread} -k4,4 -k3,3 -k1,1 -k2,2n tmp.unsorted.txt | bedtools groupby -i - -g 1,3,4 -c 6 -o collapse | sed 's/,.*//g' | cut -f 4 > tmp.unsorted2.txt
				;;
			*)
			 echo "Unknown data type"; break;;
		esac
		sort --parallel ${nthread} -k1,1 tmp.unsorted2.txt > tmp.Read1.txt

		nLine=`cat tmp.Read1.txt | wc -l `
		nDup=`echo "scale=1; ${nAssigned}-${nLine}" | bc`
		pDup=`echo "scale=1; 100*${nDup}/${nAssigned}" | bc`
		echo -e "Duplication: ${nDup} (${pDup}%)\n" | tee -a ${prefix}.Report.txt
		
		Rscript /storage/zhangyanxiaoLab/xiongxiong/scripts/paired_tag/analysis/line_reference_duplication.R ${prefix} ${library} ${nRaw} ${pDup}

		samtools view -@ ${nthread} ${prefix}_Assign2BC.bam | grep -w -F -f tmp.Read1.txt - | cat tmp.header - | samtools view -@ ${nthread} -bh | samtools sort -m 1G -@ ${nthread} > ${prefix}_UsefulReads.bam

		nUseful=`samtools view -@ ${nthread} ${prefix}_UsefulReads.bam | wc -l `
		pUseful=`echo "scale=1; 100*${nUseful}/${nRaw}" | bc`
		echo -e "Useful reads: ${nUseful} (${pUseful}%)\n" | tee -a ${prefix}.Report.txt

		awk '{print substr($1,length($1)-18,8)}' tmp.Read1.txt | sed 's/:/\t/g' > tmp.bc
		cut -f 1 tmp.bc | sort --parallel ${nthread} | bedtools groupby -g 1 -c 1 -o count -i - | awk '{print "Barcode3", $0}' OFS='\t' > tmp.bc3
		cut -f 2 tmp.bc | sort --parallel ${nthread} | bedtools groupby -g 1 -c 1 -o count -i - | awk '{print "Barcode2", $0}' OFS='\t' > tmp.bc2
		cut -f 3 tmp.bc | sort --parallel ${nthread} | bedtools groupby -g 1 -c 1 -o count -i - | awk '{print "Barcode1", $0}' OFS='\t' > tmp.bc1
		cat tmp.bc? > ${prefix}_BarcodeUsage.txt

		samtools view -@ ${nthread} ${prefix}_UsefulReads.bam | awk '{print substr($1,length($1)-18,8), $1}' OFS='\t' > tmp.CB.txt
		sort --parallel ${nthread} tmp.CB.txt | bedtools groupby -g 1 -c 2 -o count -i - | sort --parallel ${nthread} -k2,2nr > ${prefix}_CB.Counts

		nLine=`samtools view -@ ${nthread} ${prefix}_UsefulReads.bam | wc -l `
		nFrag=`echo "scale=1; $nLine/1" | bc`
		samtools view -H ${prefix}_UsefulReads.bam | awk '$1 == "@SQ" {OFS="\t";print $2,$3}' - | sed 's/.N://g' > tmp.size
		bedtools genomecov -ibam ${prefix}_UsefulReads.bam -bg | awk '{printf "%s\t%.0f\t%.0f\t%.2f\n",$1,$2,$3,$4*1000000/"'$nFrag'"}' | awk '{$4/=1;print}' OFS='\t' > tmp.unsorted.bdg
		sort --parallel ${nthread} -k1,1 -k2,2n tmp.unsorted.bdg > tmp.bdg
		bedGraphToBigWig tmp.bdg tmp.size ${prefix}_UsefulReads.bw

		pigz -p ${nthread} tmp.FullyLigated.R1.fq tmp.FullyLigated.R2.fq
		mv tmp.FullyLigated.R1.fq.gz ${prefix}_FullyLigated_R1.fastq.gz
		mv tmp.FullyLigated.R2.fq.gz ${prefix}_FullyLigated_R2.fastq.gz
		mv tmp.ReadType.txt ${prefix}_ReadType.txt

		mkdir -p ./process ./bam ./bw
    if [ "$debug" == -1 ]; then
    rm tmp*
    fi 

		mv ${prefix}_FullyLigated_mapped.bam ${prefix}_FullyLigated_unmapped.bam ${prefix}_CB.Counts ${prefix}_BarcodeUsage.txt ${prefix}_ReadType.txt *.gz ${prefix}_Bait.fa ${prefix}_Bait.bam ${prefix}_Barcode.fa ./process
		mv *.bam ./bam
		mv ${prefix}_UsefulReads.bw ./bw
		echo "======== ${prefix} Completed. ========" | tee -a ${prefix}.Report.txt
		cd ..
