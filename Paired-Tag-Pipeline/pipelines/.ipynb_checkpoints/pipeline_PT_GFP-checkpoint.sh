#!/bin/bash
##-------------------------
# @Author: XiongXiong
# @Date: 2023/3/7
##-------------------------

show_usage="Usage: bash pipeline_PT.sh -g <genome> -l <library> -p <thread>\n\n
			-g|--reference-genome \t reference genome, must be one of hg38, mm10 or ce11\n
			-l|--library-type \t library type, must be one of DNA and RNA\n
			-p|--thread \t how many threads to use. Default:8\n
			-h|--help \t\t get help document\n"

GETOPT_ARGS=`getopt -o g:l:p:h -al reference-genome:,library-type:,thread:,help -- "$@"`
eval set -- "$GETOPT_ARGS"

while [ -n "$1" ]
do
	case "$1" in
		-g|--reference-genome) reference_genome=$2; shift 2;;
		-l|--library-type) library=$2; shift 2;;
		-p|--thread) nthread=$2; shift 2;;
		-h|--help) echo -e $show_usage; exit 0;;
		--) break;;
		*) echo -e $1,$2,$show_usage; break;;
	esac
done

if [[ ${library} != "DNA" && ${library} != "RNA" ]]; then
	echo -e "Unknown data type, should be one of DNA or RNA.\n"
	echo -e $show_usage
	exit 0

elif [[ ${reference_genome} != "hg38" && ${reference_genome} != "mm10" && ${reference_genome} != "ce11" ]]; then
	echo -e "Invalid reference genome, should be one of hg38(homo sapiens), mm10(Mus musculus) or ce11(Ceanorhabditis elegans).\n"
	echo -e $show_usage
	exit 0
	
else
	
	if [[ -z ${nthread} ]]; then
		nthread=8
	fi

	rename _001.fastq.gz .fastq.gz ./fastq/*
	for prefix in `ls ./fastq/*_R1.fastq.gz | sed "s/_R1.fastq.gz//g; s/\.\/fastq\///g" | sort --parallel ${nthread} -u `; do
		
		if [ -d ${prefix} ]
		then continue
		fi

		mkdir ./${prefix}
		cd ./${prefix}

		zcat ../fastq/${prefix}'_'R2.fastq.gz | awk '{if(NR%4==1) print $1; else print $0}' > tmp.R2.fq
		awk '{if(NR%4==1) print $1}' tmp.R2.fq > tmp.R2.ReadNames
		awk 'NR%4==2' tmp.R2.fq > tmp.R2.seqs
		paste tmp.R2.ReadNames tmp.R2.seqs > tmp.R2.Reads

		zcat ../fastq/${prefix}'_'R1.fastq.gz | awk '{if(NR%4==1) print $1; else print $0}' > tmp.R1.fq

		echo -e "# reference genome: ${reference_genome}\n# library: ${library}\n# barcode: new barcode 2\n"  | tee ${prefix}.Report.txt
		echo -e "============= ${prefix} =============" | tee -a ${prefix}.Report.txt
		nRaw=`cat tmp.R2.Reads | wc -l`
		pRaw=`echo "scale=1; 100*${nRaw}/${nRaw}" | bc`
		echo "Raw: ${nRaw} (${pRaw}%)" | tee -a ${prefix}.Report.txt

		python /storage/zhangyanxiaoLab/xiongxiong/scripts/paired_tag/scripts/BarcodeExtractor.py tmp.R2.Reads ${library}

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
				nUnload=`awk 'substr($1,1,20)~/CCTGCAGG/' tmp.FullyLigated.R1.fq | wc -l `
				nUnload2=`awk 'substr($1,1,35)~/CTCAAGCACGTGGAT/' tmp.FullyLigated.R1.fq | wc -l `;;
			*)
			 break;;
		esac

		pUnload=`echo "scale=1; 100*${nUnload}/${nFL}" | bc`
		echo -e "+++ Unload: ${nUnload} (${pUnload}%)" | tee -a ${prefix}.Report.txt

		pUnload2=`echo "scale=1; 100*${nUnload2}/${nFL}" | bc`
		echo -e "+++ Unload(no bc1): ${nUnload2} (${pUnload2}%)" | tee -a ${prefix}.Report.txt

		awk '{if($2=="FullyLigated") print ">"$1"\n"$3}' tmp.ReadType.txt | sed 's/@//g' > ${prefix}'_'Barcode.fa
		bowtie /storage/zhangyanxiaoLab/xiongxiong/index/bowtie/new_BC2/PT_Barcode -p ${nthread} --no-unal -v 1 -m 1 --norc -f ${prefix}'_'Barcode.fa -S tmp.BC.sam
		samtools view -@ ${nthread} -bhS tmp.BC.sam > ${prefix}'_'Barcode.bam

		nConta=`awk '$2=="Contamination"' tmp.ReadType.txt | wc -l `
		pConta=`echo "scale=1; 100*${nConta}/${nRaw}" | bc`
		echo "Cross Contamination: ${nConta} (${pConta}%)" | tee -a ${prefix}.Report.txt

		awk '{if($2=="random") print $1}' tmp.ReadType.txt | grep -A 3 -w -F -f - tmp.R1.fq | awk '$0 != "--"' | awk '{if(NR%4==3) print "+"; else print $0}' > tmp.random.R1.fq

		nRandom=`awk '$2=="random"' tmp.ReadType.txt | wc -l `
		pRandom=`echo "scale=1; 100*${nRandom}/${nRaw}" | bc`
		echo "Random: ${nRandom} (${pRandom}%)" | tee -a ${prefix}.Report.txt

		awk '{if($2=="polyG") print $1}' tmp.ReadType.txt | grep -A 3 -w -F -f - tmp.R1.fq | awk '$0 != "--"' | awk '{if(NR%4==3) print "+"; else print $0}' > tmp.polyG.R1.fq

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
			 break;;
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
				trim_galore --no_report_file --cores 6 -a CTGTCTCTTATACA --gzip tmp.FullyLigated.QC.R1_trimmed.fq.gz 2> /dev/null
				trim_galore --no_report_file --cores 6 -a TCGAATGCTCTGGC --gzip tmp.FullyLigated.QC.R1_trimmed_trimmed.fq.gz 2> /dev/null ## trim barcode2 directly ligated with barcode3 + barcode2 + barcode1
				mv tmp.FullyLigated.QC.R1_trimmed_trimmed_trimmed.fq.gz ${prefix}'_'FullyLigated_trimmed_R1.fastq.gz

				QC=`zcat ${prefix}'_'FullyLigated_trimmed_R1.fastq.gz | wc -l `
				nQC=`echo "${QC}/4" | bc`
				pQC=`echo "scale=1; 100*${nQC}/${nFL}" | bc`
				echo -e "\nQuality control: ${nQC} (${pQC}%)" | tee -a ${prefix}.Report.txt

				bowtie2 -p ${nthread} -x /storage/zhangyanxiaoLab/share/bowtie2_index/${reference_genome} -U ${prefix}'_'FullyLigated_trimmed_R1.fastq.gz | samtools view -bh -@ ${nthread} > tmp.${prefix}'_'FullyLigated.bam
				samtools view -@ ${nthread} -bh -F 260 tmp.${prefix}'_'FullyLigated.bam > ${prefix}'_'FullyLigated_mapped.bam
				samtools view -@ ${nthread} -bh -f 4 tmp.${prefix}'_'FullyLigated.bam > ${prefix}'_'FullyLigated_unmapped.bam
				;;
			"RNA")
				trim_galore --no_report_file --cores 6 -a AAAAAAAAAAAAAAAACCTGCAGGNNNNACGAATGCTCTGGCCTNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN tmp.FullyLigated.QC.R1_trimmed.fq.gz 2> /dev/null ### trim oligo-dT primer
				trim_galore --no_report_file --cores 6 -a CCTGCAGGNNNNACGAATGCTCTGGCCTNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN tmp.FullyLigated.QC.R1_trimmed_trimmed.fq.gz 2> /dev/null ## trim N6 primer
				trim_galore --no_report_file --cores 6 -a ACGAATGCTCTGGC tmp.FullyLigated.QC.R1_trimmed_trimmed_trimmed.fq.gz 2> /dev/null ## trim barcode2 directly ligated with barcode3 + barcode2 + barcode1
				mv tmp.FullyLigated.QC.R1_trimmed_trimmed_trimmed_trimmed.fq.gz ${prefix}'_'FullyLigated_trimmed_R1.fastq.gz
				QC=`zcat ${prefix}'_'FullyLigated_trimmed_R1.fastq.gz | wc -l `
				nQC=`echo "${QC}/4" | bc`
				pQC=`echo "scale=1; 100*${nQC}/${nFL}" | bc`
				echo -e "\nQuality control: ${nQC} (${pQC}%)" | tee -a ${prefix}.Report.txt

				STAR --outSAMunmapped Within --runThreadN ${nthread} --genomeDir /storage/zhangyanxiaoLab/qihongjian/datasets/mm10_gfp/mm10 --readFilesIn ${prefix}'_'FullyLigated_trimmed_R1.fastq.gz --readFilesCommand zcat --outFileNamePrefix tmp.FullyLigated\. --outSAMtype BAM Unsorted
				samtools view -@ ${nthread} -bh -F 260 tmp.FullyLigated.Aligned.out.bam > ${prefix}'_'FullyLigated_mapped.bam
				samtools view -@ ${nthread} -bh -f 4 tmp.FullyLigated.Aligned.out.bam > ${prefix}'_'FullyLigated_unmapped.bam
				;;
			*)
			 break;;
		esac

		nMapping=`samtools view -@ ${nthread} ${prefix}'_'FullyLigated_mapped.bam | wc -l `
		pMapping=`echo "scale=1; 100*${nMapping}/${nQC}" | bc`
		echo -e "Mapping: ${nMapping} (${pMapping}%)" | tee -a ${prefix}.Report.txt

		samtools view -h -@ ${nthread} ${prefix}'_'FullyLigated_mapped.bam | awk '/^@/ || $3!~/[_,L,M,Y]/' | samtools view -bh -@ ${nthread} | samtools sort -m 4G -@ ${nthread} > tmp.FullyLigated.Dup.bam
		nQ10=`samtools view -@ ${nthread} tmp.FullyLigated.Dup.bam | wc -l `
		nOthers=`echo "scale=1; ${nMapping}-${nQ10}" | bc`
		pOthers=`echo "scale=1; 100*${nOthers}/${nMapping}" | bc`
		echo -e "Contigs | chrM: ${nOthers} (${pOthers}%)" | tee -a ${prefix}.Report.txt

		samtools view -H tmp.FullyLigated.Dup.bam > tmp.header

		samtools view -@ ${nthread} tmp.BC.sam | awk '{if($3!="*") print $1, $3}' OFS='\t' > tmp.unsorted.txt
		sort --parallel ${nthread} -k1,1 tmp.unsorted.txt > tmp.BC.txt
		samtools view -@ ${nthread} tmp.FullyLigated.Dup.bam | awk '{print substr($1, 1, length($1)-11), $0}' OFS='\t' > tmp.unsorted.txt
		sort --parallel ${nthread} -k1,1 tmp.unsorted.txt | join -j 1 tmp.BC.txt - | awk '{umi=substr($3,length($3)-9,10)} {print $1":"$2":"umi, $0}' | cut -d ' ' -f 1,5- | sed 's/ /\t/g' | cat tmp.header - | samtools view -@ ${nthread} -bh > ${prefix}_Assign2BC.bam

		nAssigned=`samtools view -@ ${nthread} ${prefix}_Assign2BC.bam | wc -l `
		pAssigned=`echo "scale=1; 100*${nAssigned}/${nQ10}" | bc`
		echo -e "Assigned: ${nAssigned} (${pAssigned}%)" | tee -a ${prefix}.Report.txt
		bedtools bamtobed -i ${prefix}_Assign2BC.bam | awk '{split($4,a,":"); inds=length(a)} {print $1, $2, a[inds], a[inds-3]":"a[inds-2]":"a[inds-1], $6, $4}' OFS='\t' > tmp.unsorted.txt

		case ${library} in
			"DNA")
				sort --parallel ${nthread} -k1,1 -k2,2n -k3,3 -k4,4 -k5,5 -k6,6 tmp.unsorted.txt | bedtools groupby -i - -g 1,2,3,4 -c 6 -o collapse | sed 's/,.*//g' | cut -f 5 > tmp.unsorted2.txt
				;;
			"RNA")
				sort --parallel ${nthread} -k1,1 -k2,2n -k3,3 -k4,4 -k5,5 -k6,6 tmp.unsorted.txt | bedtools groupby -i - -g 1,3,4 -c 6 -o collapse | sed 's/,.*//g' | cut -f 4 > tmp.unsorted2.txt
				;;
			*)
			 break;;
		esac
		sort --parallel ${nthread} -k1,1 tmp.unsorted2.txt > tmp.Read1.txt

		nLine=`cat tmp.Read1.txt | wc -l `
		nDup=`echo "scale=1; ${nAssigned}-${nLine}" | bc`
		pDup=`echo "scale=1; 100*${nDup}/${nAssigned}" | bc`
		echo -e "Duplication: ${nDup} (${pDup}%)\n" | tee -a ${prefix}.Report.txt
		
		Rscript /storage/zhangyanxiaoLab/xiongxiong/scripts/paired_tag/analysis/line_reference_duplication.R ${prefix} ${library} ${nRaw} ${pDup}

		samtools view -@ ${nthread} ${prefix}_Assign2BC.bam | grep -w -F -f tmp.Read1.txt - | cat tmp.header - | samtools view -@ ${nthread} -bh | samtools sort -m 4G -@ ${nthread} > ${prefix}_UsefulReads.bam
		
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

		pigz -p ${nthread} tmp.FullyLigated.R1.fq tmp.FullyLigated.R2.fq tmp.random.R1.fq tmp.polyG.R1.fq
		mv tmp.FullyLigated.R1.fq.gz ${prefix}'_'FullyLigated_R1.fastq.gz
		mv tmp.FullyLigated.R2.fq.gz ${prefix}'_'FullyLigated_R2.fastq.gz
		mv tmp.random.R1.fq.gz ${prefix}'_'random_R1.fastq.gz
		mv tmp.polyG.R1.fq.gz ${prefix}'_'polyG_R1.fastq.gz
		mv tmp.ReadType.txt ${prefix}'_'ReadType.txt

		mkdir -p ./process ./bam ./bw

		rm tmp*
		mv ${prefix}'_'FullyLigated_mapped.bam ${prefix}'_'FullyLigated_unmapped.bam ${prefix}_CB.Counts ${prefix}_BarcodeUsage.txt ${prefix}'_'ReadType.txt *.gz ${prefix}'_'Barcode.fa ./process
		mv *.bam ./bam
		mv ${prefix}_UsefulReads.bw ./bw
		echo "======== ${prefix} Completed. ========" | tee -a ${prefix}.Report.txt
		cd ..
	done
	echo "All samples completed."
fi
