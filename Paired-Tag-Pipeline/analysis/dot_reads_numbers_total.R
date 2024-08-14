##----------------------------
# @Author: Xiong Xiong
# @Date: 2023/1/12
##----------------------------
options(repos="https://mirrors.tuna.tsinghua.edu.cn/CRAN/")
##-----------------------------
library(ggplot2)
library(patchwork)

PlotBC <- function(counTable, dnaCutoff, rnaCutoff){
	
	dnaName <- colnames(counTable)[2]
	rnaName <- colnames(counTable)[3]
	colnames(counTable) <- c('BC', 'CountDNA', 'CountRNA')

	counTable$Type <- 'Fail'
	counTable$Type[which((counTable$CountDNA >= dnaCutoff) & (counTable$CountRNA >= rnaCutoff))] <- 'Pass'
	
	nMedDNA <- paste0('DNA:', gsub(' ', '', format(round(median(counTable$CountDNA[which(counTable$Type=='Pass')])), big.mark=',')))
	nMedRNA <- paste0('RNA:', gsub(' ', '', format(round(median(counTable$CountRNA[which(counTable$Type=='Pass')])), big.mark=',')))

	counTable$CountDNA <- log10(counTable$CountDNA)
	counTable$CountRNA <- log10(counTable$CountRNA)

	nPassFilter <- gsub(' ', '', format(length(which(counTable$Type=='Pass')), big.mark=','))

	figOut <- ggplot(counTable, aes(x=CountDNA, y=CountRNA, color=Type))+geom_point(size=.1)+
				labs(title=paste(dnaName, rnaName, sep=' & '), x='# DNA reads', y='# RNA reads')+
				annotate('text', x=5.5, y=5.5, label=nPassFilter, color='#B22222', size=3)+
				annotate('text', x=6, y=1.4, label='median reads', color='#000000', hjust=1, size=2.5)+
				annotate('text', x=6, y=.8, label=nMedDNA, color='#000000', hjust=1, size=2.5)+
				annotate('text', x=6, y=.2, label=nMedRNA, color='#000000', hjust=1, size=2.5)+
				theme_classic()+scale_color_manual(values=c('#808080', '#B22222'))+
				geom_vline(xintercept=log10(dnaCutoff), color='#000000', linetype='dashed', linewidth=.3)+
				geom_hline(yintercept=log10(rnaCutoff), color='#000000', linetype='dashed', linewidth=.3)+
				scale_x_continuous(limits=c(0, 6), breaks=c(0, log10(dnaCutoff), 6), labels=c('0', as.character(dnaCutoff), expression(10^6)), expand=expansion(mult=c(.05, .05)))+
				scale_y_continuous(limits=c(0, 6), breaks=c(0, log10(rnaCutoff), 6), labels=c('0', as.character(rnaCutoff), expression(10^6)), expand=expansion(mult=c(.05, .05)))+
				theme(plot.margin=unit(c(-.05, .02, 0, .01), 'inches'),
					plot.title=element_text(hjust=.5, vjust=-.9, colour='black', size=12),
					panel.border=element_rect(fill=NA, linewidth=.5, color='black'),
					legend.position='none',
					axis.line.x=element_blank(),
					axis.line.y=element_blank(),
					axis.ticks.x=element_line(linewidth=.3),
					axis.ticks.y=element_line(linewidth=.3),
					axis.title.y=element_text(colour='black', size=10),
					axis.text.x=element_text(colour='black', size=9),
					axis.text.y=element_text(colour='black', size=9, angle=90, hjust=0.5))
	return(figOut)
}

argSet <- commandArgs(trailingOnly = TRUE)
nPairs <- length(argSet) - 2
dnaCutoff <- as.numeric(argSet[1])
rnaCutoff <- as.numeric(argSet[2])

for(i in 1:nPairs){
	
	sampleName <- argSet[i+2]
	fileSets <- unlist(strsplit(sampleName, split=','))
	dnaPrefix <- gsub('/.*/|_CB.Counts', '', fileSets[1])
	rnaPrefix <- gsub('/.*/|_CB.Counts', '', fileSets[2])

	dna <- read.table(fileSets[1])
	colnames(dna) <- c('BC', dnaPrefix)
	rna <- read.table(fileSets[2])
	colnames(rna) <- c('BC', rnaPrefix)

	counTable <- merge.data.frame(dna, rna, by='BC')
	assign(paste0('fig', i), PlotBC(counTable, dnaCutoff, rnaCutoff))
}

figStr <- 'fig1'
figWidth <- 2.4
if(nPairs > 1){
	for(i in 2:nPairs){
		figStr <- paste(figStr, paste0('fig', i), sep='+')
	}
	figWidth <- 4.8
	figStr <- paste(figStr, 'plot_layout(ncol=2)', sep='+')
}
eval(parse(text=figStr))
ggsave('03_Sub-lib_Read_Number.pdf', width=figWidth, height=2.4*ceiling(.5*nPairs))
