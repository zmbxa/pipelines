#!/usr/bin/Rscript
##----------------------------
# @Author: Xiong Xiong
# @Date: 2023/1/11
##----------------------------
options(repos="https://mirrors.tuna.tsinghua.edu.cn/CRAN/")
##-----------------------------
library(ggplot2)

argSet <- commandArgs(trailingOnly = TRUE)
load('/storage/zhangyanxiaoLab/xiongxiong/share/Paired-Tag/ZCX_ZYX_BarcodePercent.RData')
cbCount <- data.frame('Barcode'=character(0),
					  'Count'=numeric(0),
					  'CumPerc'=numeric(0),
					  'Rank'=numeric(0),
					  'Lib'=character(0))

sampleName = argSet[1]
print(sampleName)
prefix <- strsplit(sampleName, "/")[[1]][1] 
tmpCount <- read.table(sampleName)
colnames(tmpCount) <- c('Barcode', 'Count')
tmpCount$CumPerc <- 100*cumsum(tmpCount$Count)/sum(tmpCount$Count)
tmpCount$Rank <- (1:nrow(tmpCount))/1000
tmpCount$Lib <- prefix
cbCount <- rbind(cbCount, tmpCount)

plotData <- rbind(compData, cbCount)

p <- ggplot(plotData, aes(x=Rank, y=CumPerc, color=Lib))+geom_line()+
	geom_vline(xintercept=c(2.5, 50), color='#708090', linetype='dashed')+
	theme_classic()+labs(title=NULL, x='Ranked barcode', y='Cumulative reads proportion (%)')+
	scale_color_manual(values=colorSet[1:length(unique(plotData$Lib))])+
	scale_x_continuous(limits=c(-10, 120), breaks=c(2.5, 50, 100), labels=c('2.5k', '50k', '100k'), expand=c(0, 0))+
	scale_y_continuous(limits=c(0, 105), breaks=seq(0, 100, 25), expand=c(0, 0))+
	theme(plot.margin=unit(c(.15, .01, 0, .03), 'inches'),
		plot.title=element_blank(),
		panel.grid.major.x=element_blank(),
		panel.grid.major.y=element_line(color='#C0C0C0', linewidth=.5),
		panel.border=element_rect(fill=NA, linewidth=.75, color='black'),
		legend.background=element_blank(),
		legend.position=c(.8, .5),
		legend.key=element_blank(),
		legend.key.height=unit(.13, 'inches'),
		legend.title=element_blank(),
		legend.box='',
		legend.text=element_text(size=9),
		legend.spacing.x=unit(.01, 'inches'),
		legend.spacing.y=unit(0, 'inches'),
		axis.ticks.x=element_blank(),
		axis.ticks.y=element_blank(),
		axis.line.x=element_blank(),
		axis.line.y=element_blank(),
		axis.title.x=element_text(size=12),
		axis.title.y=element_text(size=12),
		axis.text.x=element_text(size=12, colour='black'),
		axis.text.y=element_text(size=12, colour='black', angle=90, hjust=0.5))

ggsave(paste0("qc/", prefix, "_barcode_rank.pdf"), width=8, height=8)

