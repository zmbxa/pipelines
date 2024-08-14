##----------------------------
# @Author: Xiong Xiong
# @Date: 2023/2/7
##----------------------------
options(repos="https://mirrors.tuna.tsinghua.edu.cn/CRAN/")
##-----------------------------
library(ggplot2)

argSet <- commandArgs(trailingOnly=TRUE)

load('/storage/zhangyanxiaoLab/xiongxiong/share/Paired-Tag/Dupliction_GSE152020.RData')
load('/storage/zhangyanxiaoLab/xiongxiong/share/others/colorset.RData')

libLabel <- argSet[1]
libtype <- argSet[2]
nRaw <- as.numeric(argSet[3])/1000000
pDup <- as.numeric(argSet[4])

plotSet <- dupRefAll[dupRefAll$Lib == libtype, ]

p <- ggplot(plotSet, aes(x=Count, y=Dup, color=Name))+geom_line()+
		geom_point(aes(nRaw, pDup), shape=23, fill='#32CD32', color='#32CD32')+
		annotate('rect', xmin=nRaw-50, xmax=nRaw+50, ymin=pDup-10, ymax=pDup-4, fill='#FFFFFF')+
		annotate('text', x=nRaw, y=pDup-7, label=libLabel, color='#32CD32')+
		theme_classic()+labs(title=NULL, x='Sequencing depth (million)', y='Duplication (%)')+
		scale_x_continuous(limits=c(0, 370), breaks=seq(0, 360, 120), expand=expansion(mult=c(.006, 0)))+
		scale_y_continuous(limits=c(0, 96), breaks=seq(0, 90, 30), expand=c(0, 0))+
		labs(title=libtype)+
		scale_color_manual(values=colorSet[1:length(unique(plotSet$Name))])+
		theme(plot.margin=unit(c(0, .08, .01, .03), 'inches'),
			plot.title=element_text(size=12, hjust=.5, vjust=-.5),
			panel.grid.major.x=element_blank(),
			panel.grid.major.y=element_line(color='#C0C0C0', linewidth=.5),
			panel.border=element_rect(fill=NA, linewidth=.75, color='black'),
			legend.background=element_blank(),
			legend.position='none',
			legend.key=element_blank(),
			legend.title=element_blank(),
			legend.box="",
			legend.text=element_text(size=11),
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
ggsave(paste0(libLabel, '_duplication.pdf'), p, width=3, height=3)

if(nRaw<=10){
	plotSet <- dupRefAll[dupRefAll$Lib == libtype & dupRefAll$Count<=12, ]

	p <- ggplot(plotSet, aes(x=Count, y=Dup, color=Name))+geom_line()+
			geom_point(aes(nRaw, pDup), shape=23, fill='#32CD32', color='#32CD32')+
			annotate('text', x=nRaw, y=pDup-7, label=libLabel, color='#32CD32')+
			theme_classic()+labs(title=NULL, x='Sequencing depth (million)', y='Duplication (%)')+
			scale_x_continuous(limits=c(0, 12.5), breaks=seq(0, 12, 3), expand=expansion(mult=c(.006, 0)))+
			scale_y_continuous(limits=c(0, 96), breaks=seq(0, 90, 30), expand=c(0, 0))+
			labs(title=libtype)+
			scale_color_manual(values=colorSet[1:length(unique(plotSet$Name))])+
			theme(plot.margin=unit(c(0, .08, .01, .03), 'inches'),
				plot.title=element_text(size=12, hjust=.5, vjust=-.5),
				panel.grid.major.x=element_blank(),
				panel.grid.major.y=element_line(color='#C0C0C0', linewidth=.5),
				panel.border=element_rect(fill=NA, linewidth=.75, color='black'),
				legend.background=element_blank(),
				legend.position='none',
				legend.key=element_blank(),
				legend.title=element_blank(),
				legend.box="",
				legend.text=element_text(size=11),
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
	ggsave(paste0(libLabel, '_duplication_shallow.pdf'), p, width=3, height=3)
}

