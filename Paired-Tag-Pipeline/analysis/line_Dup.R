##----------------------------
# @Author: Xiong Xiong
# @Date: 2023/2/7
##----------------------------
options(repos="https://mirrors.tuna.tsinghua.edu.cn/CRAN/")
##-----------------------------
library(ggplot2)
library(patchwork)

PlotDup <- function(dupSet, libtype){
	
	plotSet <- dupSet[dupSet$Lib == libtype, ]
	
	figOut <- ggplot(plotSet, aes(x=Count, y=Dup, color=Name))+geom_line()+
		theme_classic()+labs(title=NULL, x='Assigned reads (million)', y='Duplication (%)')+
		scale_x_continuous(limits=c(0, 86), breaks=seq(0, 80, 20), expand=expansion(mult=c(.006, 0)))+
		scale_y_continuous(limits=c(0, 86), breaks=seq(0, 80, 20), expand=c(0, 0))+
		labs(title=libtype)+
		scale_color_manual(values=c('#F19135', '#304872'))+
		theme(plot.margin=unit(c(0, .01, 0, .03), 'inches'),
			plot.title=element_text(size=12, hjust=.5, vjust=-.5),
			panel.grid.major.x=element_blank(),
			panel.grid.major.y=element_line(color='#C0C0C0', linewidth=.5),
			panel.border=element_rect(fill=NA, linewidth=.75, color='black'),
			legend.background=element_blank(),
			legend.position=c(.5, .1),
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
	
	return(figOut)
}

load('/storage/zhangyanxiaoLab/xiongxiong/share/Paired-Tag/ZCX_Sub-lib_Duplication.RData')

dupFile <- read.table('Downsample_Dup.txt', header=T)

dupSet <- rbind(dupRef, dupFile)
dupSet$Count <- dupSet$Count/1000000

figDNA <- PlotDup(dupSet, 'DNA')
figRNA <- PlotDup(dupSet, 'RNA')

figDNA+figRNA+plot_layout(ncol=2)
ggsave('downsample_duplication.pdf', width=5.5, height=3)
