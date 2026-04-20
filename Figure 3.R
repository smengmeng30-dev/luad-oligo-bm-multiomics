#!/usr/bin/env Rscript

#' @title 开发规范
#' @description
#' 对于每个脚本，希望有明确的帮助文档以及代码注释
#' 帮助文档这里推荐[docopt](https://github.com/docopt/docopt.R)
#' @description
#' 1. 对于R包，默认的repos为[RSPM](https://packagemanager.rstudio.com)，可通过`getOption("repos")`查看；
#' 2. 默认情况下，我们仅维护每个 R 包的其中一个版本
#' `.libPaths()`查看 R 包路径；
#' 3. 对于需要自己维护的 R 包，可以在 terminal 找一个可写路径，
#' 使用`.libPaths(new = c("./",.libPaths()))`定义当前路径为最高优先级
#' 之后即可使用`install.packages()`下载至当前目录
#' 在下一次启动 R 会话时使用`.libPaths()`进行导入
#' @description
#' 对于每个 R 文件，建议写完后使用`styler::style_file()`格式化
#' @description
#' 默认的工作路径为每个用户 terminal 主目录
setwd("/GeneCloud003/genecloud/Org_terminal/org_140/terminal/songmm_18513955794/project/yzz_brainMeta/1021/")


library(ComplexHeatmap)
library(data.table)
# snv+cnv tree map --------------------------------------------------------

### snv+cnv
all.vari <- read.delim("../画图需求v1/input/all.treeheatmap.txt",header = T)



for(i in unique(all.vari$PID)){
  tmp <- all.vari[all.vari$PID==i,]
  tmp1 <- dcast(tmp,PID+flag~group,fun.aggregate = length)
  tmp1 <- tmp1[order(tmp1$Primary,tmp1$Metastases),]
  mat <- tmp1
  mat$p.color[mat$Primary!=0 & mat$Metastases!=0] <- 'red'
  mat$m.color[mat$Primary!=0 & mat$Metastases!=0] <- 'red'
  
  mat$p.color[mat$Primary!=0 & mat$Metastases==0] <- 'blue'
  mat$m.color[mat$Primary!=0 & mat$Metastases==0] <- 'grey'
  
  mat$p.color[mat$Primary==0 & mat$Metastases!=0] <- 'grey'
  mat$m.color[mat$Primary==0 & mat$Metastases!=0] <- 'orange'
  
  
  rownames(mat) <- mat$flag
  mat <- mat[,c(-1,-2,-3,-4)]
  mat <- as.matrix(mat)
  colors <- c("#0C7094","#E21D31","#6E9144","#EFEFEF")
  names(colors) <- c("red",'blue','orange',"grey")
  filename <- paste("../paper1.Figure/tree/",i,".heatmap.pdf",sep="")
  pdf(file = filename,width=4,height=8)
  p <- Heatmap(mat,
               show_row_names = F,
               border = F,
               col = colors,
               column_gap = unit(1, "mm"),
               heatmap_width = unit(1,"cm"),heatmap_height = unit(10,"cm"),
               column_split  = 1:ncol(mat),
               show_heatmap_legend = F,
               name='Mutation',cluster_rows = F,cluster_columns = T)
  print(p)
  dev.off()
  
  
  trunk_clonal <- rownames(mat)[which(mat[,1]=='red' & mat[,2]=='red')]
  trunk_clonal_driver <- unique(trunk_clonal)
  
  p_subclonal <- rownames(mat)[which(mat[,1]=='blue' & mat[,2]=='grey')]
  p_subclonal_driver <- unique(p_subclonal)
  
  m_subclonal <- rownames(mat)[which(mat[,1]=='grey' & mat[,2]=='orange')]
  m_subclonal_driver <- unique(m_subclonal)
  
  trunk_clonal_driver <- paste(trunk_clonal_driver,collapse = "\n")
  p_subclonal_driver <- paste(p_subclonal_driver,collapse = "\n")
  m_subclonal_driver <- paste(m_subclonal_driver,collapse = "\n")
  
  dat <- data.frame(num=c(length(trunk_clonal),
                          length(p_subclonal),length(m_subclonal)),length=0)
  
  col <- which(dat$num==max(dat$num))
  dat$length[col] <- 10
  dat$length <- dat$num*dat$length[col]/dat$num[col]
  dat$length <- ifelse(dat$length<1,1,dat$length)
  
  pdf(file = paste("../paper1.Figure/tree/",i,".tree.pdf",sep=""),width=8,height=10)
  plot.new()
  plot.window(c(-20,20),c(-20,20))
  lwd <- 3
  
  trunk_clonal_y <- dat$length[1]
  text_y <- dat$length[1]/2
  lines(x=c(0,0),y=c(0,trunk_clonal_y),type="l",col="#0C7094",lwd=lwd)
  #text(x=0.1,y=text_y,labels = trunk_clonal_driver,cex=0.7,col="#225EA8",adj=0,font=2)
  text(x=0.2,y=trunk_clonal_y,labels = paste("(",dat$num[1],")"),cex=0.7,col="#8E8A8A",adj=0,font=2)
  
  
  N_x <- dat$length[2]/2 
  N_y <- (sqrt(3)*dat$length[2]/2+trunk_clonal_y)
  lines(x=c(0,N_x),y=c(trunk_clonal_y,N_y),type="l",col="#E21D31",lwd=lwd)
  #text(x=N_x/2,y=N_y/2,labels = p_subclonal_driver,cex=0.7,col="#F39B7FFF",adj=0,font=2,srt = 30)
  text(x=N_x,y=N_y,labels = "P",cex=0.7,col="black",adj=0,font=2)
  text(x=N_x-0.2,y=N_y,labels = paste("(",dat$num[2],")"),cex=0.7,col="#8E8A8A",adj=0,font=2)
  
  S_x <- -dat$length[3]/2
  S_y <- (sqrt(3)*dat$length[3]/2+trunk_clonal_y)
  lines(x=c(0,S_x),y=c(trunk_clonal_y,S_y),type="l",col="#6E9144",lwd=lwd)
  #text(x=S_x/2,y=S_y/2,labels = m_subclonal_driver,cex=0.7,col="#91D1C2FF",adj=0,font=2,srt = 30)
  text(x=S_x,y=S_y,labels = "M",cex=0.7,col="black",adj=0,font=2)
  text(x=S_x+0.2,y=S_y,labels = paste("(",dat$num[3],")"),cex=0.7,col="#8E8A8A",adj=0,font=2)
  dev.off()
}




# trunk vs private --------------------------------------------------------
trunk.priv <- dcast(all.vari,flag+PID~group)
trunk.priv$group <- ifelse(is.na(trunk.priv$Metastases) & !is.na(trunk.priv$Primary),"Primary",ifelse(!is.na(trunk.priv$Metastases) & is.na(trunk.priv$Primary),"Metastases","Trunk"))

trunk.priv <- as.data.frame(table(trunk.priv[,c("PID","group")]))

trunk.priv <- trunk.priv %>%
   group_by(PID) %>%
  mutate(total=sum(Freq))

trunk.priv$percentage <- trunk.priv$Freq/trunk.priv$total

library(ggpubr)

ggboxplot(trunk.priv,x='group',y='percentage',
          add='jitter',
          xlab = '',
          color='group',palette = c("#0C7094","#E21D31","#6E9144"))+theme_bw()+theme(legend.position = 'none',panel.grid = element_blank())+stat_compare_means(comparisons = list(c("Trunk","Primary"),c("Trunk","Metastases"),c("Primary","Metastases")))

ggsave("../paper1.Figure/trunk&specific.boxplot.pdf",width = 4,height = 4)


