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

.libPaths("/GeneCloud003/genecloud/Org_terminal/org_140/terminal/songmm_18513955794/R4.4.1_lib/")

library(readxl)
library(reshape2)


# Figure 1A  突变+CNV+SV瀑布图 -------------------------------------------------

### snv
snv <- read_xlsx("all.snv.filter.recheck.xlsx",sheet = 1)
snv <- snv[snv$recheck=='Y',]
snv1 <- snv[,c("Tumor_Sample_Barcode","Hugo_Symbol","pHGVS","Variant_Classification")]

library(maftools)
clinical <- read.delim("../clinical.group.txt",check.names = F)

pri.maf <- snv[snv$Tumor_Sample_Barcode%in%clinical$sampleid.1021[clinical$group=='Primary'],]

pri.maf <- read.maf(pri.maf,vc_nonSyn = unique(snv$Variant_Classification),isTCGA = F)

laml.titv = titv(maf = pri.maf, plot = FALSE, useSyn = TRUE)
plotTiTv(res = laml.titv)
#计算pvalue
wilcox.test(laml.titv$TiTv.fractions$Ti,laml.titv$TiTv.fractions$Tv)
test <- melt(laml.titv$fraction.contribution)
kruskal.test(test$value ~ test[, "variable"])

meta.maf <- snv[snv$Tumor_Sample_Barcode%in%clinical$sampleid.1021[clinical$group=='Metastases'],]

meta.maf <- read.maf(meta.maf,vc_nonSyn = unique(snv$Variant_Classification),isTCGA = F)

laml.titv = titv(maf = meta.maf, plot = FALSE, useSyn = TRUE)
plotTiTv(res = laml.titv)
#计算pvalue
wilcox.test(laml.titv$TiTv.fractions$Ti,laml.titv$TiTv.fractions$Tv)
test <- melt(laml.titv$fraction.contribution)
kruskal.test(test$value ~ test[, "variable"])


## cnv
cnv <- read.delim("all.cnv.filter.txt")
cnv$sampleid <- gsub("D.*$","",cnv$sampleid)
cnv <- cnv[,c("sampleid","Gene","Status")]
cnv$type <- 'CNV'
colnames(cnv) <- colnames(snv1)
cnv$Variant_Classification <- paste(cnv$Variant_Classification,cnv$pHGVS,sep=' ')
## sv
sv <- read.delim("all.sv.filter.txt",stringsAsFactors = F)
sv$sampleid <- gsub("D.*$","",sv$sampleid)

sv <- dcast(sv,FusionGene~sampleid)
rownames(sv) <- sv$FusionGene
sv <- sv[,-1]
sv <- apply(sv, 2, function(x){x[x!=0] <- "Fusion";x[x!='Fusion'] <- "Not";return(x)})



###合并
all <- rbind(snv1,cnv)
write.table(all,file = 'heatmap.input.txt',quote = F,sep='\t',row.names = F,col.names = F)

mmatx<-matrix("Not",nrow(sv),length(setdiff(all$Tumor_Sample_Barcode,colnames(sv))))
rownames(mmatx)<-rownames(sv)
colnames(mmatx)<-setdiff(all$Tumor_Sample_Barcode,colnames(sv))
sv <- cbind(sv,mmatx)


## msi
MSI <- readRDS('./MSI.rds')
MSI$Tumor_Sample_Barcode <- gsub("D.*$","",MSI$sampleid)

## TMB
TMB <- snv[snv$isTMB=='Y',]
TMB <- as.data.frame(table(TMB$Tumor_Sample_Barcode))
TMB$Tumor_Sample_Barcode <- TMB$Var1

## clinical
clinical <- read.delim("../clinical.group.txt",check.names = F)

clinical <- merge(MSI,clinical,by.y='sampleid.1021',by.x = 'Tumor_Sample_Barcode')
clinical <- merge(clinical,TMB,by='Tumor_Sample_Barcode')

clinical$Stage <- gsub("A|B","",clinical$`Stage at first diagnosis`)
clinical$OS.month <- clinical$`OS (days)(lognumber)`/30
clinical$DFS.month <- clinical$`DFS (days)(lognumber)`/30


###### heatmap 瀑布图
pri.clin <- clinical[clinical$group=='Primary',]
all.pri <- all[all$Tumor_Sample_Barcode %in% pri.clin$Tumor_Sample_Barcode,]
sv.pri <- sv[,pri.clin$Tumor_Sample_Barcode]


p1 <- heatmap.pri(all.pri,sv.pri,pri.clin,top=20,0)

meta.clin <- clinical[clinical$group=='Metastases',]
all.meta <- all[all$Tumor_Sample_Barcode %in% meta.clin$Tumor_Sample_Barcode,]
sv.meta <- sv[,meta.clin$Tumor_Sample_Barcode]


p2 <- heatmap.pri(all.meta,sv.meta,meta.clin,top=20,1)

pdf("../paper1.Figure/heatmap.pdf",width = 12,height = 6)
p1+p2
dev.off()


test <- clinical
rownames(test) <- test$Tumor_Sample_Barcode

test[p1@column_names_param[["labels"]][p1@column_order],"PID"]
test[p2@column_names_param[["labels"]][p2@column_order],"PID"]


heatmap.pri <- function(all,sv,df,top=20,flag=0){
  ###### heatmap 瀑布图
  laml.maf1 <- dcast(all,Hugo_Symbol~Tumor_Sample_Barcode,
                     value.var = 'Variant_Classification',
                     fun.aggregate = function(x){
                       paste(x,collapse = '&',sep = "&")})                        
  
  rownames(laml.maf1) <- laml.maf1$Hugo_Symbol
  laml.maf1 <- laml.maf1[,-1]
  laml.maf1 <- apply(laml.maf1, 2, function(x){x[grep("&",x)] <- "Multi_Hit";return(x)})
  head(laml.maf1)
  
  # sort by gene frequency ---------------------------------------------------------------
  ## 根据基因频率从高到低对基因进行排序
  percentage <- rowSums(laml.maf1!="")
  percentage <- percentage[order(percentage,decreasing = T)]
  laml.maf1 <- laml.maf1[names(percentage),]
  head(laml.maf1[,1:4])
  
  
  ### 保持临床样本编号和突变编号顺序一致
  rownames(df) <- df$Tumor_Sample_Barcode
  df <- df[colnames(laml.maf1),]
  
  ### 保持sv样本编号和突变编号一致
  sv <- as.data.frame(sv[,colnames(laml.maf1)])
  ## 根据基因频率从高到低对基因进行排序
  sv$percentage <- rowSums(sv!="Not")
  sv <- sv[order(sv$percentage,decreasing = T),]
  if(flag==0) {sv <- sv[sv$percentage>=2,]}
  sv <- sv[,-ncol(sv)]
  mutcol = c("missense" = "#0C7094", #偏深蓝
             "nonsense" = "#E21D31", #红
             "frameshift" = "#E2E2E2",
             "cds-del" = "#FFA56C", #偏深绿
             "cds-ins"="#F9BE9D", #灰色
             "cds-indel" = "#FC9054", #橙色
             "splice-3" = "#6E9144", #浅绿
             "splice-5" = "#98B76C", #浅绿
             "span" = "#8FADDF", #蓝色
             "promoter"="#403A84", #紫色
             "CNV gain" = "#ED0000CC", #浅绿
             "CNV loss" = "#0099B4CC", #蓝色
             "Multi_Hit" = '#ADB6B6CC')
  
  
  mut_alter_fun = list( 
    background = function(x, y, w, h) {
      grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = "#FFFFFF", col = "#DFDFDF"))
    },
    "missense" = function(x, y, w, h) {
      grid.rect(x, y,w*0.9, h*0.9, gp = gpar(fill = as.character(mutcol['missense']), col = "#DFDFDF"))
    },
    "nonsense" = function(x, y, w, h) {
      grid.rect(x, y,w*0.9, h*0.9, gp = gpar(fill = as.character(mutcol['nonsense']), col = "#DFDFDF"))
    },
    "frameshift" = function(x, y, w, h) {
      grid.rect(x, y,w*0.9, h*0.9, gp = gpar(fill = as.character(mutcol['frameshift']), col = "#DFDFDF"))
    },
    "cds-del" = function(x, y, w, h) {
      grid.rect(x, y,w*0.9, h*0.9, gp = gpar(fill = as.character(mutcol['cds-del']), col = "#DFDFDF"))
    },
    "cds-ins" = function(x, y, w, h) {
      grid.rect(x, y,w*0.9, h*0.9, gp = gpar(fill = as.character(mutcol['cds-ins']), col = "#DFDFDF"))
    },
    "cds-indel" = function(x, y, w, h) {
      grid.rect(x, y,w*0.9, h*0.9, gp = gpar(fill = as.character(mutcol['cds-indel']), col = "#DFDFDF"))
    },
    "splice-3" = function(x, y, w, h) {
      grid.rect(x, y,w*0.9, h*0.9, gp = gpar(fill = as.character(mutcol['splice-3']), col = "#DFDFDF"))
    },
    "splice-5" = function(x, y, w, h) {
      grid.rect(x, y,w*0.9, h*0.9, gp = gpar(fill = as.character(mutcol['splice-5']), col = "#DFDFDF"))
    },
    "span" = function(x, y, w, h) {
      grid.rect(x, y,w*0.9, h*0.9, gp = gpar(fill = as.character(mutcol['span']), col = "#DFDFDF"))
    },
    "promoter" = function(x, y, w, h) {
      grid.rect(x, y,w*0.9, h*0.9, gp = gpar(fill = as.character(mutcol['promoter']), col = "#DFDFDF"))
    },
    "CNV gain" = function(x, y, w, h) {
      grid.rect(x, y,w*0.9, h*0.9, gp = gpar(fill = as.character(mutcol['CNV gain']), col = "#DFDFDF"))
    },
    "CNV loss" = function(x, y, w, h) {
      grid.rect(x, y,w*0.9, h*0.9, gp = gpar(fill = as.character(mutcol['CNV gain']), col = "#DFDFDF"))
    },
    "Multi_Hit" = function(x, y, w, h) {
      grid.rect(x, y,w*0.9, h*0.9, gp = gpar(fill = as.character(mutcol['Multi_Hit']), col = "#DFDFDF"))
    }
  )
  
  heatmap_legend_param <- list(title = "Mutation",
                               at = names(mutcol),
                               labels = names(mutcol))
  
  ## 注释信息

  
  library(ComplexHeatmap)
  library(circlize)
  
  annocolor <- list(
    Stage = c("IV" = "#C41563","III" = "#CC686C","II" = "#A9ADAA"),
    Gender = c("male" = "#EEA236FF","female" = "#868686FF"),
    Age= colorRamp2(c(30, 60, 80), c("#F7FBFF", "#6BAED6", "#08306B")),
    MSI = c("MSS" = "#DFB3B5","MSI-H" = "#CA676B","/" = '#EDF0F4'),
    Smoke = c("Y" = "#007094","N" = "#EDF0F4"),
    Type = c("Primary" = "#CA7E25","Metastases" = "#C40000"),
    Time = c("synchronous" = "#3182BDFF","asynchronous" = "#E6550DFF","unknown" = "#EDF0F4"),
    Therapy = c("surgery" = "#E28E25","chemotherapy" = "#C30300", "EGFR-TKI" = "#CA676B"),
    OS= colorRamp2(c(5, 30, 70), c("#F7FBFF", "#6BAED6", "#08306B")),
    DFS= colorRamp2(c(1, 10, 40), c("#F7FBFF", "#6BAED6", "#08306B"))
  )
  
  anno.foul <- NULL
  for(i in 1:nrow(sv)){
    annocolor[[rownames(sv)[i]]] <- c('Fusion'='#C40000','Not'="#EDF0F4")
  }

  
  if(flag==0){
    topann_All <- HeatmapAnnotation(
      `ALK:EML4` = as.character(sv["ALK:EML4",]),
      `EML4:ALK` = as.character(sv["EML4:ALK",]),
      `MET:MET` = as.character(sv["MET:MET",]),
      MSI=df$msi_status,
                                    Stage=df$Stage,
                                    Gender=df$Sex,
                                    Age=df$Age,
                                    Smoke=df$Smoking,
                                    Type=df$group,
                                    Time=df$time,
                                    Therapy=df$therapy,
                                    OS=df$OS.month,
                                    DFS=df$DFS.month,
                                    show_annotation_name = T,
                                    annotation_name_gp = gpar(fontsize = 9),
                                    col = annocolor,border = T,
                                    simple_anno_size = unit(3, "mm"),
                                    show_legend = T,
                                    gap = unit(2, "points"),
                                    #gp = gpar(col = "#DFDFDF"),
                                    annotation_name_side = "left"
    )
    
  }else{
    topann_All <- HeatmapAnnotation(
      `ALK:EML4` = as.character(sv["ALK:EML4",]),
      `EML4:ALK` = as.character(sv["EML4:ALK",]),
      `MET:MET` = as.character(sv["MET:MET",]),
      MSI=df$msi_status,
                                    Stage=df$Stage,
                                    Gender=df$Sex,
                                    Age=df$Age,
                                    Smoke=df$Smoking,
                                    Type=df$group,
                                    Time=df$time,
                                    Therapy=df$therapy,
                                    OS=df$OS.month,
                                    DFS=df$DFS.month,
                                    show_annotation_name = T,
                                    annotation_name_gp = gpar(fontsize = 9),
                                    col = annocolor,border = T,
                                    simple_anno_size = unit(3, "mm"),
                                    show_legend = T,
                                    gap = unit(2, "points"),
                                    #gp = gpar(col = "#DFDFDF"),
                                    annotation_name_side = "left"
    )
    
  }
 
  ### 频谱图制作

  ht<-oncoPrint(laml.maf1[1:20,],
                alter_fun = mut_alter_fun, col = mutcol,
                #column_title = column_title,
                bottom_annotation = topann_All, #注释信息在底部
                top_annotation = HeatmapAnnotation(TMB=anno_barplot(df$Freq,gp = gpar(fill = '#A3A2A2',col = '#A3A2A2')),annotation_name_side = "left"),
                remove_empty_columns = FALSE, #去掉空列
                remove_empty_rows = FALSE, #去掉空行
                show_row_names = TRUE,
                row_names_side = "left", #基因位置
                alter_fun_is_vectorized = FALSE,
                show_pct = FALSE,
                pct_side = "right",
                right_annotation =NULL,
                show_column_names = TRUE,
                heatmap_legend_param = heatmap_legend_param,
                show_heatmap_legend = TRUE
  )
  return(ht)
}




# Figure 1B 突变 trunk+private ----------------------------------------------
snv.trunk.priv <- read.delim("../画图需求v1/input/snv.trunk.private.txt")

library(ggpubr)
ggbarplot(snv.trunk.priv,x='PID',y='Freq',position = position_fill(),
          color = 'group',
          fill='group',palette = c("#0C7094","#E21D31","#6E9144"))+theme(axis.text.x=element_text(angle = 90,hjust = 1))

ggsave("../paper1.Figure/snv.trunk&specific.pdf",width = 6,height = 4)



# Figure 1C CNV trunk+branch ----------------------------------------------

cnv.trunk.priv <- read.delim("../画图需求v1/input/cnv.trunk.private.txt")
 
ggbarplot(cnv.trunk.priv,x='PID',y='Freq',position = position_fill(),
          color = 'group',
          fill='group',palette = c("#0C7094","#E21D31","#6E9144"))+theme(axis.text.x=element_text(angle = 90,hjust = 1))

ggsave("../paper1.Figure/cnv.trunk&specific.pdf",width = 6,height = 4)

### 计算R*=Cshare/(Cpt_private+Cmt_private+Cshare)
### The parallel progression model was characterized by a range of 0 < R* ≤ 0.5, while the linear progression model was defined as 0.5 < R* <1.0.
library(dplyr)

R_star <- dcast(cnv.trunk.priv,PID~group,value.var = 'Freq')
R_star$R <- R_star$Trunk/(R_star$Trunk+R_star$`Metastases specific`+R_star$`Primary specific`)

R_star$evolution[R_star$R <= 0.5] <- 'parallel'
R_star$evolution[0.5 < R_star$R] <- 'linear'

table(R_star$evolution)

library(ggpubr)

R_star <- R_star[order(R_star$R,decreasing = F),]
ggbarplot(R_star,x='PID',y='R',fill = 'orange')+
  geom_hline(yintercept = 0.5,linetype='longdash')+
  ylim(0,1)+
  annotate(geom = 'text',x=8,y=0.75,label='Linear Progression',colour='blue')+
  annotate(geom = 'text',x=8,y=0.25,label='Parallel Progression',colour='red')

# Figure 1D cnv trunk & snv trunk correlation -----------------------------

trunk <- read.delim("../画图需求v1/input/snv&cnv.trunk.txt")
label.x=min(trunk[,'SNV'])

ggscatter(trunk, x = 'SNV', y = 'CNV',
          ylab='trunk ratio of CNV',
          xlab='trunk ratio of SNV',
          ggtheme = theme_bw(),
          color = "#E21D31", size = 2, # Points color, shape and size
          add = "reg.line",  # Add regression line
          add.params = list(color = "#0C7094", fill = '#A3A2A2'), # Customize reg. line
          conf.int = TRUE, # Add confidence interval
          cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
          ### 相关性参数设置，参考stat_cor
          cor.coeff.args = list(method = "spearman", 
                                label.x = label.x, 
                                label.sep = "\n")
)+theme(panel.grid = element_blank())


ggsave("../paper1.Figure/cnv&snv.trunk.correlation.pdf",width = 4,height = 4)



# Figure 1E CNV burden ----------------------------------------------------

cnv.burden <- read.delim("../画图需求v1/input/cnv.all.burden.txt")
cnv.gain.burden <- read.delim("../画图需求v1/input/cnv.gain.burden.txt")
cnv.loss.burden <- read.delim("../画图需求v1/input/cnv.loss.burden.txt")

ggpaired(cnv.burden,cond1 = 'Primary',
         palette = c("#0C7094","#E21D31"),
         line.color = 'grey',
         xlab = '',
         ylab = 'CNV burden',
         cond2 = 'Metastases',color = 'condition')+
  stat_compare_means(paired = T)+theme_bw()+theme(panel.grid = element_blank(),legend.position = 'none')

ggsave("../paper1.Figure/all.cnv.burden.pdf",width = 3,height = 3)



ggpaired(cnv.gain.burden,cond1 = 'Primary',
         palette = c("#0C7094","#E21D31"),
         line.color = 'grey',
         xlab = '',
         ylab = 'CNV gain burden',
         cond2 = 'Metastases',color = 'condition')+
  stat_compare_means(paired = T)+theme_bw()+theme(panel.grid = element_blank(),legend.position = 'none')

ggsave("../paper1.Figure/cnv.gain.burden.pdf",width = 3,height = 3)



ggpaired(cnv.loss.burden,cond1 = 'Primary',
         palette = c("#0C7094","#E21D31"),
         line.color = 'grey',
         xlab = '',
         ylab = 'CNV loss burden',
         cond2 = 'Metastases',color = 'condition')+
  stat_compare_means(paired = T)+theme_bw()+theme(panel.grid = element_blank(),legend.position = 'none')

ggsave("../paper1.Figure/cnv.loss.burden.pdf",width = 3,height = 3)



# Figure 1F TMB -----------------------------------------------------------

TMB <- read.delim("../画图需求v1/input/TMB.txt")

ggpaired(TMB,cond1 = 'Primary',
         palette = c("#0C7094","#E21D31"),
         line.color = 'grey',
         xlab = '',
         ylab = 'TMB',
         cond2 = 'Metastases',color = 'condition')+
  stat_compare_means(paired = T)+theme_bw()+theme(panel.grid = element_blank(),legend.position = 'none')

ggsave("../paper1.Figure/TMB.pdf",width = 3,height = 3)


# Figure 1G MSI -----------------------------------------------------------

MSI.freq <- as.data.frame(table(clinical[,c("group","msi_status")]))

ggbarplot(MSI.freq,x='group',y='Freq',position = position_fill(),
          label = T,
          lab.pos = c("in"),
          lab.col = "white",
          fill='msi_status',palette = c("#A3A2A2","#0C7094","#E21D31"),
          )+theme_bw()+theme(panel.grid = element_blank())

ggsave("../paper1.Figure/MSI.pdf",width = 3,height = 3)

# Figure 1G TMB correlation of pri vs meta -----------------------------------------------------------
label.x=min(TMB[,'Metastases'])

ggscatter(TMB, x = 'Primary', y = 'Metastases',
          xlab='TMB of Primary',
          ylab='TMB of Metastases',
          ggtheme = theme_bw(),
          color = "#E21D31", size = 2, # Points color, shape and size
          add = "reg.line",  # Add regression line
          add.params = list(color = "#0C7094", fill = '#A3A2A2'), # Customize reg. line
          conf.int = TRUE, # Add confidence interval
          cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
          ### 相关性参数设置，参考stat_cor
          cor.coeff.args = list(method = "spearman", 
                                label.x = label.x, 
                                label.sep = "\n")
)+theme(panel.grid = element_blank())


ggsave("../paper1.Figure/TMB.correlation.pdf",width = 4,height = 4)


# Figure S1 pathway enrichment --------------------------------------------


library(readxl)
library(tidyr)
library(ggplot2)
library(cowplot)
library(dplyr)
pathway <- as.data.frame(read_excel("DB/pathway.xlsx"))
pathway$Pathway <- as.character(pathway$Pathway)
top10 <- pathway[!grepl("DDR",pathway$Pathway) & !grepl("driver",pathway$Pathway),]
top10 <- top10 %>% separate_rows(Pathway,sep="&")


pathway_enrich <- unique(merge(all,top10,by.x="Hugo_Symbol",by.y="Gene",all=F)[,c("Tumor_Sample_Barcode","Pathway")])

pathway_enrich <- merge(pathway_enrich,clinical,by='Tumor_Sample_Barcode')

pathway_enrich <- unique(pathway_enrich[,c("Tumor_Sample_Barcode","Pathway","group")])

p.p <- pathway_enrich %>% group_by(Pathway,group) %>% dplyr::summarise(freq=n())

ggplot(data=p.p)+
  geom_bar(linewidth=0.8,aes(x=Pathway,y=freq,fill=group,
                             color=group),stat="identity",position = "dodge")+
  scale_color_manual(values = c("#E21D31","#0C7094"))+
  scale_fill_manual(values = c("#E21D3166","#0C709466"))+
  labs(y="Number of Samples")+theme_cowplot()

ggsave('../paper1.Figure/top10pathway.primary_vs_metastases.barplot.pdf',width = 8,height = 4)
