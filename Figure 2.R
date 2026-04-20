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


# Figure 2A ccf plot ------------------------------------------------------
clinical <- read.delim("../clinical.group.txt")
ccf <- read.delim("../画图需求v1/input/all.ccf.txt")
ccf$cluster <- paste("cluster",ccf$CitupCluster,sep='')

ggscatter(ccf,x='primary',y='metastases',
          fill = 'cluster',shape = 21,label = "Gene",label.select = list(criteria = "`primary` > 50 | `metastases` > 50"),palette = 'lancet',repel=T,font.label=8)+
  geom_abline(intercept = 0,slope = 1,linetype='dashed')+
  theme_bw()+theme(panel.grid = element_blank())+facet_wrap(~PID)


ggsave("../paper1.Figure/CCF.plot.pdf",width = 9,height = 7)


ggscatter(ccf,x='primary',y='metastases',shape = 21,label = "Gene",label.select = list(criteria = "`primary` > 50 | `metastases` > 50"),
          palette = 'lancet',repel=T,font.label=8)+
  geom_abline(intercept = 0,slope = 1,linetype='dashed')+
  theme_bw()+theme(panel.grid = element_blank())+stat_cor()


ggsave("../paper1.Figure/CCF.plot.combine.pdf",width = 6,height = 6)


# Figure 2B clone number --------------------------------------------------

clone <- read.delim("../画图需求v1/input/clone.number.txt")

ggpaired(clone,cond1 = 'p.cloneNumber',
         palette = c("#0C7094","#E21D31"),
         line.color = 'grey',
         xlab = '',
         ylab = 'clone number',
         cond2 = 'm.cloneNumber',color = 'condition')+
  stat_compare_means(paired = T)+theme_bw()+theme(panel.grid = element_blank(),legend.position = 'none')

ggsave("../paper1.Figure/m_vs_p.clone.number.pdf",width = 3,height = 3)


# Figure 2C ITH score --------------------------------------------------


library(readxl)
library(maftools)

snv <- read_excel("all.snv.filter.recheck.xlsx")
snv <- snv[snv$recheck=='Y',]

maf <- read.maf(snv,vc_nonSyn = unique(snv$Variant_Classification))

#ITH <- math.score(maf,vafCol = 'caseAF',vafCutOff = 0)

Heterogeneity_result = inferHeterogeneity(maf = maf,tsb = unique(snv$Tumor_Sample_Barcode),vafCol = 'caseAF')
Heterogeneity_result2 = Heterogeneity_result$clusterData
Heterogeneity_result3 = Heterogeneity_result2[!duplicated(Heterogeneity_result2$Tumor_Sample_Barcode),c("Tumor_Sample_Barcode","MATH")]


Heterogeneity_result3 <- merge(Heterogeneity_result3,clinical,by='Tumor_Sample_Barcode')


library(ggpubr)

Heterogeneity_result4 <- dcast(Heterogeneity_result3,PID~group,value.var = 'MATH')

ggpaired(Heterogeneity_result4,cond1 = 'Primary',cond2 = 'Metastases',
         color = 'condition',line.color = 'grey',
         ylab = 'ITH score',xlab='')+stat_compare_means(paired = T)+
  scale_color_manual(values = c("#0C7094","#E21D31"))+theme_bw()+theme(legend.position = 'none',panel.grid = element_blank())


ggsave("../paper1.Figure/ITH.pdf",width = 3,height = 3)

# Figure 2D clonal和subclonal分布 --------------------------------------------------
library(ggpubr)
library(data.table)
clonal_subclonal <- read.delim("../画图需求v2/clonal_subclonal.txt")

clonal_subclonal <- melt(clonal_subclonal,measure.vars=6:7)
clonal_subclonal$sample <- paste(clonal_subclonal$PID,clonal_subclonal$variable,sep="_")
clonal_subclonal$evo1 <- ifelse(clonal_subclonal$evo %in% c("trunk clonal"),"clonal","subclonal")
clonal_subclonal <- clonal_subclonal[!is.na(clonal_subclonal$value),]
clonal_subclonal$evo2 <- ifelse(as.numeric(gsub("%","",clonal_subclonal$value))>50,'clonal','subclonal')

total <- unique(clonal_subclonal[,c("sample","Gene",'evo1')])

total <- as.data.frame(table(total[,c("Gene","evo1")]))
summ <- total %>%
  group_by(Gene) %>%
  mutate(sum=sum(Freq))

summ <- unique(summ[order(summ$sum,decreasing = T),c("Gene","sum")])
show_gene <- summ$Gene[summ$sum>=5]

show_gene_bar <- total[total$Gene %in% as.character(show_gene),]
show_gene_bar$Gene <- factor(show_gene_bar$Gene,levels = as.character(show_gene))
show_gene_bar <- show_gene_bar[order(show_gene_bar$Gene),]

ggbarplot(show_gene_bar,x='Gene',y='Freq',fill='evo1',color = 'white',palette = c("#0C7094","#E21D31"),position = position_dodge(),xlab = '')+theme_bw()+theme(panel.grid = element_blank())

ggsave("../paper1.Figure/clonal_subclonal.distribution.pdf",width = 8,height = 3)



total.pri <- unique(clonal_subclonal[clonal_subclonal$variable=='Primary',c("sample","Gene",'evo1')])

total.pri <- as.data.frame(table(total.pri[,c("Gene","evo1")]))
summ.pri <- total.pri %>%
  group_by(Gene) %>%
  mutate(sum=sum(Freq))

summ.pri <- unique(summ.pri[order(summ.pri$sum,decreasing = T),c("Gene","sum")])
show_gene.pri <- summ.pri$Gene[summ.pri$sum>=3]

show_gene_bar.pri <- total.pri[total.pri$Gene %in% as.character(show_gene.pri),]
show_gene_bar.pri$Gene <- factor(show_gene_bar.pri$Gene,levels = as.character(show_gene.pri))
show_gene_bar.pri <- show_gene_bar.pri[order(show_gene_bar.pri$Gene),]

ggbarplot(show_gene_bar.pri,x='Gene',y='Freq',fill='evo1',color = 'white',palette = c("#0C7094","#E21D31"),position = position_dodge(),xlab = '')+theme_bw()+theme(panel.grid = element_blank())

ggsave("../paper1.Figure/pri.clonal_subclonal.distribution.pdf",width = 8,height = 3)


total.meta <- unique(clonal_subclonal[clonal_subclonal$variable=='Metastases',c("sample","Gene",'evo1')])

total.meta <- as.data.frame(table(total.meta[,c("Gene","evo1")]))
summ.meta <- total.meta %>%
  group_by(Gene) %>%
  mutate(sum=sum(Freq))

summ.meta <- unique(summ.meta[order(summ.meta$sum,decreasing = T),c("Gene","sum")])
show_gene.meta <- summ.meta$Gene[summ.meta$sum>=3]

show_gene_bar.meta <- total.meta[total.meta$Gene %in% as.character(show_gene.meta),]
show_gene_bar.meta$Gene <- factor(show_gene_bar.meta$Gene,levels = as.character(show_gene.meta))
show_gene_bar.meta <- show_gene_bar.meta[order(show_gene_bar.meta$Gene),]

ggbarplot(show_gene_bar.meta,x='Gene',y='Freq',fill='evo1',color = 'white',palette = c("#0C7094","#E21D31"),position = position_dodge(),xlab = '')+theme_bw()+theme(panel.grid = element_blank())

ggsave("../paper1.Figure/meta.clonal_subclonal.distribution.pdf",width = 8,height = 3)


write.table(show_gene_bar,file = '../paper1.Figure/clonal_subclonal.stat.txt',quote = F,sep='\t',row.names = F)

write.table(show_gene_bar.pri,file = '../paper1.Figure/pri.clonal_subclonal.stat.txt',quote = F,sep='\t',row.names = F)

write.table(show_gene_bar.meta,file = '../paper1.Figure/meta.clonal_subclonal.stat.txt',quote = F,sep='\t',row.names = F)


total <- unique(clonal_subclonal[,c("sample","Gene",'evo1')])
total <- as.data.frame(table(total[,c("Gene","evo1")]))

clonal_num <- length(unique(clonal_subclonal$sample[clonal_subclonal$evo1=='clonal']))
subclonal_num <- length(unique(clonal_subclonal$sample[clonal_subclonal$evo1=='subclonal']))


total1 <- dcast(total,Gene~evo1,value.var = 'Freq')

res <- NULL

for(i in 1:nrow(total1)){
  gene <- total1[i,'Gene']
  clonal <- total1[i,'clonal']
  subclonal <- total1[i,'subclonal']
  others_clonal <- clonal_num-total1[i,'clonal']
  others_subclonal <- subclonal_num-total1[i,'subclonal']
  
  ######## asyn的原发 vs 脑转 （clonal vs subclonal）
  tmp <- fisher.test(matrix(c(clonal,subclonal,others_clonal,others_subclonal),nrow = 2,byrow = T))
  
  res <- rbind(res,data.frame(gene=gene,
                              clonal=clonal,
                              others_clonal=others_clonal,
                              subclonal=subclonal,
                              others_subclonal=others_subclonal,
                              pvalue=tmp$p.value,
                              OR=as.numeric(tmp$estimate)))
  
  
}


write.table(res,file = '../paper1.Figure/clonal_subclonal.fisher.txt',quote = FALSE,sep='\t',row.names = FALSE)





total <- unique(clonal_subclonal[,c("sample","Gene",'evo1','variable')])
total <- as.data.frame(table(total[,c("Gene","evo1",'variable')]))
summ <- total %>%
  group_by(Gene) %>%
  mutate(sum=sum(Freq))

summ <- unique(summ[order(summ$sum,decreasing = T),c("Gene","sum")])
show_gene <- summ$Gene[summ$sum>=4]

show_gene_bar <- total[total$Gene %in% as.character(show_gene),]
show_gene_bar$Gene <- factor(show_gene_bar$Gene,levels = as.character(show_gene))
show_gene_bar <- show_gene_bar[order(show_gene_bar$Gene),]


shown_top <- dcast(show_gene_bar,Gene~variable+evo1,value.var = 'Freq')
shown_top[,2:3] <- shown_top[,2:3]/rowSums(shown_top[,2:3])
shown_top[,4:5] <- shown_top[,4:5]/rowSums(shown_top[,4:5])


test <- unname(t(shown_top[,2:3]))
colnames(test) <- shown_top$Gene

pdf("../paper1.Figure/clonal_subclonal.distribution.pdf",width = 8,height = 4)

bp <- barplot(test, width = 0.5, space = 1.5,
              col = c('#3182bd90','#de2d2690'),
              border = FALSE, las = 2, cex.axis = 0.9, ylab = 'Percent')
barplot(t(shown_top[,4:5]), width = 0.5,
        space = c(2.5,rep(1.5, nrow(shown_top))), col = c('#3182bd','#de2d26'),
        border = TRUE, add = TRUE, axes = FALSE, axisnames = FALSE)
legend(x=0.5,y=1,'topright', legend = c('Primary clonal', 'Primary subclonal', 'Metastases clonal', 'Metastases subclonal'),
       fill = c('#3182bd90','#de2d2690','#3182bd','#de2d26'), border = TRUE, bty = 'n', ncol = 2, cex = 0.7)

dev.off()
