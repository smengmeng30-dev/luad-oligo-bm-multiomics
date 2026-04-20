
# 原发和脑转整体的甲基化水平比较 ---------------------------------------------------------
methy_level <-read.delim("all.cpg_human.xls",check.names = F)
clincal <- read.delim("../clinical.group.txt")

mc_Mal <- merge(methy_level,clincal,by.x = 'Sample',by.y='sampleid.gm')

colnames(mc_Mal) <- gsub('\\(.*\\)','',colnames(mc_Mal)) # 去除括号及括号内%
mc_Mal$mCpG <- as.numeric(mc_Mal$mCpG)/100 # 暂时只用了mCpG值
mc_Mal$mCHG <- as.numeric(mc_Mal$mCHG)/100 # 暂时只用了mCpG值
mc_Mal$mCHH <- as.numeric(mc_Mal$mCHH)/100 # 暂时只用了mCpG值
mc_Mal$mC <- as.numeric(mc_Mal$mC)/100 # 暂时只用了mCpG值

mc_Mal$group <- factor(mc_Mal$group,levels = c('Primary','Metastases'))

my_comparisons <- c('Primary','Metastases')


ggboxplot(mc_Mal,x='group',y='mCpG',
          add='jitter',color='group',
          palette = c("#0C7094","#E21D31"),
          xlab = '',
          ylab ="Rates of mCpG / total C")+stat_compare_means()+theme_bw()+theme(panel.grid = element_blank(),legend.position = 'none')

ggsave(file="../paper1.Figure/mCpG_global_meth.pdf",width=3,height=3)


ggboxplot(mc_Mal,x='group',y='mCHG',
          add='jitter',color='group',
          palette = c("#0C7094","#E21D31"),
          ylab ="Rates of mCHG / total C")+stat_compare_means()+theme_bw()+theme(panel.grid = element_blank(),legend.position = 'none')

ggsave(file="../paper1.Figure/mCHG_global_meth.pdf",width=3,height=3)


ggboxplot(mc_Mal,x='group',y='mCHH',
          add='jitter',color='group',
          palette = c("#0C7094","#E21D31"),
          ylab ="Rates of mCHH / total C")+stat_compare_means()+theme_bw()+theme(panel.grid = element_blank(),legend.position = 'none')

ggsave(file="../paper1.Figure/mCHH_global_meth.pdf",width=3,height=3)


ggboxplot(mc_Mal,x='group',y='mC',
          add='jitter',color='group',
          palette = c("#0C7094","#E21D31"),
          ylab ="Rates of mC / total C")+stat_compare_means()+theme_bw()+theme(panel.grid = element_blank(),legend.position = 'none')

ggsave(file="../paper1.Figure/mC_global_meth.pdf",width=3,height=3)



### paired

paired_mal <- dcast(mc_Mal,PID~group,value.var = 'mCpG')

ggpaired(paired_mal,cond1 = 'Primary',cond2 = 'Metastases',
         color='condition',palette = c("#0C7094","#E21D31"),
         ylab='mCpG',
         xlab = "",
         line.color = 'lightgrey')+stat_compare_means(paired = T)+theme_bw()+theme(panel.grid = element_blank(),legend.position = 'none')

ggsave(file="../paper1.Figure/mCpG_global_meth.paired.pdf",width=3,height=3)


paired_mal <- dcast(mc_Mal,PID~group,value.var = 'mCHG')
ggpaired(paired_mal,cond1 = 'Primary',cond2 = 'Metastases',
         color='condition',palette = c("#0C7094","#E21D31"),
         ylab='mCHG',
         xlab = "",
         line.color = 'lightgrey')+stat_compare_means(paired = T)+theme_bw()+theme(panel.grid = element_blank(),legend.position = 'none')

ggsave(file="../paper1.Figure/mCHG_global_meth.paired.pdf",width=3,height=3)

paired_mal <- dcast(mc_Mal,PID~group,value.var = 'mCHH')
ggpaired(paired_mal,cond1 = 'Primary',cond2 = 'Metastases',
         color='condition',palette = c("#0C7094","#E21D31"),
         ylab='mCHH',
         xlab = "",
         line.color = 'lightgrey')+stat_compare_means(paired = T)+theme_bw()+theme(panel.grid = element_blank(),legend.position = 'none')

ggsave(file="../paper1.Figure/mCHH_global_meth.paired.pdf",width=3,height=3)

paired_mal <- dcast(mc_Mal,PID~group,value.var = 'mC')
ggpaired(paired_mal,cond1 = 'Primary',cond2 = 'Metastases',
         color='condition',palette = c("#0C7094","#E21D31"),
         ylab='mC',
         xlab = "",
         line.color = 'lightgrey')+stat_compare_means(paired = T)+theme_bw()+theme(panel.grid = element_blank(),legend.position = 'none')

ggsave(file="../paper1.Figure/mC_global_meth.paired.pdf",width=3,height=3)

# cpg聚类 ---------------------------------------------------------------------

library(data.table)
cpg.all <- as.data.frame(fread("metilene.DMR/DMC/metilene_case_ctrl.input"))
group <- read.delim("metilene.DMR/metastases_vs_primary.txt")
colnames(cpg.all)[3:ncol(cpg.all)] <- c(unlist(strsplit(group$metastases,",")),unlist(strsplit(group$primary,",")))

clinical <- read.delim("../clinical.group.txt")

rownames(cpg.all) <- paste(cpg.all$chrom,cpg.all$pos,sep=':')
cpg.all <- cpg.all[,c(-1,-2)]

cpg.all1 <- cpg.all[rowSums(!is.na(cpg.all)) >= 40*0.9 & rowSums(cpg.all,na.rm = TRUE)!=0,]

dim(cpg.all1)

high_var_cpg <- selectHighVari(cpg.all1,top=T)

p <- plot(high_var_cpg,"../paper1.Figure/Top1000.highVariableCpG.pdf")

saveRDS(high_var_cpg,file = '../paper1.Figure/Top1000.highVariableCpG.rds')

gc()

g1 <- cpg.all[,clinical$sampleid.gm[clinical$group=='Metastases']]
g11 <- g1[rowSums(!is.na(g1))>=ncol(g1)*0.9 & rowSums(g1,na.rm = TRUE)!=0 ,]
g11_high <- selectHighVari(g11,top=T)

plot(g11_high,"../paper1.Figure/Top1000.meta.highVariableCpG.pdf")

gc()

g3 <- cpg.all[,clinical$sampleid.gm[clinical$group=='Primary']]
g31 <- g3[rowSums(!is.na(g3))>=ncol(g3)*0.9 & rowSums(g3,na.rm = TRUE)!=0,]
g31_high <- selectHighVari(g31,top=T)

plot(g31_high,"../paper1.Figure/Top1000.pri.highVariableCpG.pdf")

gc()

top <- cpg.all1[rownames(cpg.all1) %in% unique(c(rownames(g11_high),rownames(g31_high))),]

plot(top,"../paper1.Figure/Top2000.highVariableCpG.pdf")

gc()



##############一致性聚类

selectHighVari <- function(cpg,top=T){
  
  ###### 高变
  variances <- apply(cpg, 1, function(x){var(x,na.rm = TRUE)})  # 计算每个位点的方差
  # 计算前1%的阈值
  threshold <- quantile(variances, 0.99)
  
  # 筛选出变异程度最高的前1%的CpG位点
  high_var_cpg <- cpg[variances > as.numeric(threshold), ]
  print(dim(high_var_cpg))
  
  top.res <- NULL
  if(top){
    top.res <- high_var_cpg[1:1000,]
  }else{
    top.res <- high_var_cpg
  }
  return(top.res)
}

plot <- function(dat,out2){
  
  library(ComplexHeatmap)
  group <- read.delim("../variation.group.txt")
  clinical <- read.delim("../clinical.group.txt")
  group <- merge(clinical,group,by.y='吉因加编码',by.x='sampleid.1021')
  rownames(group) <- group$sampleid.gm
  
  attri <- na.omit(unique(c(group$PID.x,group$group,group$time,group$EGFR.19del.21L858R,group$ALK.Fusion,group$MET.Fusion,group$KRAS,group$TP53,group$RBM10,group$CDKN2A.B.Loss,group$MYC.Gain,group$EGFR.gain,group$RB1.loss,group$ZFHX3)))
  
  
  col <- c(ggsci::pal_npg(palette = 'nrc')(10),ggsci::pal_futurama('planetexpress')(12),"#F3A9AB","#F39E2E","#BCBCC2","#EA500D","#90B5D8","#005599","#7A609F","#E46FA1","#264C53","#5BA98E","#009BC9","#009045","#B34345","red","blue","red")[1:length(attri)]
  
  names(col) <- attri
  
  col.list <- list(PID=col[na.omit(unique(group$PID.x))],Site=col[na.omit(unique(group$group))],Time=col[na.omit(unique(group$time))],EGFR=col[na.omit(unique(group$EGFR.19del.21L858R))],EGFR.gain=col[na.omit(unique(group$EGFR.gain))],ALK=col[na.omit(unique(group$ALK.Fusion))],MET=col[na.omit(unique(group$MET.Fusion))],KRAS=col[na.omit(unique(group$KRAS))],TP53=col[na.omit(unique(group$TP53))],RBM10=col[na.omit(unique(group$RBM10))],`CDKN2A/B`=col[na.omit(unique(group$CDKN2A.B.Loss))],MYC=col[na.omit(unique(group$MYC.Gain))],RB1.loss=col[na.omit(unique(group$RB1.loss))],ZFHX3=col[na.omit(unique(group$ZFHX3))])
  
  
  group <- group[colnames(dat),]
  print(identical(rownames(group),colnames(dat)))
  
  pdf(file = out2,width = 10,height = 6)
  p <- Heatmap(as.matrix(dat),show_column_names  = F,
               show_row_names = F,
               #row_km = 2,
               #clustering_method_rows = 'pam',
               bottom_annotation = columnAnnotation(PID=group$PID.x,Site=group$group,Time=group$time,EGFR=group$EGFR.19del.21L858R,EGFR.gain=group$EGFR.gain,ALK=group$ALK.Fusion,MET=group$MET.Fusion,KRAS=group$KRAS,TP53=group$TP53,RBM10=group$RBM10,`CDKN2A/B`=group$CDKN2A.B.Loss,MYC=group$MYC.Gain,RB1.loss=group$RB1.loss,ZFHX3=group$ZFHX3,
                                                    col=col.list,simple_anno_size = unit(2.5, "mm"),annotation_name_gp = gpar(fontsize = 5)),
               #left_annotation = rowAnnotation(CpG=top.anno1$annot.type,col=list(CpG=col[40:43]),simple_anno_size = unit(2.5, "mm")),
               #right_annotation =rowAnnotation(anno=anno_mark(at=anno_at,labels=anno_mark,labels_gp = gpar(fontsize = 5)))
  )
  print(p)
  dev.off()
  return(p)
}


tf_plot <- function(path,qvalue=0.1,pvalue=0.05){
  hyper <- read.delim(paste(path,"/motif/hyper/knownResults.txt",sep=""),check.names = F)
  
  hyper$`log10(q.value)` <- -log10(hyper$`q-value (Benjamini)`)
  hyper$motifRank <- 1:nrow(hyper)
  hyper$TF <- gsub("\\(.*","",hyper$`Motif Name`)
  hyper$group <- 'hyper'
  
  hypo <- read.delim(paste(path,"/motif/hypo/knownResults.txt",sep=""),check.names = F)
  
  hypo$`log10(q.value)` <- log10(hypo$`q-value (Benjamini)`)
  hypo$motifRank <- 1:nrow(hypo)
  hypo$TF <- gsub("\\(.*","",hypo$`Motif Name`)
  hypo$group <- 'hypo'
  
  
  all <- rbind(hyper[,c("motifRank","TF","group","log10(q.value)",'q-value (Benjamini)','P-value')],hypo[,c("motifRank","TF","group","log10(q.value)",'q-value (Benjamini)','P-value')])
  
  library(ggpubr)
  
  all$sizeColor <- ifelse(all$`q-value (Benjamini)`<qvalue & all$`P-value`<pvalue,'red','grey')
  
  
  show_label <- all[all$`q-value (Benjamini)`<qvalue & all$`P-value`<pvalue,]
  
  p <- ggscatter(all,x='motifRank',y='log10(q.value)',
                 color = 'sizeColor',
                 palette = c('grey','red'),
                 #label = 'TF',
                 #label.select = list(criteria = "`q-value (Benjamini)` <0.1")
  )+ggrepel::geom_text_repel(data = show_label,aes(x=motifRank,y=`log10(q.value)`,label=TF))
  
  return(list(p,show_label))
}


enricher2 <- function(genelist){
  library(msigdbr)
  library(clusterProfiler)
  keggdb = msigdbr(species = "Homo sapiens", category = "C2",subcategory = "KEGG")
  keggdb <- keggdb[,c('gs_name','gene_symbol')]
  colnames(keggdb) <- c('TERM','GENE')
  
  ekegg.up <- clusterProfiler::enricher(genelist,
                                        TERM2GENE=keggdb,pvalueCutoff = 0.05,
                                        pAdjustMethod = "BH",qvalueCutoff = 0.2)
  
  return(ekegg.up)
  
}

# Tp53模式分析 ---------------------------------------------------------------------
### TP53聚成了两类
p@column_names_param$labels[column_order(p)]
## 一簇是unlist(ht[[1]])；另一簇是unlist(ht[[2]])
ht <- column_dend(p)
c1.tp53.wild <- colnames(high_var_cpg)[unlist(ht[[1]])]
c2.tp53.mut <- colnames(high_var_cpg)[unlist(ht[[2]])]

group <- read.delim("../variation.group.txt")
clinical <- read.delim("../clinical.group.txt")
group <- merge(clinical,group,by.y='吉因加编码',by.x='sampleid.1021')
rownames(group) <- group$sampleid.gm

### 第二个cluster中都是TP53突变型
raw.tp53.mut <- group$sampleid.gm[!is.na(group$TP53) & group$TP53=='TP53']
### 去掉第一个cluster中的Tp53突变型，因为第一个cluster认为是TP53野生型
c1.tp53.wild <- setdiff(c1.tp53.wild,raw.tp53.mut)

### 根据cluster1和cluster2构建cpg位点矩阵
dat.final <- high_var_cpg[,c(c2.tp53.mut,c1.tp53.wild)]
dat.final$p.value <- apply(dat.final,1,function(x) {
  res <- wilcox.test(na.omit(as.numeric(x[1:length(c2.tp53.mut)])),na.omit(as.numeric(x[(length(c2.tp53.mut)+1):length(x)])))
  return(res$p.value)})

dat.final$TP53.mut.mean <- rowMeans(dat.final[,c2.tp53.mut],na.rm = T)
dat.final$TP53.wild.mean <- rowMeans(dat.final[,c1.tp53.wild],na.rm = T)
dat.final$mean.meth.diff <- dat.final$TP53.mut.mean-dat.final$TP53.wild.mean

#### dat.final进行cpg位点注释

devtools::install_github('rcavalcante/annotatr',force = T)
BiocManager::install('TxDb.Hsapiens.UCSC.hg19.knownGene')
library(GenomicRanges)
library(annotatr)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(readxl)

mat <- dat.final
mat$ID <- rownames(mat)
rownames(mat) <- NULL
mat <- tidyr::separate(mat,col='ID',into=c("chr","end"),sep=':')
mat$end <- as.numeric(mat$end)
mat$start <- mat$end-1
mat$chr <- paste("chr",mat$chr,sep="")

mat_regions <- GRanges(mat)
annotations <- readRDS("metilene.DMR/DMR/annotations.db.rds")

dm.all_annotated = annotate_regions(
  regions = mat_regions,
  annotations = annotations,
  ignore.strand = TRUE,
  quiet = FALSE)
# A GRanges object is returned
print(dm.all_annotated)
# Coerce to a data.frame
df_dm.all_annotated = data.frame(dm.all_annotated)
saveRDS(dm.all_annotated,file = '../paper1.Figure/TP53.group.cpg.anno.rds')



# DMC中发生最高频甲基化的基因 ---------------------------------------------------------------------
library(ggpubr)
library(data.table)
library(GenomicRanges)


clinical <- read.delim("../clinical.group.txt")
dmc <- readRDS("metilene.DMR/DMC/DMC.all.annotation.rds")
dmc <- as.data.frame(dmc)

### 差异DMC
dmc1 <- dmc[dmc$p.value<0.05 & abs(dmc$mean.meth.diff)>0.1,]
dmc1$flag <- paste(gsub("^chr","",dmc1$seqnames),dmc1$end,sep=':')

#dmc11 <- unique(dmc1$flag[dmc1$annot.type %in% c('hg19_cpg_islands',"hg19_genes_promoters","hg19_genes_1to5kb")])

#dmc2 <- unique(dmc1[dmc1$flag %in% dmc11 & !is.na(dmc1$annot.symbol),c("flag","annot.symbol","mean.meth.diff")])

### 筛选出脑转移高甲基化，且均值不小于0.1， 或者原发高甲基化，且均值不小于0.1的DMC

#dmc11 <- dmc1[(dmc1$mean.meth.diff>0 & dmc1$metastases>0.1) | (dmc1$mean.meth.diff<0 & dmc1$primary>0.1),]

dmc2 <- unique(dmc1[!is.na(dmc1$annot.symbol),c("flag","annot.symbol","mean.meth.diff")])

### 分别找到hyper和hypo发生高频甲基化的基因
gene.dis.hyper <- as.data.frame(table(dmc2$annot.symbol[dmc2$mean.meth.diff>0]))
gene.dis.hyper <- gene.dis.hyper[order(gene.dis.hyper$Freq,decreasing = T),]
gene.dis.hyper$group <- 'hyper'

gene.dis.hypo <- as.data.frame(table(dmc2$annot.symbol[dmc2$mean.meth.diff<0]))
gene.dis.hypo <- gene.dis.hypo[order(gene.dis.hypo$Freq,decreasing = T),]

gene.dis.hypo$group <- 'hypo'


all.gene.dis <- rbind(gene.dis.hyper[1:20,],gene.dis.hypo[1:20,])

all.gene.dis <- all.gene.dis[order(all.gene.dis$Freq,decreasing = T),]

all.gene.dis$Var1 <- factor(all.gene.dis$Var1,levels = all.gene.dis$Var1[!duplicated(all.gene.dis$Var1)])

ggbarplot(all.gene.dis,x='Var1',y='Freq',color='group',fill='group',position = position_dodge(),palette = c("#D32331","#0E6C8C"),xlab = '')+theme_bw()+theme(panel.grid = element_blank(),axis.text.x = element_text(angle = 90,hjust = 1))

ggsave("../paper1.Figure/highFreq.gene.pdf",width = 8,height = 4)



library(msigdbr)
keggdb = msigdbr(species = "Homo sapiens", category = "C2",subcategory = "KEGG")
keggdb <- keggdb[,c('gs_name','gene_symbol')]
colnames(keggdb) <- c('TERM','GENE')
CAMC <- keggdb$GENE[keggdb$TERM=='KEGG_CELL_ADHESION_MOLECULES_CAMS']
NLRI <- keggdb$GENE[keggdb$TERM=='KEGG_NEUROACTIVE_LIGAND_RECEPTOR_INTERACTION']

intersect(CAMC,gene.dis.hyper$Var1[1:30])
intersect(CAMC,gene.dis.hypo$Var1[1:30])

intersect(NLRI,gene.dis.hyper$Var1[1:30])
intersect(NLRI,gene.dis.hypo$Var1[1:30])

ekegg.up <- clusterProfiler::enricher(gene.dis.hyper$Var1[1:50],
                                      TERM2GENE=keggdb,pvalueCutoff = 1,
                                      pAdjustMethod = "BH",qvalueCutoff = 1)

kegg_res <- ekegg.up@result
kegg_res$generatio <- apply(kegg_res,1,generatio)

p <- enrichplot(kegg_res,"Pathways enriched in metastases (hyper in metastases)")
print(p)

ggsave("../paper1.Figure/highFreqGene.meta.pathway.pdf",width = 6,height = 4)

ekegg.down <- clusterProfiler::enricher(gene.dis.hypo$Var1[1:50],
                                        TERM2GENE=keggdb,pvalueCutoff = 0.05,
                                        pAdjustMethod = "BH",qvalueCutoff = 0.2)

kegg_res <- ekegg.down@result

kegg_res$generatio <- apply(kegg_res,1,generatio)

p <- enrichplot(kegg_res,"Pathways enriched in primary (hyper in primary)")

print(p)

ggsave("../paper1.Figure/highFreqGene.pri.pathway.pdf",width = 8,height = 4)


#### 再从这些高频甲基化基因中查看经常发生甲基化的基因区域
high <- unique(dmc1[!is.na(dmc1$annot.symbol) & dmc1$annot.symbol %in% unique(all.gene.dis$Var1),c("flag","annot.symbol","mean.meth.diff","annot.type")])

high$group <- ifelse(high$mean.meth.diff>0,"hyper","hypo")

high1 <- as.data.frame(table(high[,c("group","annot.symbol","annot.type")]))



# 提取NLGN1的DMC -------------------------------------------------------------

NLGN1 <- dmc1[!is.na(dmc1$annot.symbol) & dmc1$annot.symbol=='NLGN1',]
NLGN1$flag <-  paste(NLGN1$seqnames,NLGN1$start,NLGN1$end,sep=":")

all.cpg <- fread("all.CpG.all.txt",check.names = F)

motif.hyper.input <- unique(NLGN1[NLGN1$mean.meth.diff > 0.1,c("flag","seqnames","start","end")])
motif.hyper.input.strand <- as.data.frame(all.cpg[all.cpg$flag %in% motif.hyper.input$flag,])
rownames(motif.hyper.input.strand) <- motif.hyper.input.strand$flag
motif.hyper.input$strand <- motif.hyper.input.strand[motif.hyper.input$flag,'strand']

motif.hypo.input <- unique(NLGN1[NLGN1$mean.meth.diff < -0.1,c("flag","seqnames","start","end")])
motif.hypo.input.strand <- as.data.frame(all.cpg[all.cpg$flag %in% motif.hypo.input$flag,])
rownames(motif.hypo.input.strand) <- motif.hypo.input.strand$flag
motif.hypo.input$strand <- motif.hypo.input.strand[motif.hypo.input$flag,'strand']


write.table(motif.hyper.input,file = "../paper1.Figure/NLGN1.motif.hyper.input.txt",quote = F,sep='\t',row.names = F)

write.table(motif.hypo.input,file = "../paper1.Figure/NLGN1.motif.hypo.input.txt",quote = F,sep='\t',row.names = F)

### 预测得到的NPAS4转录因子与NLGN1的关系

hyper.tf <- read.delim("../paper1.Figure/NLGN1/hyper/knownResults.txt")
hyper.tf$tf <- toupper(gsub("\\(.*","",hyper.tf$Motif.Name))

hypo.tf <- read.delim("../paper1.Figure/NLGN1/hypo/knownResults.txt")
hypo.tf$tf <- toupper(gsub("\\(.*","",hypo.tf$Motif.Name))

TPM <- read.delim("../RNA-seq/All.TPM_geneSymbol.txt",check.names = F)
group <- read.delim("../clinical.group.txt")
meta.id <- group$sampleid.rna[group$group=='Metastases']
pri.id <- group$sampleid.rna[group$group=='Primary']

for(i in unique(c(hyper.tf$tf[1:20],hypo.tf$tf[1:20]))){
  if(!i %in% TPM$GeneSymbol){next}
  TPM1 <- TPM[TPM$GeneSymbol %in% c(i,"NLGN1"),]
  TPM.meta <- TPM1[,c("GeneSymbol",meta.id)]
  TPM.pri <- TPM1[,c("GeneSymbol",pri.id)]
  
  t1 <- cor.test(as.numeric(TPM.meta[1,2:ncol(TPM.meta)]),as.numeric(TPM.meta[2,2:ncol(TPM.meta)]))
  t2 <- cor.test(as.numeric(TPM.pri[1,2:ncol(TPM.pri)]),as.numeric(TPM.pri[2,2:ncol(TPM.pri)]))
  t3 <- cor.test(as.numeric(TPM1[1,2:ncol(TPM1)]),as.numeric(TPM1[2,2:ncol(TPM1)]))
  
  if(!is.na(t1$p.value) & t1$p.value<0.05){print(i);print(t1$estimate)}
}

library(RTN)
library(snow)
library(ComplexHeatmap)
library(ClassDiscovery)
library(RColorBrewer)
library(gplots)

standarize.fun <- function(indata=NULL, halfwidth=NULL, centerFlag=T, scaleFlag=T) {  
  outdata=t(scale(t(indata), center=centerFlag, scale=scaleFlag))
  if (!is.null(halfwidth)) {
    outdata[outdata>halfwidth]=halfwidth
    outdata[outdata<(-halfwidth)]= -halfwidth
  }
  return(outdata)
}


# 加载基因表达以及样本数值信息
tpm <- read.table("../RNA-seq/geneSymbol_tpm.xls",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
pheno <- read.table("../clinical.group.txt",sep = "\t", check.names = F,stringsAsFactors = F,header = T)

# 加载MIBC特异性的调控子
tfs <- motif
tfs$regulon <- toupper(tfs$TF)
# 取共有的基因名
regulatoryElements <- intersect(tfs$regulon, rownames(tpm))

# 运行TNI构建程序
# we used the R package “RTN” to reconstruct transcriptional regulatory networks (regulons)
tpm.log <- tni.constructor(expData = as.matrix(log2(tpm + 1)), # 样图计算时候没有取对数, 
                           regulatoryElements = regulatoryElements)

# 通过置换以及bootstrap计算reference regulatory network.
# mutual information analysis and Spearman rank-order correlation deduced the possible associations between a regulator and all potential target from the transcriptome expression profile, and permutation analysis was utilized to erase associations with an FDR > 0.00001. Bootstrapping strategy removed unstable associations through one thousand times of resampling with consensus bootstrap greater than 95%. 
# 这里量力而行设置多核，或者直接单核运算
options(cluster=snow::makeCluster(spec = 4, "SOCK")) # 打开4核并行计算（不确定是不是4核，不过我windows只用4，服务器我开12）
tpm.log <- tni.permutation(tpm.log, pValueCutoff = 0.05, nPermutations = 100)
tpm.log <- tni.bootstrap(tpm.log, nBootstraps = 100)
stopCluster(getOption("cluster")) # 关闭并行计算

# 计算DPI-filtered regulatory network
# Data processing inequality filtering eliminated the weakest associations in triangles of two regulators and common targets
tpm.log1 <- tni.dpi.filter(tpm.log, eps = 0, sizeThreshold = TRUE, minRegulonSize = 5)

tni.regulon.summary(tpm.log1)

# 保存TNI对象以便后续分析
save(tpm.log1, file="../paper1.Figure/tpm.log.RData")

# load("rtni_tcgaBLCA.RData")
# 计算每个样本的regulon活性
# Individual regulon activity was estimated by two-sided GSEA
tpm.gsea2 <- tni.gsea2(tpm.log1, regulatoryElements = regulatoryElements)
regact <- tni.get(tpm.gsea2, what = "regulonActivity")

# 保存活性对象
save(regact,file = "../paper1.Figure/regact.RData") 



# DMC motif --------------------------------------------------------------

p1 <- tf_plot("metilene.DMR/DMC/",0.05,0.05)

p1[[1]]+theme(legend.position = 'none')

ggsave("../paper1.Figure/meta_vs_pri.dmc.motif.pdf",width = 5,height = 5)

motif <- p1[[2]]

toprank <- motif[!duplicated(motif$TF),]
hyper <- enricher2(toupper(toprank$TF[toprank$group=='hyper']))
hyper <- hyper@result

hypo <- enricher2(toupper(toprank$TF[toprank$group=='hypo']))
hypo <- hypo@result

common <- unique(motif[,c("TF","group")])
motif1 <- motif[motif$TF %in% names(which(table(common$TF)==1)),]


write.table(motif1,file='../paper1.Figure/sig.TF.txt',quote = F,sep='\t',row.names = F)


# DMC甲基化注释 ----------------------------------------------------------------
devtools::install_github('rcavalcante/annotatr',force = T)
BiocManager::install('TxDb.Hsapiens.UCSC.hg19.knownGene')
library(annotatr)
library(GenomicRanges)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(readxl)

dmc.anno <- readRDS("metilene.DMR/DMC/DMC.all.annotation.rds")
dmc.anno <- dmc.anno[abs(dmc.anno$mean.meth.diff)>0.1 & dmc.anno$p.value<0.05,]
dmc.anno$regulated <- ifelse(dmc.anno$mean.meth.diff>0,'hyper','hypo')

# See the GRanges column of dm_annotaed expanded
print(head(dmc.anno))


annots_order = c(
  'hg19_genes_1to5kb',
  'hg19_genes_promoters',
  'hg19_genes_5UTRs',
  'hg19_genes_exons',
  'hg19_genes_intronexonboundaries',
  'hg19_genes_introns',
  'hg19_genes_3UTRs',
  'hg19_genes_intergenic')

annotated_regions = as.data.frame(dmc.anno, row.names = NULL)
annotated_regions = subset_order_tbl(tbl = annotated_regions, 
                                     col = "annot.type", col_order = annotation_order)
annotated_regions = dplyr::distinct(dplyr::ungroup(annotated_regions), 
                                    across(c("seqnames", "start", "end", "annot.type")), 
                                    .keep_all = TRUE)

ggplot(annotated_regions, aes_string(x = "annot.type")) + 
  geom_bar(aes_string(fill = "annot.type"), position = "dodge") + 
  theme_bw() +scale_fill_manual(values = c("#E55C27","#F1A646","#21579D","#9BBADD","#F2B1B0","#7C65A4","#E37AA9",'black')) + theme(panel.grid = element_blank(),axis.text.x = element_text(angle = 30, 
                                                                    hjust = 1), legend.title = element_blank(), legend.position = "none", 
                                         legend.key = element_rect(color = "white"))+xlab("")+ylab("Count")

ggsave("../paper1.Figure/meta_vs_pri.DMC.anno.bar.pdf",width = 5,height = 4)



# The orders for the x-axis labels.
x_order = c(
  'hg19_genes_1to5kb',
  'hg19_genes_promoters',
  'hg19_genes_5UTRs',
  'hg19_genes_exons',
  'hg19_genes_introns',
  'hg19_genes_3UTRs',
  'hg19_genes_intergenic')
# The orders for the fill labels.
fill_order = c(
  'hyper',
  'hypo',
  'none')
dm_vs_kg_cat = plot_categorical(
  annotated_regions = dmc.anno, x='annot.type', fill='regulated',
  x_order = x_order, fill_order = fill_order, position='fill',
  legend_title = 'DM Status',
  x_label = 'knownGene Annotations',
  y_label = 'Proportion')+theme_bw()

dm_vs_kg_cat+scale_fill_manual(values = c("#E21D31","#0C7094"))

ggsave("../paper1.Figure/meta_vs_pri.DMC.anno.pdf",width = 5,height = 4)


dmc.anno.df <- as.data.frame(dmc.anno)
num <- as.data.frame(table(dmc.anno.df$regulated))

num$percentage <- round(num$Freq/sum(num$Freq),2)*100

labs <- paste0(num$Var1, " (", num$percentage, "%)")
ggpie(num,x='percentage', label = labs,fill = "Var1", color = "white",
      palette = c("#E21D31","#0C7094"))+theme(legend.position = 'none')

ggsave("../paper1.Figure/DMC.hyper.hypo.distribution.pdf",width = 3,height = 3)

# DMC甲基化&表达共分析 -------------------------------------------------------
df1 <- as.data.frame(readRDS("metilene.DMR/DMC/DMC.all.annotation.rds"))

df11 <- df1[abs(df1$mean.meth.diff)>0.1 & df1$p.value<0.05,]
#df11 <- df1[abs(df1$mean.meth.diff)>0.1 & df1$p.value<0.05 & df1$q.value<0.2,]

meta_vs_pri <- read.delim("../RNA-seq/metastases_versus_primary.deseq.xls")
colnames(meta_vs_pri)[1] <- 'GeneID'
meta_vs_pri <- meta_vs_pri[abs(meta_vs_pri$log2FoldChange)>=1 & meta_vs_pri$pvalue<0.05,]

p10 <- co_analysis(df11,meta_vs_pri,show_label = F)
p10[[1]]

ggsave("../paper1.Figure/co_expr.pdf",width=6,height=6)

pro_gb <- p10[[2]]

write.table(p10[[2]],file = '../paper1.Figure/deg_dmc.coanalysis.txt',quote = F,sep='\t',row.names = F)

up <- pro_gb$annot.symbol[(pro_gb$group=='gene body' & pro_gb$mean.meth.diff>0 & pro_gb$log2FoldChange>0)|
                            (pro_gb$group=='promoter' & pro_gb$mean.meth.diff<0 & pro_gb$log2FoldChange>0)]

down <- pro_gb$annot.symbol[(pro_gb$group=='gene body' & pro_gb$mean.meth.diff<0 & pro_gb$log2FoldChange<0)|
                              (pro_gb$group=='promoter' & pro_gb$mean.meth.diff>0 & pro_gb$log2FoldChange<0)]


library(msigdbr)
keggdb = msigdbr(species = "Homo sapiens", category = "C2",subcategory = "KEGG")
keggdb <- keggdb[,c('gs_name','gene_symbol')]
colnames(keggdb) <- c('TERM','GENE')
NLRI <- keggdb$GENE[keggdb$TERM=='KEGG_NEUROACTIVE_LIGAND_RECEPTOR_INTERACTION']

intersect(NLRI,up)
intersect(NLRI,down)

ekegg.up <- clusterProfiler::enricher(unique(down),
                                      TERM2GENE=keggdb,pvalueCutoff = 1,
                                      pAdjustMethod = "BH",qvalueCutoff = 1)

kegg_res <- ekegg.up@result[ekegg.up@result$p.adjust < 0.05,]

kegg_res$generatio <- apply(kegg_res,1,generatio)

p <- enrichplot(kegg_res,"Pathways enriched in Primary (genes upregulated in primary but hyper in metastases)")
print(p)

ggsave("../paper1.Figure/co.meta.pathway.pdf",width = 6,height = 4)

ekegg.down <- clusterProfiler::enricher(unique(up),
                                        TERM2GENE=keggdb,pvalueCutoff = 0.05,
                                        pAdjustMethod = "BH",qvalueCutoff = 0.2)

kegg_res <- ekegg.down@result[ekegg.down@result$p.adjust < 0.05,]

kegg_res$generatio <- apply(kegg_res,1,generatio)

p <- enrichplot(kegg_res,"Pathways enriched in metastases (genes upregulated in metastases but hyper in primary)")

print(p)

ggsave("../paper1.Figure/co.primary.pathway.pdf",width = 8,height = 4)



co_analysis <- function(df_dm_annotated,deg,show_label=F){
  
  exp <- read.delim("../RNA-seq/all_samples_cluster_TPM.xls",check.names = F)
  
  deg <- merge(deg,exp[,1:2],by='GeneID')
  
  combined <- merge(df_dm_annotated,deg,by.x='annot.symbol',by.y='Symbol')
  
  combined$group[combined$annot.type=="hg19_genes_promoters"] <- 'promoter'
  
  combined$group[combined$annot.type %in% c('hg19_genes_introns','hg19_genes_exons','hg19_genes_intronexonboundaries','hg19_genes_5UTRs')] <- 'gene body'
  
  combined <- combined[!is.na(combined$group),]
  
  combined1 <- unique(combined[,c("annot.symbol","seqnames","start","end","mean.meth.diff","log2FoldChange",'group')])
  
  combined1$flag <- paste(combined1$seqnames,combined1$start,combined1$end,combined1$annot.symbol,sep = ":")
  
  combined1 <- combined1[(combined1$mean.meth.diff*combined1$log2FoldChange>0 & combined1$group=='gene body') | (combined1$mean.meth.diff*combined1$log2FoldChange<0 & combined1$group=='promoter'),]
  
  combined1$label <- ifelse((abs(combined1$log2FoldChange)>4) | abs(combined1$mean.meth.diff)>0.2,combined1$annot.symbol,"")
  
  p <- ggscatter(combined1,x='mean.meth.diff',y='log2FoldChange',fill = 'group',palette = 'npg',shape  = 21)+
    geom_vline(xintercept=c(0.1,-0.1),linetype=2)+
    geom_hline(yintercept=c(-1,1),linetype=2)
  
  if(show_label){
    
    p <- p +   ggrepel::geom_text_repel(aes(label=label),max.overlaps = 100)
    
  }
  
  return(list(p,combined1))
}


generatio <- function(x){
  ratio <- strsplit(x[3],"/")
  count <- as.integer(ratio[[1]][1])
  total <- as.integer(ratio[[1]][2])
  return(count/total)
}


enrichplot <- function(data,title){
  library(ggplot2)
  library(dplyr)
  #library(ggthemes)
  
  #data <- data %>% group_by(cell) %>% top_n(n=5,wt=-p.adjust)
  
  p <- ggplot(data)+geom_point(aes(y=generatio,x=reorder(Description,generatio),fill=-log10(p.adjust),size=Count),shape=21,colour='black')+
    
    ggtitle(title)+
    
    coord_flip()+
    
    scale_fill_gradient(low = 'white',high = '#af2934')+
    
    theme_light()+
    
    theme(axis.text.x = element_text(angle = 90,hjust=1,vjust=0.5),axis.text=element_text(size = 10,color = "black"))+
    
    theme(axis.text.y = element_text(size = 10),axis.title = element_text(size = 10))+
    
    #facet_grid(rows=vars(Class),scales = "free_y",space = "free_y")+
    
    theme(strip.background=element_rect(fill = c("blue")))
  return(p)
}



# 比较NLGN1的甲基化和表达的相关性 ------------------------------------------------------
NLGN1.dmc <- df1[df1$annot.symbol=='NLGN1',]
promoter <- NLGN1.dmc[NLGN1.dmc$annot.type=='hg19_genes_promoters',]
genbody <- NLGN1.dmc[NLGN1.dmc$annot.type %in% c('hg19_genes_introns','hg19_genes_exons','hg19_genes_intronexonboundaries','hg19_genes_5UTRs'),]
promoter.dmc <- as.data.frame(all.cpg[all.cpg$flag %in% unique(paste(promoter$seqnames,promoter$start,promoter$end,sep=":")),])
genbody.dmc <- as.data.frame(all.cpg[all.cpg$flag %in% unique(paste(genbody$seqnames,genbody$start,genbody$end,sep=":")),])

rownames(genbody.dmc) <- genbody.dmc$flag
genbody.dmc <- as.data.frame(t(genbody.dmc[,c(-1,-42,-43)]))

genbody.dmc$mean <- rowMeans(genbody.dmc,na.rm = T)

group <- read.delim("../clinical.group.txt")
genbody.dmc <- merge(genbody.dmc,group,by.x='row.names',by.y='sampleid.gm')
tpm <- read.delim("../RNA-seq/All.TPM_geneSymbol.xls",check.names = F,row.names = 1)
nlgn1.tpm <- as.data.frame(t(tpm['NLGN1',,drop=F]))
genbody.dmc <- merge(genbody.dmc,nlgn1.tpm,by.x='sampleid.rna',by.y='row.names')

cor.test(genbody.dmc$mean,genbody.dmc$NLGN1)
cor.test(genbody.dmc$mean[genbody.dmc$group=='Primary'],genbody.dmc$NLGN1[genbody.dmc$group=='Primary'])
cor.test(genbody.dmc$mean[genbody.dmc$group=='Metastases'],genbody.dmc$NLGN1[genbody.dmc$group=='Metastases'])




rownames(promoter.dmc) <- promoter.dmc$flag
promoter.dmc <- as.data.frame(t(promoter.dmc[,c(-1,-42,-43)]))

promoter.dmc$mean <- rowMeans(promoter.dmc,na.rm = T)

group <- read.delim("../clinical.group.txt")
promoter.dmc <- merge(promoter.dmc,group,by.x='row.names',by.y='sampleid.gm')
tpm <- read.delim("../RNA-seq/All.TPM_geneSymbol.xls",check.names = F,row.names = 1)
nlgn1.tpm <- as.data.frame(t(tpm['NLGN1',,drop=F]))
promoter.dmc <- merge(promoter.dmc,nlgn1.tpm,by.x='sampleid.rna',by.y='row.names')

cor.test(promoter.dmc$mean,promoter.dmc$NLGN1)
cor.test(promoter.dmc$mean[promoter.dmc$group=='Primary'],promoter.dmc$NLGN1[promoter.dmc$group=='Primary'])
cor.test(promoter.dmc$mean[promoter.dmc$group=='Metastases'],promoter.dmc$NLGN1[promoter.dmc$group=='Metastases'])

# DMC cpg types basis -----------------------------------------------------
library(data.table)


dmc.tp53 <- cpg.bais('../paper1.Figure/TP53.group.cpg.anno.rds',NULL,"TP53+ .vs. TP53-")

dmc.res <- cpg.bais("metilene.DMR/DMC/DMC.all.annotation.rds",NULL,"meta_vs_pri")

egfr.positive.res <- cpg.bais("EGFR.DMC/EGFR_positive/DMC.all.annotation.rds",NULL,"EGFR.positive.pri_vs_meta")

meta.posVSneg.res <- cpg.bais("EGFR.DMC/EGFR_positive_meta/DMC.all.annotation.rds",NULL,"meta.EGFR.positive_vs_negative")

pri.posVSneg.res <- cpg.bais("EGFR.DMC/EGFR_positive_pri/DMC.all.annotation.rds",NULL,"pri.EGFR.positive_vs_negative")


all <- rbind(dmc.tp53,dmc.res,egfr.positive.res,meta.posVSneg.res,pri.posVSneg.res)

write.table(all,file='../paper1.Figure/all.fisher.txt',quote=F,sep='\t',row.names=F)

all1 <- dcast(all,cpg_type~group+tag,value.var = 'OR')
rownames(all1) <- all1$cpg_type
all1 <- all1[,-1]


all2 <- dcast(all,cpg_type~group+tag,value.var = 'pvalue')
rownames(all2) <- all2$cpg_type
all2 <- all2[,-1]

library(ComplexHeatmap)
library(circlize)

all11 <- t(log10(all1))
all12 <- t(all1)
all21 <- t(all2)

pdf(file = '../paper1.Figure/cpg.basis.pdf',width = 6,height = 6)
Heatmap(all11,
        border = F,
        name = 'log10(OR)',
        row_split = 1:nrow(all11),
        col=colorRamp2(c(min(log10(all$OR)),0,max(setdiff(log10(all$OR),Inf))),c("#3980AD","white","#D97166")),
        rect_gp = gpar(col = 'black'),
        cluster_rows = F,cluster_columns = F,
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(sprintf("%.1f", all12[i, j]), x, y, gp = gpar(fontsize = 10))
          if(!is.na(all21[i,j]) & all21[i,j]<0.05){
            grid.rect(x = x, y = y, width = width, height = height, 
                      gp = gpar(col = "black", fill = NA))
          }
        }
        
)
dev.off()


cpg.bais <- function(all.input,dmc.input=NULL,tag=NULL){
  
  all <- as.data.frame(readRDS(all.input))
  all.count <- unique(all[all$annot.type %in% c("hg19_cpg_inter","hg19_cpg_shores","hg19_cpg_shelves","hg19_cpg_islands"),c("seqnames","start","end",'annot.type')])
  
  all.count <- as.data.frame(table(all.count$annot.type))
  
  if(is.null(dmc.input)){
    dmr <- all[all$p.value<0.05 & abs(all$mean.meth.diff)>0.1,]
  }else{
    dmr <- as.data.frame(readRDS(dmc.input))
  }
  
  hypo <- unique(dmr[dmr$mean.meth.diff<0 & dmr$annot.type %in% c("hg19_cpg_inter","hg19_cpg_shores","hg19_cpg_shelves","hg19_cpg_islands"),c("seqnames","start","end",'annot.type')])
  
  hypo.count <- as.data.frame(table(hypo$annot.type))
  
  ### hyper
  hyper <- unique(dmr[dmr$mean.meth.diff>0 & dmr$annot.type %in% c("hg19_cpg_inter","hg19_cpg_shores","hg19_cpg_shelves","hg19_cpg_islands"),c("seqnames","start","end",'annot.type')])
  hyper.count <- as.data.frame(table(hyper$annot.type))
  
  ### 相比全部cpg位点，hypo的cpg位点富集到island，还是shlef，shore？
  all.res <- NULL
  for(i in unique(c(hypo.count$Var1,hyper.count$Var1))){
    
    res1 <- NULL
    res2 <- NULL
    
    if(i %in% hypo.count$Var1){
      ### 分别为该cpg类型的hypo DMC数目、非该cpg类型的hypo DMC数目、该cpg类型的所有cpg数目、非该cpg类型的所有cpg数目
    res1 <- fisher.test(matrix(c(hypo.count$Freq[hypo.count$Var1==i],
                                 sum(hypo.count$Freq[hypo.count$Var1!=i]),
                                 all.count$Freq[all.count$Var1==i],
                                 sum(all.count$Freq[all.count$Var1!=i])),byrow = T,nrow = 2))
    
    }
    
    if(i %in% hyper.count$Var1){
      res2 <- fisher.test(matrix(c(hyper.count$Freq[hyper.count$Var1==i],
                                   sum(hyper.count$Freq[hyper.count$Var1!=i]),
                                   all.count$Freq[all.count$Var1==i],
                                   sum(all.count$Freq[all.count$Var1!=i])),byrow = T,nrow = 2))
    }
    
    if(!is.null(res1)){
      tmp.res <- data.frame(group=c('hypo'),
                            OR=c(res1$estimate),
                            pvalue=c(res1$p.value),
                            cpg_type=c(i),
                            dmc.cpg=c(hypo.count$Freq[hypo.count$Var1==i]),
                            dmc.cpg.non=c(sum(hypo.count$Freq[hypo.count$Var1!=i])),
                            all.cpg=c(all.count$Freq[all.count$Var1==i]),
                            all.cpg.non=c(sum(all.count$Freq[all.count$Var1!=i])))
      all.res <- rbind(all.res,tmp.res)
    }
    if(!is.null(res2)){
      tmp.res <- data.frame(group=c('hyper'),
                            OR=c(res2$estimate),
                            pvalue=c(res2$p.value),
                            cpg_type=c(i),
                            dmc.cpg=c(hyper.count$Freq[hyper.count$Var1==i]),
                            dmc.cpg.non=c(sum(hyper.count$Freq[hyper.count$Var1!=i])),
                            all.cpg=c(all.count$Freq[all.count$Var1==i]),
                            all.cpg.non=c(sum(all.count$Freq[all.count$Var1!=i])))
      
      all.res <- rbind(all.res,tmp.res)
    }
  }
  all.res$tag <- tag
  return(all.res)
}


# EGFR variation vs expr vs methylation -----------------------------------


# TF活性-----------------------------------------------------------

BiocManager::install("RTN")
BiocManager::install("ClassDiscovery")

library(RTN)
library(snow)
library(ComplexHeatmap)
library(ClassDiscovery)
library(RColorBrewer)
library(gplots)

standarize.fun <- function(indata=NULL, halfwidth=NULL, centerFlag=T, scaleFlag=T) {  
  outdata=t(scale(t(indata), center=centerFlag, scale=scaleFlag))
  if (!is.null(halfwidth)) {
    outdata[outdata>halfwidth]=halfwidth
    outdata[outdata<(-halfwidth)]= -halfwidth
  }
  return(outdata)
}


# 加载基因表达以及样本数值信息
tpm <- read.table("../RNA-seq/geneSymbol_tpm.xls",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
pheno <- read.table("../clinical.group.txt",sep = "\t", check.names = F,stringsAsFactors = F,header = T)

# 加载MIBC特异性的调控子
tfs <- motif
tfs$regulon <- toupper(tfs$TF)
# 取共有的基因名
regulatoryElements <- intersect(tfs$regulon, rownames(tpm))

# 运行TNI构建程序
# we used the R package “RTN” to reconstruct transcriptional regulatory networks (regulons)
tpm.log <- tni.constructor(expData = as.matrix(log2(tpm + 1)), # 样图计算时候没有取对数, 
                           regulatoryElements = regulatoryElements)

# 通过置换以及bootstrap计算reference regulatory network.
# mutual information analysis and Spearman rank-order correlation deduced the possible associations between a regulator and all potential target from the transcriptome expression profile, and permutation analysis was utilized to erase associations with an FDR > 0.00001. Bootstrapping strategy removed unstable associations through one thousand times of resampling with consensus bootstrap greater than 95%. 
# 这里量力而行设置多核，或者直接单核运算
options(cluster=snow::makeCluster(spec = 4, "SOCK")) # 打开4核并行计算（不确定是不是4核，不过我windows只用4，服务器我开12）
tpm.log <- tni.permutation(tpm.log, pValueCutoff = 0.05, nPermutations = 100)
tpm.log <- tni.bootstrap(tpm.log, nBootstraps = 100)
stopCluster(getOption("cluster")) # 关闭并行计算

# 计算DPI-filtered regulatory network
# Data processing inequality filtering eliminated the weakest associations in triangles of two regulators and common targets
tpm.log1 <- tni.dpi.filter(tpm.log, eps = 0, sizeThreshold = TRUE, minRegulonSize = 5)

tni.regulon.summary(tpm.log1)

# 保存TNI对象以便后续分析
save(tpm.log1, file="../paper1.Figure/tpm.log.RData")

# load("rtni_tcgaBLCA.RData")
# 计算每个样本的regulon活性
# Individual regulon activity was estimated by two-sided GSEA
tpm.gsea2 <- tni.gsea2(tpm.log1, regulatoryElements = regulatoryElements)
regact <- tni.get(tpm.gsea2, what = "regulonActivity")

# 保存活性对象
save(regact,file = "../paper1.Figure/regact.RData") 




# 设置颜色
clust.col <- c("#0C7094","#E21D31")
blue <- "#5bc0eb"
gold <- "#ECE700"

annCol <- pheno[order(pheno$group),] # 构建样本注释信息，并对亚型进行排序
rownames(annCol) <- annCol$sampleid.rna

regulon <- regact$differential[rownames(annCol),]

plotdata <- standarize.fun(t(regulon),halfwidth = 1.5) # 标准化regulon的活性

annColors <- list()
annColors[["group"]] <- c("Primary" = clust.col[1],
                          "Metastases" = clust.col[2]
)
hcg <- hclust(distanceMatrix(as.matrix(regulon), "euclidean"), "ward.D")




hm <- pheatmap(plotdata[hcg$order,],
               border_color = NA, # 热图单元格无边框
               #color = colorpanel(64,low=blue,mid = "black",high=gold),
               cluster_rows = T, # 行不聚类
               cluster_cols = F, # 列聚类
               show_rownames = T, # 显示行名
               show_colnames = F, # 不显示列名
               gaps_col = cumsum(table(annCol$group))[1:2], # 亚型分割
               #cellwidth = 0.8, # 固定单元格宽度
               #cellheight = 10, # 固定单元格高度
               name = "Regulon", # 图例名字
               annotation_col = annCol[,"group",drop = F], # 样本注释
               annotation_colors = annColors["group"]) # 样本注释的对应颜色

pdf("../paper1.Figure/regulon.heatmap.pdf", width = 8,height = 6)
draw(hm) # 输出热图
invisible(dev.off())


############# 加上样本属性，突变、原发转移类型等

group <- read.delim("../clinical.group.txt")
var <- read.delim("../variation.group.txt")
group <- merge(group,var,by.x='sampleid.1021',by.y='吉因加编码')
rownames(group) <- group$sampleid.rna

group <- group[rownames(regact$differential),]


attri <- unique(c(group$group,group$time,group$therapy,group$EGFR.19del.21L858R,group$ALK.Fusion,group$MET.Fusion,group$KRAS,group$TP53,group$RBM10,group$CDKN2A.B.Loss,group$MYC.Gain,group$EGFR.gain,group$RB1.loss,group$ZFHX3))

col <- c(ggsci::pal_npg(palette = 'nrc')(10),ggsci::pal_futurama('planetexpress')(12),"#F3A9AB","#F39E2E","#BCBCC2","white","#EA500D","#90B5D8","#005599","#7A609F","#E46FA1","#264C53","#5BA98E","#009BC9","#009045","#B34345","red","blue","red")[1:length(attri)]

names(col) <- attri

col=list(group=col[1:2],time=col[3:5],therapy=col[6:9],`EGFR.19del.21L858R`=col[10:11],ALK.Fusion=col[12],MET.Fusion=col[13],KRAS=col[14:15],TP53=col[16],RBM10=col[17],`CDKN2A.B.Loss`=col[18],MYC.Gain=col[19],EGFR.gain=col[20],RB1.loss=col[21],ZFHX3=col[22])

pdf("../paper1.Figure/regulon.heatmap2.pdf", width = 10,height = 8)
pheatmap(t(regact$differential),
         cluster_rows = T,cluster_cols = T,
         show_colnames = F,
         fontsize_row = 6,
         clustering_method = "ward.D2", 
         clustering_distance_rows = "correlation",
         clustering_distance_cols = "correlation",
         border_color = NA,
         annotation_col = group[,c("group","time","therapy","EGFR.19del.21L858R","ALK.Fusion" ,"MET.Fusion","KRAS","TP53","RBM10","CDKN2A.B.Loss","MYC.Gain","EGFR.gain","RB1.loss","ZFHX3"),drop = F],
         annotation_colors = col
)

dev.off()


############# 发现原发和转移呈现不同的TF活性模式，再进一步比较在两组中显著的TF活性
dat <- as.data.frame(regact$differential[rownames(annCol),])
dat$group <- annCol$group
dat$ID <- rownames(dat)
dat <- melt(dat,measure.var=1:31)

##两组间显著的TF活性

final.tfs <- c("ATF1","FOXF1","GATA6","MEF2A","SOX10","SOX21","SP1","TBX5")

ggboxplot(dat[dat$variable %in% final.tfs,],x='variable',y='value',xlab = '',ylab='Regulon score',
          color='group',palette = c("#0C7094","#E21D31"),add = 'jitter')+stat_compare_means(aes(group=group,label = ..p.format..))

ggsave("../paper1.Figure/TF.regulonScore.pdf",width = 7,height = 4)



load("../paper1.Figure/tpm.log.RData")

#权重的绝对值表示 MI 值，而符号 （+/-） 表示基于调节器与其目标之间的 Pearson 相关性的预测作用方式。

regulons <- tni.get(tpm.log1, what = "regulons.and.mode", idkey = "ID")

intersect(NLRI,names(regulons$ATF1))
intersect(CAMC,names(regulons$ATF1))
intersect("NLGN1",names(regulons$ATF1))

head(regulons$ATF1)

intersect(NLRI,names(regulons$FOXF1))
intersect(CAMC,names(regulons$FOXF1))
intersect("NLGN1",names(regulons$ATF1))

head(regulons$FOXF1)

intersect(NLRI,names(regulons$GATA6))
intersect(CAMC,names(regulons$GATA6))
intersect("NLGN1",names(regulons$GATA6))


head(regulons$GATA6)
intersect(NLRI,names(regulons$MEF2A))
intersect(CAMC,names(regulons$MEF2A))
intersect("NLGN1",names(regulons$MEF2A))

head(regulons$MEF2A)
intersect(NLRI,names(regulons$TBX5))
intersect(CAMC,names(regulons$TBX5))
intersect("NLGN1",names(regulons$TBX5))


head(regulons$TBX5)

## SOX10的靶基因TBXA2R与DFS有关，
intersect(NLRI,names(regulons$SOX10))
intersect(CAMC,names(regulons$SOX10))
intersect("NLGN1",names(regulons$SOX10))

head(regulons$SOX10)
intersect(NLRI,names(regulons$SOX21))
intersect(CAMC,names(regulons$SOX21))
intersect("NLGN1",names(regulons$SOX21))

head(regulons$SOX21)
intersect(NLRI,names(regulons$SP1))
intersect(CAMC,names(regulons$SP1))
intersect("NLGN1",names(regulons$SP1))


head(regulons$SP1)


g<-tni.graph(tpm.log1, regulatoryElements = final.tfs)
library(RedeR)
rdp <- RedPort()
calld(rdp,checkcalls=TRUE)
addGraph(rdp, g, layout=NULL)
addLegend.color(rdp, g, type="edge")
addLegend.shape(rdp, g)
relax(rdp, ps = TRUE)


deg <- read.delim("../RNA-seq/metastases_versus_primary.deseq.sig.xls")
t1 <- intersect(gene.dis.hypo$Var1[1:20],intersect(deg$gene_symbol,NLRI))
t2 <- intersect(gene.dis.hyper$Var1[1:20],intersect(deg$gene_symbol,NLRI))

t1 <- intersect(gene.dis.hypo$Var1[1:100],NLRI)
t2 <- intersect(gene.dis.hyper$Var1[1:100],NLRI)

t1 <- gene.dis.hypo$Var1[1:20]
t2 <- gene.dis.hyper$Var1[1:20]


for(i in colnames(test)){
  
  #print(intersect(NLRI,names(regulons[[i]])))
  #print(intersect(CAMC,names(regulons[[i]])))
  inter <- intersect("NLGN1",names(regulons[[i]]))
  if(length(inter)!=0){
    print(i)
    print("*****************")
    print(inter)
  }
}
