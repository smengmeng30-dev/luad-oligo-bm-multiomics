
library(readxl)


# A) TIL infiltration -----------------------------------------------------

gsva_matrix_norm <- read.delim("ssgsea.immune_score_normlization.tsv",check.names = F)


scores=read.table("estimate_score.gct",skip = 2,header = T,check.names = F)
rownames(scores)=scores[,1]
scores=scores[,3:ncol(scores)]
colnames(scores) <- gsub("^X","",colnames(scores))
scores <- as.data.frame(t(scores))
scores <- scores[colnames(gsva_matrix_norm),]

identical(rownames(scores),colnames(gsva_matrix_norm))

immuneCluster <- read.delim("ConsensusClusterPlus/cluster.txt")
identical(rownames(immuneCluster),colnames(gsva_matrix_norm))

group <- read.delim("../clinical.group.txt",stringsAsFactors = F)
rownames(group) <- group$sampleid.rna
group <- group[colnames(gsva_matrix_norm),]
group$SID <- paste(group$PID,group$group,sep="_")

identical(rownames(group),colnames(gsva_matrix_norm))

colnames(gsva_matrix_norm) <- group$SID

library(ComplexHeatmap)


pdf(file="../paper1.Figure/immuneRelated.heatmap.pdf",width=12,height=5)

p <- Heatmap(gsva_matrix_norm,cluster_rows = T,
             show_column_names  = T,
             cluster_columns = cluster_within_group(gsva_matrix_norm,immuneCluster$cluster),
             #column_split = 2,
             #column_km = 2,
             col = circlize::colorRamp2(c(0,0.25,0.5,0.75,1), c('#09386c','#7ab5d5', "white",'#d8634f','#6a011f')),
             #clustering_method_columns = 'ward.D2',
             #clustering_distance_rows='euclidean',
             #clustering_distance_columns='euclidean',
             #clustering_method_rows = 'ward.D2',
             rect_gp = gpar(col = 'white'),
             name='enrichment\nscore',
             #column_split = c(rep('SCLC-LCC',14),rep('SCLC-NSCLC',18)),
             #column_gap=unit(3,'mm'),
             #border = TRUE,
             show_column_dend = F,
             top_annotation = columnAnnotation(group=group$group,
                                               #Time=group$time,
                                               cluster=immuneCluster$cluster,
                                               StromalScore=scores$StromalScore,
                                               ImmuneScore=scores$ImmuneScore,
                                               ESTIMATEScore=scores$ESTIMATEScore,
                                               
                                               col=list(group=c(Metastases="#E21D31",
                                                                Primary="#0C7094"),
                                                        #Time=c(asynchronous='#00468BFF',synchronous='#6a011f',unknown='grey'),
                                                        StromalScore=circlize::colorRamp2(c(-2000, -1000, 0, 1000, 2000),
                                                                                          c("#FCFBFD","#DADAEB", "#9E9AC8","#6A51A3","#3F007D")),
                                                        ImmuneScore=circlize::colorRamp2(c(-2000, -1000, 0, 1000, 2000),
                                                                                         c("#F7FCF5",'#C7E9C0', "#74C476",'#238B45','#00441B')),
                                                        ESTIMATEScore=circlize::colorRamp2(c(-2000, -1000, 0, 1000, 2000),
                                                                                           c('#F7FBFF','#C6DBEF', "#6BAED6",'#2171B5','#08306B')),
                                                        cluster=c('1'='#70ad9e','2'='#d1a7cf','3'='#8ad5f4')),
                                               simple_anno_size = unit(3, "mm")
             ))
print(p)
dev.off()



# B) each T cell ----------------------------------------------------------


gsva_matrix <- read.delim("ssgsea.immune_score.tsv",check.names = F)

gsva_matrix1 <- as.data.frame(gsva_matrix)
gsva_matrix1$cell <- rownames(gsva_matrix1)
gsva_matrix1 <- melt(gsva_matrix1)

gsva_matrix1 <- merge(gsva_matrix1,group,by.x='variable',by.y='row.names')

gsva_matrix1$type <- paste(gsva_matrix1$time,gsva_matrix1$group,sep='\n')


p <- list()
for(i in unique(gsva_matrix1$cell)){
  tmp <- gsva_matrix1[gsva_matrix1$cell==i,]
  res <- wilcox.test(tmp$value[tmp$group=='Metastases'],tmp$value[tmp$group=='Primary'])
  
  if(!is.na(res$p.value) & res$p.value<0.05){
    p[[i]] <- ggboxplot(tmp,x='group',y='value',color = 'group',palette = c("#0C7094","#E21D31"),add='jitter',xlab = '',ylab = i)+stat_compare_means()+theme_bw()+theme(legend.position = 'none',panel.grid = element_blank())
  }
}

library(gridExtra)
p1 <- marrangeGrob(p,nrow=2,ncol=2,top = '')
ggsave("../paper1.Figure/immunescore.boxplot.pdf",p1,width = 5,height = 5)



# C) estimate score -------------------------------------------------------


scores2 <- merge(scores,group,by.x='row.names',by.y='row.names')
scores2 <- melt(scores2,measure.vars = 2:4)

p <- list()
for(i in unique(scores2$variable)){
  tmp <- scores2[scores2$variable==i,]
  res <- wilcox.test(tmp$value[tmp$group=='Metastases'],tmp$value[tmp$group=='Primary'])
  p[[i]] <- ggboxplot(tmp,x='group',y='value',color = 'group',palette = c("#0C7094","#E21D31"),add='jitter',xlab = '',ylab = i)+stat_compare_means()+theme_bw()+theme(legend.position = 'none',panel.grid = element_blank())
}

library(gridExtra)
p1 <- p[['ESTIMATEScore']]
ggsave("../paper1.Figure/estimatescore.boxplot.pdf",p1,width = 3,height = 3)
p2 <- p[['ImmuneScore']]
ggsave("../paper1.Figure/ImmuneScore.boxplot.pdf",p2,width = 3,height = 3)





snv <- read_xlsx("../1021/all.snv.filter.recheck.xlsx",sheet = 1)
snv <- snv[snv$recheck=='Y',]
TMB <- snv[snv$isTMB=='Y',]
TMB <- as.data.frame(table(TMB$Tumor_Sample_Barcode))


dim(gsva_matrix1)

gsva_matrix11 <- gsva_matrix1[gsva_matrix1$cell %in% c("Natural killer cell","Activated B cell","Type 1 T helper cell","Immature  B cell"),]
gsva_matrix11 <- dcast(gsva_matrix11,sampleid.1021+sampleid.rna~cell,value.var = 'value')

output <- merge(TMB,gsva_matrix11,by.x='Var1',by.y='sampleid.1021')

scores3 <- scores2[scores2$variable %in% c("ImmuneScore"),]
scores3 <- dcast(scores3,PID+group+sampleid.rna~variable,value.var = 'value')

output <- merge(output,scores3,by='sampleid.rna')



## collagen胶原蛋白比较
library(GSVA)
collagen <- read.delim("../RNA-seq/collagen.txt")
gene.list <- list()
gene.list[['collagen']] <- collagen$name

expr <- read.delim("../RNA-seq/geneSymbol_tpm.xls",check.names = F,row.names = 1)

expr <- log(expr)

gsva_matrix<- gsva(as.matrix(expr), gene.list,method='ssgsea',kcdf='Gaussian')

gsva_matrix <- as.data.frame(t(gsva_matrix))

output <- merge(output,gsva_matrix,by.x='sampleid.rna',by.y='row.names')

output <- output[,c("PID","group","Freq","ImmuneScore","Activated B cell","Immature  B cell","Natural killer cell","Type 1 T helper cell","collagen")]

colnames(output)[3] <- 'TMB'

write.table(output,file = '../paper1.Figure/tmb_immunescore_collagen.txt',quote = F,sep='\t',row.names = F)

# D) TIL levels -----------------------------------------------------------


# DEG 火山图 -----------------------------------------------------------------


deg <- read.delim("metastases_versus_primary.deseq.xls",stringsAsFactors = F)
deg$regulated[deg$log2FoldChange>=1 & deg$pvalue<0.05] <- "up"
deg$regulated[deg$log2FoldChange <= -1 & deg$pvalue<0.05] <- "down"
deg$regulated[is.na(deg$regulated)]<- "not"
deg$pvalue <- -log10(deg$pvalue)
label <- deg[deg$regulated!="not",]
top_down <- label[order(abs(label$log2FoldChange),label$pvalue,decreasing = T),]
label <- top_down[c(1:20,which(top_down$gene_symbol=='NLGN1')),]

ggplot(deg,aes(x=log2FoldChange,y=pvalue))+geom_point(aes(color=regulated))+
  scale_color_manual(values=c(alpha("#00468BFF", 0.7),alpha("black", 0.1),alpha("#AD002AFF", 0.7)))+
  #scale_x_continuous(breaks = c(-27,-10,0,10,27),limits = c(-27,27),expand=c(0,0))+
  theme_bw()+theme(panel.grid = element_blank())+
  geom_vline(xintercept=c(-1,1),linetype=2)+
  geom_hline(yintercept=-log10(0.05),linetype=2)+
  ggrepel::geom_label_repel(dat=label,aes(x=log2FoldChange,
                                          y=pvalue,label=gene_symbol),size=2,label.size=unit(0.03,'mm'),
                            max.overlaps=30)+ylab("-log10(pvalue)")



ggsave("../paper1.Figure/deg.plot.pdf",width = 6,height = 5)



# GSEA --------------------------------------------------------------------
library(clusterProfiler)
fc <- read.delim("metastases_versus_primary.deseq.xls")
gene.list <- fc$log2FoldChange
names(gene.list) <- fc$ensembl

gene.list <- sort(gene.list,decreasing = T)

library(msigdbr)

GOdb = msigdbr(species = "Homo sapiens", category = "C5",subcategory = "GO:BP")
GOdb <- GOdb[,c('gs_name','ensembl_gene')]
GOdb$gs_name <- gsub("^GOBP_","",GOdb$gs_name)
GOdb$gs_name <- tolower(gsub("_"," ",GOdb$gs_name))
colnames(GOdb) <- c("TREM","GENE")

keggdb = msigdbr(species = "Homo sapiens", category = "C2",subcategory = "KEGG")
keggdb <- keggdb[,c('gs_name','ensembl_gene')]
keggdb$gs_name <- gsub("^KEGG_","",keggdb$gs_name)
keggdb$gs_name <- tolower(gsub("_"," ",keggdb$gs_name))
colnames(keggdb) <- c("TREM","GENE")

hallmarkdb = msigdbr(species = "Homo sapiens", category = "H")
hallmarkdb <- hallmarkdb[,c('gs_name','ensembl_gene')]
hallmarkdb$gs_name <- gsub("^HALLMARK_","",hallmarkdb$gs_name)
hallmarkdb$gs_name <- tolower(gsub("_"," ",hallmarkdb$gs_name))
colnames(hallmarkdb) <- c("TREM","GENE")



gogmt <- clusterProfiler::GSEA(gene.list, TERM2GENE = GOdb,
                               nPermSimple = 100000,
                               pvalueCutoff = 0.01,eps=1e-30)

write.table(as.data.frame(gogmt@result),file='../paper1.Figure/meta_vs_pri.go.gsea.txt',quote = F,sep='\t',row.names = F)

up.id <- gogmt@result$ID[gogmt@result$NES>0][1:5]
down.id <- gogmt@result$ID[gogmt@result$NES<0][1:5]



enrichplot::gseaplot2(gogmt, geneSetID = c(up.id,down.id),
                      color = c("#00468BFF","#ED0000FF","#42B540FF","#0099B4FF","#925E9FFF","#FDAF91FF","#AD002AFF","#ADB6B6FF","#1B1919FF","#00468BCC"),
                      pvalue_table = F, subplots = 1:2
)

ggsave("../paper1.Figure/go.gsea.pdf",width = 5,height = 4)


kegggmt <- clusterProfiler::GSEA(gene.list, TERM2GENE = keggdb,
                                 nPermSimple = 100000,
                                 pvalueCutoff = 1,eps=1e-20)

dat <- kegggmt@result
dat <- dat[dat$pvalue<0.05 & dat$qvalues<0.1,]

write.table(dat,file='../paper1.Figure/meta_vs_pri.kegg.gsea.txt',quote = F,sep='\t',row.names = F)

enrichplot::gseaplot2(kegggmt, geneSetID = dat$ID,
                      color = c("#00468BFF","#ED0000FF","#42B540FF","#0099B4FF","#925E9FFF","#FDAF91FF","#AD002AFF","#ADB6B6FF"),
                      pvalue_table = F, subplots = 1:2
)

ggsave("../paper1.Figure/kegg.gsea.pdf",width = 5,height = 4)

hgmt <- clusterProfiler::GSEA(gene.list, TERM2GENE = hallmarkdb,
                               nPermSimple = 100000,
                               pvalueCutoff = 1,eps=1e-20)

dat <- hgmt@result
dat <- dat[dat$pvalue<0.05 & dat$qvalues<0.05,]

write.table(dat,file='../paper1.Figure/meta_vs_pri.hallmark.gsea.txt',quote = F,sep='\t',row.names = F)

enrichplot::gseaplot2(hgmt, geneSetID = dat$ID,
                      color = c("#00468BFF","#ED0000FF","#42B540FF"),
                      pvalue_table = F, subplots = 1:2
)
ggsave("../paper1.Figure/hallmark.gsea.pdf",width = 5,height = 4)

#devtools::install_github("junjunlab/GseaVis")
library(GseaVis)
library(org.Hs.eg.db)


go.path <- "07_GSEA/metastases_versus_primary/GO/metastases_versus_primary.GO.Gsea.1686327707512/"
kegg.path <- "07_GSEA/metastases_versus_primary/kegg/metastases_versus_primary.kegg.Gsea.1686326944430/"
hallmark.path <- "07_GSEA/metastases_versus_primary/hallmark/metastases_versus_primary.hallmark.Gsea.1686326929899/"

##################### go

go.res <- readGseaFile(filePath = go.path)
go.id <- go.res$meta$ID[go.res$meta$p.adjust<0.05]
# new style plot
# plot
gseaNb(filePath = go.path,geneSetID = go.id,newGsea=F,
       curveCol = c("#00468BFF","#ED0000FF","#42B540FF","#0099B4FF","#925E9FFF","#FDAF91FF","#AD002AFF"),
       subPlot = 2)

ggsave(filename = '../paper1.Figure/go.gsea.plot.pdf',width = 8,height = 4)

##################### hallmark

h.res <- readGseaFile(filePath = hallmark.path)
h.id <- h.res$meta$ID[h.res$meta$p.adjust<0.05]

gseaNb(filePath = hallmark.path,
         geneSetID = h.id,newGsea=F,
       curveCol = c("#00468BFF","#ED0000FF","#42B540FF","#0099B4FF","#925E9FFF","#FDAF91FF","#AD002AFF","#ADB6B6FF","#1B1919FF"),
       subPlot = 2)

ggsave(filename = '../paper1.Figure/hallmark.gsea.plot.pdf',width = 8,height = 4)


##################### kegg

kegg.res <- readGseaFile(filePath = kegg.path)
kegg.id <- kegg.res$meta$ID[kegg.res$meta$p.adjust<0.05]

gseaNb(filePath = kegg.path,
       geneSetID = kegg.id,curveCol = c("#00468BFF","#ED0000FF","#42B540FF","#0099B4FF","#925E9FFF","#FDAF91FF"),
       subPlot = 2)

ggsave(filename = '../paper1.Figure/kegg.gsea.plot.pdf',width = 8,height = 4)



# quanTIseq ---------------------------------------------

library(ComplexHeatmap)

quanTIseq <- read.delim("../paper1.Figure/quanTIseq.txt",row.names = 1)
clinical <- read.delim("../clinical.group.txt")
rownames(clinical) <- clinical$sampleid.rna
clinical <- clinical[colnames(quanTIseq1),]

col <- c("#0C7094","#E21D31")
names(col) <- unique(clinical$group)

quanTIseq <- quanTIseq[,-1]

quanTIseq1 <- as.matrix(t(scale(quanTIseq)))


pdf('../paper1.Figure/quanTIseq.heatmap.pdf',width = 6,height = 4)

Heatmap(quanTIseq1,
        cluster_columns = cluster_within_group(quanTIseq1,clinical$group),
        top_annotation = columnAnnotation(group=clinical$group,col=list(group=col)),show_column_names = F)

dev.off()


quanTIseq2 <- merge(quanTIseq,clinical,by.x='row.names',by.y='row.names')

quanTIseq2 <- melt(quanTIseq2,measure.vars = 2:12)

quanTIseq2 <- dcast(quanTIseq2,PID+variable~group,value.var = 'value')

ggpaired(quanTIseq2,cond1 = 'Primary',cond2 = 'Metastases',
         color = 'condition',
         
         fill = 'condition',palette = col,line.color = 'grey',
)+scale_fill_manual(values = alpha(col,0.3))+
  facet_wrap(~variable,nrow=1)+stat_compare_means(paired = T,aes(label = paste0("p = ", ..p.format..)))+
  ### set background of facet labels 
  theme(strip.background = element_rect(fill='white',colour = 'white'))

ggsave(filename = '../paper1.Figure/quanTIseq.boxplot.pdf',width = 9,height = 4)

# immune signature --------------------------------------------------------

library(tidyr)
library(ComplexHeatmap)

immunescore1 <- read.delim('../画图需求v2/input/immuneSignature.heatmap.txt',check.names = F,row.names = 1)

clinical <- read.delim("../clinical.group.txt")
rownames(clinical) <- clinical$sampleid.rna
clinical <- clinical[colnames(immunescore1),]

col <- c("#0C7094","#E21D31")
names(col) <- unique(clinical$group)

pdf('../paper1.Figure/immuneSignature.heatmap.pdf',width = 8,height = 4)
Heatmap(immunescore1,
        cluster_columns = cluster_within_group(immunescore1,clinical$group),
        top_annotation = columnAnnotation(group=clinical$group,col=list(group=col)),show_column_names = F)
dev.off()


library(dplyr)

immunescore2 <- read.delim("../paper1.Figure/immunescore.boxplot.txt")
p_values <- immunescore2 %>%
  group_by(variable) %>%
  summarise(p_value = wilcox.test(value ~ group,paired=T)$p.value)

p_values <- as.data.frame(p_values)

immunescore2 <- dcast(immunescore2,PID+variable~group,value.var = 'value')

ggpaired(immunescore2,cond1 = 'Primary',cond2 = 'Metastases',
         color = 'condition',
         xlab = '',ylab = 'Mean expression',
         fill = 'condition',palette = col,line.color = 'grey',
)+scale_fill_manual(values = alpha(col,0.3))+
  facet_wrap(~variable,nrow=1)+
  ### set background of facet labels 
  theme(strip.background = element_rect(fill='white',colour = 'white'))+
  scale_y_log10()

ggsave(filename = '../paper1.Figure/immuneSignature.boxplot.pdf',width = 9,height = 4)


library(gg.gap)
gg.gap(plot=p,
         segments=list(c(50,100),c(200,800)),
         ylim=c(0,1000))



# cell death --------------------------------------------------------------

library(GSVA)
CellDeath_gtm_ref <- CellDeath_gtm_ref[CellDeath_gtm_ref$gene!='',]

cellDeath.list <- split(as.matrix(CellDeath_gtm_ref)[, 2], CellDeath_gtm_ref[, 1])

eset_stad <- read.delim("All.TPM_geneSymbol.txt",row.names = 1,check.names = F)
gsva_matrix <- gsva(
    as.matrix(eset_stad),
    cellDeath.list,
    method = 'ssgsea',
    kcdf = 'Gaussian',
    abs.ranking = TRUE
  )

head(gsva_matrix)

gsva_matrix1 <- as.data.frame(gsva_matrix)
gsva_matrix1$Signature <- rownames(gsva_matrix1)
gsva_matrix1 <- melt(gsva_matrix1)

group <- read.delim("../clinical.group.txt")

gsva_matrix1 <- merge(gsva_matrix1,group,by.x='variable',by.y='sampleid.rna')


library(ggpubr)

ggboxplot(gsva_matrix1,x='group',y='value',
          color='group',palette = c("#0C7094","#E21D31"),
          facet.by = 'Signature',add='jitter',xlab = '')+
  stat_compare_means(aes(group=group))+facet_wrap(~Signature,scales = 'free',ncol=4,nrow=4)+theme(legend.position = 'none')


ggsave("../paper1.Figure/CellDeath.pri_vs_meta.pdf",width = 10,height = 10)
