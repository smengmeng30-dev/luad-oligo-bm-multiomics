clinical <- read.delim("../clinical.group.txt")


snv <- read_excel("all.snv.filter.recheck.xlsx")
snv <- snv[snv$recheck=='Y',]


snv$flag <- paste(snv$Chromosome,snv$Start_Position,snv$End_Position,snv$Reference_Allele,snv$Tumor_Seq_Allele2,sep=":")

snv <- merge(snv,clinical,by.x='Tumor_Sample_Barcode',by.y='sampleid.1021')

#BiocManager::install('sigminer')
#BiocManager::install('quadprog')
library(maftools)
library(sigminer)
library(NMF)
library(BSgenome.Hsapiens.UCSC.hg19)
snv1 <- snv
snv1$Tumor_Sample_Barcode <- snv1$group
snv1 <- unique(snv1[,c("Tumor_Sample_Barcode","Chromosome", "Start_Position","End_Position","Reference_Allele","Tumor_Seq_Allele2","Variant_Classification",'Variant_Type',"Hugo_Symbol","cHGVS","pHGVS")])

all.maf  <- read.maf(snv1,vc_nonSyn = unique(snv1$Variant_Classification))
all.signature <- COSMIC_V3(all.maf,mode = "legacy")

sim <- all.signature[[1]]
absolute_all <- as.data.frame(all.signature[[2]])
relative_all <- as.data.frame(all.signature[[3]])
relative_contri <- all.signature[[4]]
absolute_contri <- all.signature[[5]]



rownames(absolute_all) <- absolute_all$sample
absolute_all <- absolute_all[,-1]
absolute_all <- as.data.frame(t(absolute_all))

write.table(absolute_all,file = '../paper1.Figure/all.signature.absolute.bygroup.txt',quote = F,sep = '\t',row.names = T)

show <- absolute_all[rowSums(absolute_all)!=0,]

absolute_all1 <- absolute_all[rownames(show),]
absolute_all1$sig <- rownames(absolute_all1)
absolute_all1 <- melt(absolute_all1)

ggbarplot(absolute_all1,x='variable',y='value',fill='sig',color = 'sig',
          xlab = '',
          position = position_fill(),palette = c("#E55C27","#F1A646","#21579D","#9BBADD","#F2B1B0"))+theme_bw()+theme(panel.grid = element_blank())

ggsave("../paper1.Figure/pri_vs_meta.signature.pdf",width = 4,height = 4)

######### shared/specific signature
snv2 <- snv


snv.final <- dcast(snv2,flag~group,fun.aggregate = length)

snv.final <- snv.final[snv.final$Metastases!=0 | snv.final$Primary!=0,]

snv.final$type[snv.final$Metastases>0 & snv.final$Primary>0] <- 'Trunk'
snv.final$type[snv.final$Metastases>0 & snv.final$Primary==0] <- 'Metastases specific'
snv.final$type[snv.final$Metastases==0 & snv.final$Primary>0] <- 'Primary specific'

snv2 <- unique(snv2[,c('flag',"Chromosome", "Start_Position","End_Position","Reference_Allele","Tumor_Seq_Allele2","Variant_Classification",'Variant_Type',"Hugo_Symbol","cHGVS","pHGVS")])

snv2 <- merge(snv2,snv.final[,c('flag',"type")],by='flag')

snv2$Tumor_Sample_Barcode <- snv2$type

snv2 <- unique(snv2[,c('Tumor_Sample_Barcode',"Chromosome", "Start_Position","End_Position","Reference_Allele","Tumor_Seq_Allele2","Variant_Classification",'Variant_Type',"Hugo_Symbol","cHGVS","pHGVS")])

maf.bygroup  <- read.maf(snv2,vc_nonSyn = unique(snv2$Variant_Classification))

maf.signature <- COSMIC_V3(maf.bygroup,mode = "legacy")

sim <- maf.signature[[1]]
absolute_all <- as.data.frame(maf.signature[[2]])
relative_all <- as.data.frame(maf.signature[[3]])
relative_contri <- maf.signature[[4]]
absolute_contri <- maf.signature[[5]]


rownames(absolute_all) <- absolute_all$sample
absolute_all <- absolute_all[,-1]
absolute_all <- as.data.frame(t(absolute_all))

write.table(absolute_all,file = '../paper1.Figure/all.signature.absolute.byTrunkBranch.txt',quote = F,sep = '\t',row.names = T)

show <- absolute_all[rowSums(absolute_all)!=0,]

absolute_all1 <- absolute_all[rownames(show),]
absolute_all1$sig <- rownames(absolute_all1)
absolute_all1 <- melt(absolute_all1)

ggbarplot(absolute_all1,x='variable',y='value',fill='sig',color = 'sig',
          xlab = '',
          position = position_fill(),palette = c("#E55C27","#F1A646","#21579D","#9BBADD","#F2B1B0","#7C65A4","#E37AA9"))+theme_bw()+theme(panel.grid = element_blank())

ggsave("../paper1.Figure/trunk&branch.signature.pdf",width = 4,height = 4)
