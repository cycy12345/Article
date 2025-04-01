setwd("D:/WORK/Project/32糖尿病视网膜病变表观遗传/")
rm(list = ls())
library(GEOquery)
library(Biobase)
library(oligo)
library(pd.hg.u133.plus.2)
library(pd.clariom.d.human)
library(CyDataPro)
dir.create("data")
BiocManager::install(c( 'oligo' ),ask = F,update = F)
Me <- fread("data/GSE121820_Differential_DNA_methylation.csv")
Me <- Me %>% select(`Gene Symbol`,`p-value(region)`,chromosome,`transcript start`)
colnames(Me) <- c("Gene","Pvalue","Chr","BP")
Mee <- Me %>% distinct(Gene,.keep_all = T)
Mee<- Mee %>% filter(Chr !="chr5_h2_hap1")
Mee<- Mee %>% filter(Chr !="chr6_cox_hap1")
Mee<- Mee %>% filter(Chr !="chrX")
Mee<- Mee %>% filter(Chr !="chrY")
Mee<- Mee %>% filter(Chr !="chrM")

data <- table(Mee$Pvalue) %>% as.data.frame.array() %>% as.data.frame()

data$Pvalue <- rownames(data)
colnames(data)[1]<- "GeneNumber"
data$Pvalue <- as.numeric(data$Pvalue) %>% format(scientific=T,digits = 3)
str(data)
data$Type <- ifelse(as.numeric(data$Pvalue )< 0.001,"Yes","No")

data<- data %>%  arrange(Pvalue)
data$Pvalue
data$Pvalue<- factor(data$Pvalue,levels = rev(c("8.79000e-06","6.41533e-04","1.93339e-03","4.13042e-04",
                                                "4.04253e-04","1.88945e-03","3.77889e-04","1.90702e-03")))
dir.create("DEG/Methy")
p<-ggplot(data,aes(x=GeneNumber,y=Pvalue,fill=Type))+
  geom_bar(stat = "identity")+
  scale_fill_manual(values = c("gray","red"))+
  labs(x="Gene Number",y="Pvalue",fill="pvalue < 1e-3")+
  xlim(0,5500)+
  geom_text(aes(label=GeneNumber),size=5)+
  theme_bw()+
  geom_hline(yintercept=3.5,size=2,
             lty="dashed")+
  mytheme
p
save_plot(p,filename = "DEG/Methy/GeneNumber_barplot",style = "g",width = 8,height = 7.5)
tmp <- Mee %>% filter(Pvalue<1e-3)
length(unique(tmp$Gene))
fwrite(tmp,file = "DEG/Methy/Methy_DEG_Gene.csv")
dir='./data/GSE189005_RAW///'

od=getwd()
setwd(dir)
celFiles <- list.celfiles(listGzipped = T)
celFiles
affyRaw <- read.celfiles( celFiles )
setwd("../../")
eset <- rma(affyRaw)

dat=exprs(eset) %>% as.data.frame()
class(dat)
dat <- dat %>% rownames_to_column("ID")
boxplot(dat)
group_list = substring(celFiles,12,14)

gpl <- fread("../自身免疫病_骨关节炎_MR20/data/GEO/GPL23126-131.txt",skip = 19)
gpl<- gpl %>% filter(unigene !="---")
edit(GEODownload)
colnames(gpl)
Symbol_col = "gene_assignment"
Pro_ID <- select(gpl, c("ID", Symbol_col))
colnames(Pro_ID) = c("ID", "Gene_Symbol")
Pro_ID$ID <- as.character(Pro_ID$ID)
Pro_ID <- Pro_ID[Pro_ID$Gene_Symbol != "", ]
Pro_ID$Gene_Symbol <- sapply(Pro_ID$Gene_Symbol, function(x) {
  unlist(str_split(x, " // "))[2]
})
exp_symbol <- merge(Pro_ID, dat, by = "ID") %>% na.omit()
exp_unique <- avereps(exp_symbol[, -c(1, 2)], ID = exp_symbol$Gene_Symbol) %>% 
  as.data.frame()
exp_log<- dectect_log(exp_unique)     
colnames(exp_log)<-substring(colnames(exp_log),1,10)
Group <- ifelse(group_list=="T2D","T2DR","Control")
Group <- factor(Group,levels = c("Control","T2DR"))
pvalue <- 0.05
logFC_cutoff <- 1
design <- model.matrix(~Group)
exper <- normalizeBetweenArrays(exp_log)
boxplot(exper,outline=F,notch=T,color=Group,las=2)
fit <- lmFit(exper,design)
efit <-eBayes(fit,trend = F)
DEG <- topTable(efit,coef =2,n=Inf)
K1 <- (DEG$P.Value<pvalue)&(DEG$logFC < -logFC_cutoff)
K2 <- (DEG$P.Value<pvalue)&(DEG$logFC > logFC_cutoff)
DEG$change <- ifelse(K1,"Down",ifelse(K2,"Up","not"))
table(DEG$change)
dir.create("DEG/mRNA",recursive = T)
fwrite(DEG,file = "DEG/mRNA/DEG_result.csv",row.names = T)
save(exper,Group,file = "data/GSE189005_mRNA.Rdata")

devtools::install_github("BioSenior/ggVolcano")
library(ggVolcano)
library(ggsci)
library(randomcoloR)
library(CyDataPro)
library(pheatmap)
library(Cycolors)
Disc9(9)
cc <- distinctColorPalette(10)
randomColor(count = 10) 

"#8CE472" "#7BDDC3" "#9787D7" "#8EBAD5" "#D57A66" "#DCD266" "#A34DE2" "#E0ADCE" "#D2D9BD" "#E15CB8"
rm(list = ls())
DEG <- read.csv("DEG/LncRNA///DEG_result.csv",row.names = 1)
K1 <- (DEG$P.Value<0.05)&(DEG$logFC < -0.7)
K2 <- (DEG$P.Value<0.05)&(DEG$logFC > 0.7)
DEG$change <- ifelse(K1,"Down",ifelse(K2,"Up","not"))
table(DEG$change)
deg <- DEG
deg$Gene <- rownames(deg)
colnames(deg)
"#9787D7",'#E15CB8'
data <- add_regulate(deg, log2FC_name = "logFC",
                     fdr_name = "P.Value",log2FC = 0.7, fdr = 0.05)
p<-ggvolcano(data, x = "log2FoldChange", y = "padj",y_lab = "-Log10Pvalue",
          x_lab = "Log2FC",
          fills = c("#9787D7","#b4b4d8","#E15CB8"),
          colors = c("#9787D7","#b4b4d8","#E15CB8"),
          label = "Gene", label_number = 20, output = FALSE)+
  mytheme
  
p
save_plot(p,filename = "DEG/LncRNA//VolcanoPlot",style = "g")
load("data/GSE185011_LncRNA.Rdata")
Dexp <- exper[rownames(deg[deg$change !="not",]),]
annocol <- data.frame(group=Group) 
ann_colors=list(group=c(Control='#8EBAD5',T2DR='#D57A66'))
rownames(annocol) <- colnames(Dexp)
annocol <- arrange(annocol,group)
heatmap_exp <- Dexp[,rownames(annocol)]
p <- pheatmap(heatmap_exp, show_colnames = F, show_rownames = F,
              scale = "row",
              cluster_cols = F,
              annotation_col = annocol,
              breaks = seq(-3, 3, length.out = 100),annotation_colors = ann_colors,
              color = colorRampPalette(c('#9787D7','white','#E15CB8'))(100),labels_col = "x",
              fontsize = 10)
p
save_plot(p,filename = "DEG/LncRNA///DEG_heatmap",style = "xx")
fwrite(DEG,file = "DEG/LncRNA//DEG_result.csv",row.names = T)
Df <- deg %>% filter(change !="not")
DM <- Mee %>% filter(Pvalue < 0.0001)

intersect(Df$Gene,DM$Gene) %>% length()
library(pd.mirna.4.0)
library(limma)
rm(list = ls())
dir='./data/GSE189002_RAW///'

od=getwd()
setwd(dir)
celFiles <- list.celfiles(listGzipped = T)
celFiles
affyRaw <- read.celfiles( celFiles )
setwd("../../")
eset <- rma(affyRaw)

dat=exprs(eset) %>% as.data.frame()
class(dat)
boxplot(dat)
dat <- dat %>% rownames_to_column("ID")


group_list = substring(celFiles,12,14)

gpl <- fread("../自身免疫病_骨关节炎_MR20/data/GEO/GPL19117-74051.txt",skip = 17)
table(gpl$`Species Scientific Name`)
gpl<- gpl %>% filter(`Species Scientific Name` =="Homo sapiens")

edit(GEODownload)
colnames(gpl)
Symbol_col = "miRNA_ID"
Pro_ID <- dplyr::select(gpl, c("ID", Symbol_col))
colnames(Pro_ID) = c("ID", "Gene_Symbol")
Pro_ID$ID <- as.character(Pro_ID$ID)
Pro_ID <- Pro_ID[Pro_ID$Gene_Symbol != "", ]
Pro_ID$Gene_Symbol <- sapply(Pro_ID$Gene_Symbol, function(x) {
  unlist(str_split(x, " // "))[2]
})
exp_symbol <- merge(Pro_ID, dat, by = "ID") %>% na.omit()
exp_unique <- avereps(exp_symbol[, -c(1, 2)], ID = exp_symbol$Gene_Symbol) %>% 
  as.data.frame()
exp_log<- dectect_log(exp_unique)     
colnames(exp_log)<-substring(colnames(exp_log),1,10)
Group <- ifelse(group_list=="T2D","T2DR","Control")
Group <- factor(Group,levels = c("Control","T2DR"))
pvalue <- 0.05
logFC_cutoff <- 0.7
design <- model.matrix(~Group)
exper <- normalizeBetweenArrays(exp_log)
boxplot(exper,outline=F,notch=T,color=Group,las=2)
fit <- lmFit(exper,design)
efit <-eBayes(fit,trend = F)
DEG <- topTable(efit,coef =2,n=Inf)
K1 <- (DEG$P.Value<pvalue)&(DEG$logFC < -logFC_cutoff)
K2 <- (DEG$P.Value<pvalue)&(DEG$logFC > logFC_cutoff)
DEG$change <- ifelse(K1,"Down",ifelse(K2,"Up","not"))
table(DEG$change)
dir.create("DEG/miRNA",recursive = T)
fwrite(DEG,file = "DEG/miRNA/DEG_result.csv",row.names = T)
save(exper,Group,file = "data/GSE189002_miRNA.Rdata")


library(GEOquery)
library(clusterProfiler)
library(org.Hs.eg.db)
urld <- "https://www.ncbi.nlm.nih.gov/geo/download/?format=file&type=rnaseq_counts"
path <- paste(urld, "acc=GSE185011", "file=GSE185011_raw_counts_GRCh38.p13_NCBI.tsv.gz", sep="&");
path
options(timeout = 1000)
tbl <- data.table::fread(path, header=T, colClasses="integer")
tbl<- as.data.frame(tbl)
gset = getGEO("GSE185011", destdir = './data/', getGPL = F) 
pd = pData(gset[[1]])
pd<-pd %>% filter(str_detect(title,"HC") | str_detect(title,"DR"))
Group<- ifelse(str_detect(pd$title,"HC"),"Control","T2DR")
id <- bitr(tbl$GeneID, fromType = "ENTREZID", 
           toType = c('SYMBOL','GENETYPE'),
           OrgDb = org.Hs.eg.db)
colnames(id)[1]<- "GeneID"
class(tbl)
tbll <-  merge(tbl, id, by = "GeneID")
table(tbll$GENETYPE)
exp_symbol<- tbll %>% filter(GENETYPE =="ncRNA")
exp_symbol<- exp_symbol[,-28]
exp_unique <- avereps(exp_symbol[, -c(1, 27)], ID = exp_symbol$SYMBOL) %>% 
  as.data.frame()
exp <- exp_unique[,rownames(pd)]
keep <- rowSums( exp >= 2 ) >= 2
expr <- exp[keep, ]
pvalue <- 0.05
logFC_cutoff <- 0.7
design <- model.matrix(~Group)

exp_log<- dectect_log(expr)
boxplot(exp_log,outline=F,notch=T,color=Group,las=2)
exper <- normalizeBetweenArrays(exp_log)

boxplot(exper,outline=F,notch=T,color=Group,las=2)
fit <- lmFit(exper,design)
efit <-eBayes(fit,trend = F)
DEG <- topTable(efit,coef =2,n=Inf)
K1 <- (DEG$P.Value<pvalue)&(DEG$logFC < -logFC_cutoff)
K2 <- (DEG$P.Value<pvalue)&(DEG$logFC > logFC_cutoff)
DEG$change <- ifelse(K1,"Down",ifelse(K2,"Up","not"))
table(DEG$change)

dir.create("DEG/LncRNA")
save(exper,Group,file = "data/GSE185011_LncRNA.Rdata")
"#71ABB6" "#5EB75B" "#4C8BC0" "#E83133" "#F0B57D" "#C75C64" "#D3E1AE" "#4B5AA1" "#A360AD"
deg <- DEG 
deg$Gene <- rownames(deg)
data <- add_regulate(deg, log2FC_name = "logFC",
                     fdr_name = "P.Value",log2FC = 1, fdr = 0.05)
p<-ggvolcano(data, x = "log2FoldChange", y = "padj",y_lab = "-Log10Pvalue",
             x_lab = "Log2FC",
             fills = c("#4C8BC0","#b4b4d8","#C75C64"),
             colors = c("#4C8BC0","#b4b4d8","#C75C64"),
             label = "Gene", label_number = 20, output = FALSE)+
  mytheme

p
save_plot(p,filename = "DEG/LncRNA//VolcanoPlot",style = "g")
fwrite(DEG,file = "DEG/LncRNA/DEG_result.csv",row.names = T)

Dexp <- exper[rownames(deg[deg$change !="not",]),]
annocol <- data.frame(group=Group) 
ann_colors=list(group=c(Control='#8EBAD5',T2DR='#D57A66'))
rownames(annocol) <- colnames(Dexp)
annocol <- arrange(annocol,group)
heatmap_exp <- Dexp[,rownames(annocol)]
p <- pheatmap(heatmap_exp, show_colnames = F, show_rownames = F,
              scale = "row",
              cluster_cols = F,
              annotation_col = annocol,
              breaks = seq(-3, 3, length.out = 100),annotation_colors = ann_colors,
              color = colorRampPalette(c('#4C8BC0','white','#C75C64'))(100),labels_col = "x",
              fontsize = 10)
p
save_plot(p,filename = "DEG/LncRNA//DEG_heatmap",style = "xx")

"#9787D7",'#E15CB8',"#8EBAD5",'#D57A66',"#269846","#b4b4d8","#e94234","#4C8BC0","#C75C64"
library(clusterProfiler)
library(org.Hs.eg.db)
library(CyDataPro)
data <- fread("DEG/mRNA/DEG_result.csv") %>% filter(change !="not")

Erich_GO_KEGG(data$V1,out = "DEG/mRNA/diffGene",color = c("#269846","#e94234"),plotStyle = "bar",width = 8,height = 6,method = "GO")
Erich_GO_KEGG(data$V1,out = "DEG/mRNA/diffGene",color = c("#269846","#e94234"),plotStyle = "bar",width = 8,height = 6,method = "KEGG")
gene_df <- bitr(data$V1,
                fromType = "SYMBOL",
                toType = "ENTREZID",
                OrgDb = org.Hs.eg.db)
enrich <- enrichKEGG(gene = gene_df$ENTREZID,
                     organism = "hsa",
                     keyType = "kegg",
                     pvalueCutoff = 1,minGSSize = 30,
                     qvalueCutoff = 1,
                     use_internal_data = F)
enrich_kegg<-setReadable(enrich,OrgDb = org.Hs.eg.db,keyType = "ENTREZID")
result <- enrich_kegg@result
data <- result[c(2:4,6:17),]
data <- data %>% arrange(desc(pvalue),desc(Count))
data$Description <- factor(data$Description,levels = data$Description)
y <- str_wrap(data$Description,width = 60)
plot1 <- ggplot(data,aes(Count,Description,fill=-log10(pvalue)))+
  geom_bar(stat = "identity",width = 0.8,alpha=0.6)+
  labs(x="Gene Counts",y="",title = "KEGG enrichment")+
  geom_text(aes(x=0.03,label=Description),hjust=0,size=5,color="black",fontface="bold")+
  scale_fill_gradient(low = color[1],high = color[2])+
  theme_classic()+
  theme(axis.text.y = element_blank())+mytheme


write.table(result,file = paste0(out,"_Kegg.txt"),sep = "\t",quote = F,row.names = F)


ggsave(filename = paste0(out,"_Kegg_barP.pdf"),plot1,width = 8,height = 6)
ggsave(filename = paste0(out,"_Kegg_barP.png"),plot1,width = 8,height = 6,dpi = 300)





data <- fread("DEG/Methy/Methy_DEG_Gene.csv")
Erich_GO_KEGG(data$Gene,out = "DEG/Methy//MethyGene",color = c("gray","red"),plotStyle = "bar",width = 8,height = 6,method = "GO")
Erich_GO_KEGG(data$Gene,out = "DEG/Methy//MethyGene",color = c("gray","red"),plotStyle = "bar",width = 8,height = 6,method = "KEGG")

data <- fread("DEG/LncRNA/DEG_result.csv") %>% filter(change !="not")
Erich_GO_KEGG(data$V1,out = "DEG/LncRNA//diffGene",color = c("#9787D7","#E15CB8"),plotStyle = "bar",width = 8,height = 6,method = "GO")
Erich_GO_KEGG(data$V1,out = "DEG/LncRNA//diffGene",color = c("#9787D7","#E15CB8"),plotStyle = "bar",width = 8,height = 6,method = "KEGG")

data <- fread("DEG/miRNA//DEG_result.csv") %>% filter(change !="not")
miR <- fread("../../rawdata/Starbase/mirRNA_target.txt")
fitR <- miR %>% filter(miRNAname %in% data$V1)
gene <- unique(fitR$geneName)
Erich_GO_KEGG(gene,out = "DEG/miRNA///diffGene",color = c("#4C8BC0","#C75C64"),plotStyle = "bar",width = 8,height = 6,method = "GO")
Erich_GO_KEGG(gene,out = "DEG/miRNA///diffGene",color = c("#4C8BC0","#C75C64"),plotStyle = "bar",width = 8,height = 6,method = "KEGG")


fwrite(data.frame(miRNA = data$V1),file = "diff_miRNA.txt")
rm(list = ls())
Methy <- read.csv("DEG/Methy/Methy_DEG_Gene.csv")
Methy <- unique(Methy$Gene)
mRNA <- read.csv("DEG/mRNA/DEG_result.csv") %>% filter(change !="not")
library(ggvenn)
data <- list(Methylation=Methy,
             Diff_mRNA = mRNA$X)

p<-ggvenn(data,show_percentage = F,fill_color = c("#269846","#e94234"),
       stroke_color = "white",
       set_name_color = c("#269846","#e94234"),
       text_color = "white",
       text_size = 10)
p
dir.create("NetWork/Methy_mRNA",recursive = T)
save_plot(p,filename = "NetWork/Methy_mRNA/venn_plot",style = "g")
inner <- intersect(Methy,mRNA$X)
fwrite(data.frame(Gene=inner),file = "NetWork/Methy_mRNA/inner_Gene.csv")

Erich_GO_KEGG(inner,out = "NetWork/Methy_mRNA/diffGene",color = c("gray","red"),plotStyle = "point",width = 12,height = 6,method = "GO")
Erich_GO_KEGG(inner,out = "NetWork/Methy_mRNA/diffGene",color = c("gray","red"),plotStyle = "point",width = 10,height = 6,method = "KEGG")

DEG <- read.csv("DEG/miRNA/DEG_result.csv") %>% filter(change !="not")
miR <- fread("../../rawdata/Starbase/mirRNA_target.txt")
miRD <- miR %>% filter(miRNAname %in% DEG$X)
length(unique(miRD$miRNAname))
mRNA <- read.csv("DEG/mRNA/DEG_result.csv") %>% filter(change !="not")

mir_mRNA <- miRD %>% filter(geneName %in% mRNA$X)
length(unique(mir_mRNA$geneName))
mir_mRNA<- mir_mRNA %>% dplyr::select(miRNAname,geneName)
colnames(mir_mRNA)<- c("miRNA","mRNA")
dir.create("NetWork/miRNA_mRNA")
fwrite(mir_mRNA,file = "NetWork/miRNA_mRNA/miRNA_mRNA_paried.csv")

Erich_GO_KEGG(unique(mir_mRNA$mRNA),out = "NetWork/miRNA_mRNA/diffGene",color = c("#4C8BC0","#C75C64"),plotStyle = "point",width = 12,height = 6,method = "GO")
Erich_GO_KEGG(unique(mir_mRNA$mRNA),out = "NetWork/miRNA_mRNA/diffGene",color = c("#4C8BC0","#C75C64"),plotStyle = "point",width = 12,height = 6,method = "KEGG")

DEG <- read.csv("DEG/LncRNA//DEG_result.csv") %>% filter(change !="not")
Lnc <- fread("../../../WinSoftwares/Download/lnctard2.0/lnctard2.0.txt")
lncD <- Lnc %>% filter(Regulator %in% DEG$X)
length(unique(miRD$miRNAname))
mRNA <- read.csv("DEG/mRNA/DEG_result.csv") %>% filter(change !="not")

lnc_mRNA <- lncD %>% filter(Target %in% mRNA$X)
length(unique(lnc_mRNA$Target))
lnc_mRNA<- lnc_mRNA %>% dplyr::select(Regulator,Target)
colnames(lnc_mRNA)<- c("LncRNA","mRNA")
dir.create("NetWork/LncRNA_mRNA")
fwrite(lnc_mRNA,file = "NetWork/LncRNA_mRNA//LncRNA_mRNA_paried.csv")

Erich_GO_KEGG(unique(lnc_mRNA$mRNA),out = "NetWork/LncRNA_mRNA//diffGene",color = c("#9787D7","#E15CB8"),plotStyle = "point",width = 12,height = 6,method = "GO")
Erich_GO_KEGG(unique(lnc_mRNA$mRNA),out = "NetWork/LncRNA_mRNA//diffGene",color = c("#9787D7","#E15CB8"),plotStyle = "point",width = 12,height = 6,method = "KEGG")

library(data.table)
library(tidyverse)
mr <- fread("DEG/mRNA/DEG_result.csv") %>% filter(change !="not")
colnames(xx)
xx<- fread("data/GSE121820_Differential_DNA_methylation.csv") %>% filter(`p-value(region)` < 0.001)
mr <- intersect(xx$Gene,mr$V1)

me <- fread("DEG/Methy/Methy_DEG_Gene.csv")
mrna <- fread("DEG/mRNA/DEG_result.csv") %>% filter(change !="not")
mir <- fread("NetWork/miRNA_mRNA/miRNA_mRNA_paried.csv")

length(unique(fi$mRNA))
mi <- fread("NetWork/miRNA_mRNA/miRNA_mRNA_paried.csv")
lc <- fread("NetWork/LncRNA_mRNA/LncRNA_mRNA_paried.csv")
fi <- mi %>% filter(mRNA %in% mr)
  fwrite(fi,file = "test.csv")
fii <- fi %>% filter(mRNA%in% lc$mRNA)

library(ggvenn)
data <- list(Methylation=me$Gene,
             Diff_mRNA = mRNA$X,
             miRNA=mir$mRNA)

p<-ggvenn(data,show_percentage = F,fill_color = c("#269846","#e94234","#9787D7"),
          stroke_color = "white",
          set_name_color = c("#269846","#e94234","#9787D7"),
          text_color = "white",
          text_size = 10,stroke_size = 1)
p
dir.create("NetWork/Methy_mRNA_miRNA",recursive = T)
save_plot(p,filename = "NetWork/Methy_mRNA_miRNA//venn_plot",style = "g")
inn1 <- intersect(Methy,mRNA$X)
inn2 <- intersect(inn1,mir$mRNA)

Up <- mrna %>% filter(V1 %in% inn2)
load("data/GSE189005_mRNA.Rdata")
library(pheatmap)
Dexp <- exper[inn2,]
annocol <- data.frame(group=Group) 
ann_colors=list(group=c(Control='#8EBAD5',T2DR='#D57A66'))
rownames(annocol) <- colnames(Dexp)
annocol <- arrange(annocol,group)
heatmap_exp <- Dexp[,rownames(annocol)]
p <- pheatmap(heatmap_exp, show_colnames = F, show_rownames = T,
              scale = "row",
              cluster_cols = F,
              annotation_col = annocol,fontsize_row = 12,
              border_color = NA,
              breaks = seq(-2, 2, length.out = 100),annotation_colors = ann_colors,
              color = colorRampPalette(c('#00A087B2','white','#e94234'))(100),labels_col = "x",
              fontsize = 10)
p
save_plot(p,filename = "NetWork/Methy_mRNA_miRNA/innerGene_heatmap",style = "xx",width = 7,height = 8)
fwrite(data.frame(Gene=inn2),file = "NetWork/Methy_mRNA_miRNA/inner_Gene.csv")

gene_df <- bitr(inn2,
                fromType = "SYMBOL",
                toType = "ENTREZID",
                OrgDb = org.Hs.eg.db)
enrich_GO <- enrichGO(gene = gene_df$ENTREZID,
                      OrgDb = org.Hs.eg.db,
                      keyType = "ENTREZID",
                      ont = "ALL",
                      pvalueCutoff = 1,
                      qvalueCutoff = 1,
                      readable = F,
                      minGSSize = 10)
enrich_all<-setReadable(enrich_GO,OrgDb = org.Hs.eg.db,keyType = "ENTREZID")
barplot(enrich_all)
result <- enrich_all@result
GO <- result%>% top_n(n=-10,wt=pvalue)
y <- str_wrap(GO$Description,width = 40)
p<-ggplot(GO, aes(x = reorder(Description,-log10(pvalue)), -log10(pvalue),fill=ONTOLOGY)) + 
  geom_bar(stat = 'identity',alpha = 0.7) + 
  scale_x_discrete(labels=y)+
  coord_flip() + 
  labs(x = "",
       y="-log10(pvalue)")+
  theme_bw() + #去除背景色
  theme(panel.grid =element_blank(),
        axis.text = element_text(size=12,color="black",face = "bold"),
        axis.title = element_text(size = 12,face = "bold",color = "black"),
        axis.line = element_line(linewidth = 1.5),
        axis.ticks.length = unit(0.15,"cm"),
        axis.ticks = element_line(linewidth = 1),
        text = element_text(size = 10,colour = "black",face = "bold"))+
  theme(panel.border = element_rect(size = 0.6))+

  scale_fill_manual(values = c("#e94234","#FDBF6F","gray20"))#设置颜色
save_plot(p,filename = "NetWork/Methy_mRNA_miRNA/innerGene_GO",style = "g")
fwrite(result,file = "NetWork/Methy_mRNA_miRNA/innerGene_GO.csv")

enrich <- enrichKEGG(gene = gene_df$ENTREZID,
                     organism = "hsa",
                     keyType = "kegg",
                     pvalueCutoff = 1,
                     qvalueCutoff = 1,
                     use_internal_data = F)
enrich_kegg<-setReadable(enrich,OrgDb = org.Hs.eg.db,keyType = "ENTREZID")
result <- enrich_kegg@result
GO <- result%>% top_n(n=-10,wt=pvalue)
y <- str_wrap(GO$Description,width = 40)
p<-ggplot(GO, aes(x = reorder(Description,log10(pvalue)), -log10(pvalue),fill= -log10(pvalue))) + 
  geom_bar(stat = 'identity',alpha = 0.7) + 
  scale_x_discrete(labels=y)+
  coord_flip() + 
  labs(x = "",
       y="-log10(pvalue)")+
  theme_bw() + #去除背景色
  theme(panel.grid =element_blank(),
        axis.text = element_text(size=12,color="black",face = "bold"),
        axis.title = element_text(size = 12,face = "bold",color = "black"),
        axis.line = element_line(linewidth = 1.5),
        axis.ticks.length = unit(0.15,"cm"),
        axis.ticks = element_line(linewidth = 1),
        text = element_text(size = 10,colour = "black",face = "bold"))+
  theme(panel.border = element_rect(size = 0.6))+
  
  scale_fill_gradient2(high = "#E31A1C",low = "black",mid = "#FDC086",midpoint = 1.479054)#设置颜色
p

save_plot(p,filename = "NetWork/Methy_mRNA_miRNA/innerGene_KEGG",style = "g")
fwrite(result,file = "NetWork/Methy_mRNA_miRNA/innerGene_KEGG.csv")




library(scRNAstat) 
library(Seurat)
library(ggplot2)
library(clustree)
library(cowplot)
library(dplyr)
library(data.table)
library(stringr)
library(harmony)

fs = list.files('data/GSE204880_RAW//',pattern = '^GSM')
samples=str_split(fs,'_',simplify = T)[,1]

lapply(unique(samples),function(x){
  y=fs[grepl(x,fs)]
  folder=paste0("data/GSE204880_RAW/", str_split(y[1],'_',simplify = T)[,1])
  dir.create(folder,recursive = T)
  #为每个样本创建子文件夹
  file.rename(paste0("data/GSE204880_RAW/",y[1]),file.path(folder,"barcodes.tsv.gz"))
  #重命名文件，并移动到相应的子文件夹里
  file.rename(paste0("data/GSE204880_RAW/",y[2]),file.path(folder,"features.tsv.gz"))
  file.rename(paste0("data/GSE204880_RAW/",y[3]),file.path(folder,"matrix.mtx.gz"))
})
rm(list = ls())
gc()
dir='data/GSE204880_RAW//' 
samples=list.files( dir )
samples 
sceList = lapply(samples,function(pro){ 
   #pro=samples[1] 
  print(pro)  
  sce = CreateSeuratObject(counts =  Read10X(file.path(dir,pro )) ,
                           project = gsub('^GSM[0-9]*_','',
                                          gsub('filtered_feature_bc_matrix','',pro) )  ,# pro, #
                           min.cells = 5,
                           min.features = 200 )
  print(sce)
  
  cid = tail(colnames(sce)[order(sce$nFeature_RNA)],10000)
  sce$nFeature_RNA[cid]
  sce=sce[,colnames(sce) %in% cid]
  print(sce)
  return(sce)
})

sce.all=merge(x=sceList[[1]],
              y=sceList[ -1 ],
              add.cell.ids =  gsub('^GSM[0-9]*_','',
                                   gsub('filtered_feature_bc_matrix','',samples))    )


as.data.frame(sce.all@assays$RNA@counts[1:10, 1:2])
head(sce.all@meta.data, 10)
table(sce.all$orig.ident) 

table(sce.all@meta.data$orig.ident)
sce.all$group<-ifelse(grepl("GSM6199208|GSM6199209",sce.all$orig.ident),
                      "DK_C",ifelse(grepl("GSM6199210|GSM6199211",sce.all$orig.ident),"DK",
                                    ifelse(grepl("GSM6199212",sce.all$orig.ident),"DR","DR_C")))

table(sce.all@meta.data$group)
###### step2:QC质控 ######
dir.create("./1-QC")
setwd("./1-QC")
# 如果过滤的太狠，就需要去修改这个过滤代码
source('../scRNA_scripts/qc.R')
sce.all.filt = basic_qc(sce.all)
print(dim(sce.all))
print(dim(sce.all.filt))
setwd('../')
input_sce = sce.all
basic_qc <- function(input_sce){
  #计算线粒体基因比例
  mito_genes=rownames(input_sce)[grep("^MT-", rownames(input_sce),ignore.case = T)] 
  print(mito_genes) #可能是13个线粒体基因
  #input_sce=PercentageFeatureSet(input_sce, "^MT-", col.name = "percent_mito")
  input_sce=PercentageFeatureSet(input_sce, features = mito_genes, col.name = "percent_mito")
  fivenum(input_sce@meta.data$percent_mito)
  
  #计算核糖体基因比例
  ribo_genes=rownames(input_sce)[grep("^Rp[sl]", rownames(input_sce),ignore.case = T)]
  print(ribo_genes)
  input_sce=PercentageFeatureSet(input_sce,  features = ribo_genes, col.name = "percent_ribo")
  fivenum(input_sce@meta.data$percent_ribo)
  
  #计算红血细胞基因比例
  Hb_genes=rownames(input_sce)[grep("^Hb[^(p)]", rownames(input_sce),ignore.case = T)]
  print(Hb_genes)
  input_sce=PercentageFeatureSet(input_sce,  features = Hb_genes,col.name = "percent_hb")
  fivenum(input_sce@meta.data$percent_hb)
  
  #可视化细胞的上述比例情况
  feats <- c("nFeature_RNA", "nCount_RNA", "percent_mito", "percent_ribo", "percent_hb")
  feats <- c("nFeature_RNA", "nCount_RNA")
  p1=VlnPlot(input_sce, group.by = "orig.ident", features = feats, pt.size = 0, ncol = 2) + 
    NoLegend()
  p1 
  w=length(unique(input_sce$orig.ident))/3+5;w
  ggsave(filename="Vlnplot1.pdf",plot=p1,width = w,height = 5)
  
  feats <- c("percent_mito", "percent_ribo", "percent_hb")
  p2=VlnPlot(input_sce, group.by = "orig.ident", features = feats, pt.size = 0, ncol = 3, same.y.lims=T) + 
    scale_y_continuous(breaks=seq(0, 100, 5)) +
    NoLegend()
  p2 
  w=length(unique(input_sce$orig.ident))/2+5;w
  ggsave(filename="Vlnplot2.pdf",plot=p2,width = w,height = 5)
  
  p3=FeatureScatter(input_sce, "nCount_RNA", "nFeature_RNA", group.by = "orig.ident", pt.size = 0.5)
  ggsave(filename="Scatterplot.pdf",plot=p3)
  
  #根据上述指标，过滤低质量细胞/基因
  #过滤指标1:最少表达基因数的细胞&最少表达细胞数的基因
  # 一般来说，在CreateSeuratObject的时候已经是进行了这个过滤操作
  # 如果后期看到了自己的单细胞降维聚类分群结果很诡异，就可以回过头来看质量控制环节
  # 先走默认流程即可
  if(F){
    selected_c <- WhichCells(input_sce, expression = nFeature_RNA > 200)
    selected_f <- rownames(input_sce)[Matrix::rowSums(input_sce@assays$RNA@counts > 0 ) > 3]
    input_sce.filt <- subset(input_sce, features = selected_f, cells = selected_c)
    dim(input_sce) 
    dim(input_sce.filt) 
  }
  
  input_sce.filt =  input_sce
  # par(mar = c(4, 8, 2, 1))
  # 这里的C 这个矩阵，有一点大，可以考虑随抽样 
  C=subset(input_sce.filt,downsample=100)@assays$RNA@counts
  dim(C)
  C=Matrix::t(Matrix::t(C)/Matrix::colSums(C)) * 100
  
  most_expressed <- order(apply(C, 1, median), decreasing = T)[50:1]
  pdf("TOP50_most_expressed_gene.pdf",width=14)
  boxplot(as.matrix(Matrix::t(C[most_expressed, ])),
          cex = 0.1, las = 1, 
          xlab = "% total count per cell", 
          col = (scales::hue_pal())(50)[50:1], 
          horizontal = TRUE)
  dev.off()
  rm(C)
  
  #过滤指标2:线粒体/核糖体基因比例(根据上面的violin图)
  selected_mito <- WhichCells(input_sce.filt, expression = percent_mito < 30)
  selected_ribo <- WhichCells(input_sce.filt, expression = percent_ribo > 3)
  selected_hb <- WhichCells(input_sce.filt, expression = percent_hb < 1 )
  length(selected_hb)
  length(selected_ribo)
  length(selected_mito)
  
  input_sce.filt <- subset(input_sce.filt, cells = selected_mito)
  input_sce.filt <- subset(input_sce.filt, cells = selected_ribo)
  input_sce.filt <- subset(input_sce.filt, cells = selected_hb)
  dim(input_sce.filt)
  
  table(input_sce.filt$orig.ident) 
  
  #可视化过滤后的情况
  feats <- c("nFeature_RNA", "nCount_RNA")
  p1_filtered=VlnPlot(input_sce.filt, group.by = "orig.ident", features = feats, pt.size = 0, ncol = 2) + 
    NoLegend()
  w=length(unique(input_sce.filt$orig.ident))/3+5;w 
  p1
  ggsave(filename="Vlnplot1_filtered.pdf",plot=p1_filtered,width = w,height = 5)
  
  feats <- c("percent_mito", "percent_ribo", "percent_hb")
  p2_filtered=VlnPlot(input_sce.filt, group.by = "orig.ident", features = feats, pt.size = 0, ncol = 3) + 
    NoLegend()
  p2
  w=length(unique(input_sce.filt$orig.ident))/2+5;w 
  ggsave(filename="Vlnplot2_filtered.pdf",plot=p2_filtered,width = w,height = 5) 
  return(input_sce.filt) 
  
}

#harmony整合多个单细胞样品
input_sce <- input_sce.filt
print(dim(input_sce))
input_sce <- NormalizeData(input_sce, 
                           normalization.method = "LogNormalize",
                           scale.factor = 1e4) 
input_sce <- FindVariableFeatures(input_sce,)
input_sce <- ScaleData(input_sce)
input_sce <- RunPCA(input_sce, features = VariableFeatures(object = input_sce))
seuratObj <- RunHarmony(input_sce, "orig.ident")
names(seuratObj@reductions)
seuratObj <- RunUMAP(seuratObj,  dims = 1:15, 
                     reduction = "harmony")
 p = DimPlot(seuratObj,reduction = "umap",label=T ) 
p
 # ggsave(filename='umap-by-orig.ident-after-harmony',plot = p)
input_sce=seuratObj
input_sce <- FindNeighbors(input_sce, reduction = "harmony",
                           dims = 1:15) 
input_sce.all=input_sce

#设置不同的分辨率，观察分群效果(选择哪一个？)
for (res in c(0.01, 0.05, 0.1, 0.2, 0.3, 0.5,0.8,1)) {
  input_sce.all=FindClusters(input_sce.all, #graph.name = "CCA_snn", 
                             resolution = res, algorithm = 1)
}
colnames(input_sce.all@meta.data)
apply(input_sce.all@meta.data[,grep("RNA_snn",colnames(input_sce.all@meta.data))],2,table)

p1_dim=plot_grid(ncol = 3, DimPlot(input_sce.all, reduction = "umap", group.by = "RNA_snn_res.0.01") + 
                   ggtitle("louvain_0.01"), DimPlot(input_sce.all, reduction = "umap", group.by = "RNA_snn_res.0.1") + 
                   ggtitle("louvain_0.1"), DimPlot(input_sce.all, reduction = "umap", group.by = "RNA_snn_res.0.2") + 
                   ggtitle("louvain_0.2"))
ggsave(plot=p1_dim, filename="Dimplot_diff_resolution_low.pdf",width = 14)

p1_dim=plot_grid(ncol = 3, DimPlot(input_sce.all, reduction = "umap", group.by = "RNA_snn_res.0.8") + 
                   ggtitle("louvain_0.8"), DimPlot(input_sce.all, reduction = "umap", group.by = "RNA_snn_res.1") + 
                   ggtitle("louvain_1"), DimPlot(input_sce.all, reduction = "umap", group.by = "RNA_snn_res.0.3") + 
                   ggtitle("louvain_0.3"))
ggsave(plot=p1_dim, filename="Dimplot_diff_resolution_high.pdf",width = 18)

p2_tree=clustree(input_sce.all@meta.data, prefix = "RNA_snn_res.")
p2_tree
ggsave(plot=p2_tree, filename="Tree_diff_resolution.pdf")
table(input_sce.all@active.ident) 
saveRDS(input_sce.all, "sce.all_int.rds")
return(input_sce.all)

sel.clust = "RNA_snn_res.0.5"
sce.all.int <- SetIdent(input_sce.all, value = sel.clust)
table(sce.all.int@active.ident) 
last_markers <- c("Pdgfrb", "Cd34", "Itga8", # mesangial (MC)
                  "Pecam1", "Emcn", "Flt1", # glomerular capillary endothelial cells (GC-EC) 
                  "Nphs1", "Nphs2", # podocytes
                  "Cd79a", "Cd79b", "Ly6d", "Mzb1", #B cells 
                  "Bcl11b", "Cxcr6", "Fam189b", # T cells
                  "S100a9", "S100a8", "Retnlg", #neutrophils
                  "Cybb", "Coro1a", "Pld4", #macrophages
                  "Pdgfrb", "Myl9", "Fhl2", # retinal pericytes (RPC)
                  "Pecam1", "Cdh5", "Kdr" # retinal endothelial cells (EC-R) 
)

p_umap=DimPlot(sce.all.int, reduction = "umap",raster = F,label = T,repel = T) 
p_umap 
dir.create("Singlecell/annotation",recursive = T)
genes_to_check=str_to_upper(last_markers) 
DotPlot(sce.all.int ,features = unique(c(last_markers)))  + 
  coord_flip() + 
  theme(axis.text.x=element_text(angle=45,hjust = 1))

h=length( genes_to_check )/6+3;h
ggsave(paste('Singlecell/annotation/check_for_gene','.pdf'))

p_all_markers = DotPlot(sce.all.int , features =  unique(last_markers))  + 
  coord_flip() + 
  theme(axis.text.x=element_text(angle=45,hjust = 1)) 
p_all_markers+p_umap


colnames(sce.all.int@meta.data) 
table(sce.all.int$RNA_snn_res.0.1)
celltype=data.frame(ClusterID=0:18 ,
                    celltype= 0:18)
celltype[celltype$ClusterID %in% c( 0),2]='EC-R' 
celltype[celltype$ClusterID %in% c( 1),2]='GC-EC' 
celltype[celltype$ClusterID %in% c( 2 ),2]= 'GC-EC' # 'myeloids'   
celltype[celltype$ClusterID %in% c( 3  ),2]='EC-R' # 'myeloids'  
celltype[celltype$ClusterID %in% c( 4 ),2]='EC-R' # 'myeloids' 
celltype[celltype$ClusterID %in% c( 5 ),2]='GC-EC'
celltype[celltype$ClusterID %in% c( 6),2]='RPC'
celltype[celltype$ClusterID %in% c( 7),2]='macrophages'
celltype[celltype$ClusterID %in% c( 8),2]='EC-R'
celltype[celltype$ClusterID %in% c( 9),2]='neutrophils'
celltype[celltype$ClusterID %in% c( 10),2]='other'
celltype[celltype$ClusterID %in% c( 11),2]='GC-EC'
celltype[celltype$ClusterID %in% c( 12),2]='B cells'
celltype[celltype$ClusterID %in% c( 13),2]='GC-EC'
celltype[celltype$ClusterID %in% c( 14),2]='podocytes'
celltype[celltype$ClusterID %in% c( 15),2]='other'
celltype[celltype$ClusterID %in% c(16),2]='other' 
celltype[celltype$ClusterID %in% c(17),2]='other' 
celltype[celltype$ClusterID %in% c(18),2]='RPC' 

head(celltype)
table(celltype$celltype)
sce.all.int@meta.data$celltype = "NA"

for(i in 1:nrow(celltype)){
  sce.all.int@meta.data[which(sce.all.int@meta.data$RNA_snn_res.0.5 == celltype$ClusterID[i]),'celltype'] <- celltype$celltype[i]}
Idents(sce.all.int)=sce.all.int$celltype

sel.clust = "celltype"
sce.all.int <- SetIdent(sce.all.int, value = sel.clust)
table(sce.all.int@active.ident) 

p_all_markers=DotPlot(sce.all.int, features =  unique(last_markers),
                      assay='RNA' ,group.by = 'celltype' )  + coord_flip()+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5))
p_umap=DimPlot(sce.all.int, reduction = "umap", group.by = "celltype",label = T,label.box = T)
p_all_markers+p_umap
ggsave('Singlecell/annotation/markers_umap_by_celltype.pdf',width = 16,height = 8)

p_umap=DimPlot(sce.all.int, reduction = "umap", group.by = "celltype",label = T,label.box = T)+mytheme
p_umap
save_plot(p_umap,filename = "Singlecell/annotation/anno_umap",style = "g")
save(sce.all.int,file = "data/Single.Rdata")


rm(list = ls())
load("data/Single.Rdata")

gene <- c("Tfrc","Ap2m1","Ap2a1","Dab2","Ppp1cb")

DotPlot(sce.all.int , features =  unique(gene))  + 
  coord_flip() + 
  theme(axis.text.x=element_text(angle=45,hjust = 1)) 

DR <- subset(sce.all.int,group==c("DK"))
Idents(DR) <- DR$celltype

gene_cell_exp <- AverageExpression(DR,
                                   features = gene,
                                   group.by = 'celltype',
                                   slot = 'data') 
gene_cell_exp <- as.data.frame(gene_cell_exp$RNA)
library(ComplexHeatmap)
df <- data.frame(colnames(gene_cell_exp))
colnames(df) <- 'Cell_Type'
top_anno = HeatmapAnnotation(df = df,#细胞名/cluster
                             border = T,
                             show_annotation_name = F,
                             gp = gpar(col = 'black'),
                             col = list(Cell_Type = c('B cells'="#9ECABE",
                                                      'EC-R'="#F6F5B4",
                                                      'podocytes'="#2F528F",
                                                      "GC-EC"="#E3AD68",
                                                      "macrophages"="#ACD45E",
                                                      "neutrophils"="#D95F02",
                                                      "RPC"="#CAB2D6",
                                                      "other" = "gray")))#颜色设置

#数据标准化缩放一下
dir.create("Singlecell/Expression/T2DN/",recursive = T)
marker_exp <- t(scale(t(gene_cell_exp),scale = T,center = T)) %>% na.omit()
pdf("SingleCell/Expression/T2DN///KeyGene_agvExpression_heatmap.pdf",width = 8,height = 6)
Heatmap(marker_exp,
        cluster_rows = T,
        cluster_columns = T,
        show_column_names = F,
        show_row_names = T,
        column_title = NULL,
        heatmap_legend_param = list(
          title=' '),
        col = colorRampPalette(c("blue","white","red"))(100),
        border = 'gray',
        rect_gp = gpar(col = "gray", lwd = 1),
        row_names_gp = gpar(fontsize = 15),
        column_names_gp = gpar(fontsize = 20),
        top_annotation = top_anno)
dev.off()
p<-DotPlot(DR, features = gene,group.by = "celltype")+RotatedAxis()+
  scale_x_discrete("")+scale_y_discrete("")+mytheme
p
save_plot(p,filename = "Singlecell/Expression/T2DN//KeyGene_agvExpression_Dotplot",style = "g")

dir.create("SingleCell/Expression/T2DN///GeneExpressionPlot/",recursive = T)
i=gene[1]
for(i in gene){
  p2<- FeaturePlot(DR, features = i, min.cutoff = "q9",reduction = "umap")+
    scale_colour_gradient(low = "gray",high = "red")+mytheme+theme(title = element_text(size = 20))
  p2
  save_plot(p2, filename = paste0("Singlecell/Expression/T2DN/////GeneExpressionPlot//",i),style = "g")
}

table(sce.all.int$group)
DR <- subset(sce.all.int,group==c("DK_C"))
Idents(DR) <- DR$celltype

gene_cell_exp <- AverageExpression(DR,
                                   features = gene,
                                   group.by = 'celltype',
                                   slot = 'data') 
gene_cell_exp <- as.data.frame(gene_cell_exp$RNA)
# marker_exp <- t(scale(t(gene_cell_exp),scale = T,center = T)) %>% na.omit()
 gene_cell_exp <- gene_cell_exp %>% rownames_to_column("Gene")
Control <- melt(gene_cell_exp)
Control$Group <- "Control"
table(all$variable)
T2DR <- melt(gene_cell_exp)
T2DR$Group <-"T2DN"
all <- rbind(T2DR,Control)

for(i in gene){
  tmp <- all %>% filter(Gene==i)
  
  p<-ggplot(tmp,aes(x=variable,y=value,fill=Group))+
    geom_bar(stat = "identity",position = position_dodge(),width = 0.7,alpha=0.7)+
    scale_fill_manual(values = c("#33A02C","#E31A1C"))+
    labs(title = tmp$Gene[1],x="Cell",y="AverageExpression")+
    theme_classic()+theme(plot.title = element_text(hjust = 0.5),
                          axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1))+
    theme(axis.title = element_text(size = 15,face = "bold"),
          axis.line = element_line(linewidth = 1.5),
          axis.ticks.length = unit(0.15,"cm"),
          axis.ticks = element_line(linewidth = 1.5),
          axis.text = element_text(size = 18,face = "bold",colour = "black",),
          text = element_text(size = 18,colour = "black",face = "bold"))
  p
  save_plot(p,filename = paste0("Singlecell/Expression/T2DN////Case_Normal/",i,"barplot"),style = "g")
}
dir.create("Singlecell/Expression/T2DN/Case_Normal")  
