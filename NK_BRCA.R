rm(list = ls())
load("BC_rawList.Rdata")
sce <- merge(BC[[1]],BC[-1],add.cell.ids = names(BC),project='BC')
sce@meta.data$group <- str_split(sce@meta.data$orig.ident,"_",simplify = T)[,1]
sce <- PercentageFeatureSet(sce,pattern = '^MT-',col.name = 'percent.MT')
sce <- PercentageFeatureSet(sce,pattern = '^RP[SL]',col.name = 'percent.RP')
VlnPlot(sce,features = 'percent.MT',pt.size = 0)+NoLegend()
VlnPlot(sce,features = 'percent.RP',pt.size = 0)
p <-VlnPlot(sce,features = c("nFeature_RNA", "nCount_RNA", "percent.MT"),ncol = 3,pt.size = 0,group.by = "group")
p

p <-VlnPlot(sce,features = c("nFeature_RNA", "nCount_RNA", "percent.MT"),ncol = 3,pt.size = 0,group.by = "group")
save_plot(p,filename = "SingleCell/QC/MT_rate_villionPlot",style = "g",width = 15)
plot1 <- FeatureScatter(sce, feature1 = "nCount_RNA", feature2 = "percent.MT",group.by = "group")+theme(legend.position = "none")
plot2 <- FeatureScatter(sce, feature1 = "nCount_RNA", feature2 = "nFeature_RNA",group.by = "group")
p <- plot1 + plot2
p
save_plot(p,filename = "SingleCell/QC/Feature_Count_Cor_plot",style = "g",width = 10)
####过滤细胞
dim(sce)
sce <- subset(sce,subset = nCount_RNA>1000 & nFeature_RNA>100 & percent.MT<25 &percent.RP<30)
dim(sce)
# ###过滤基因
sce <- sce[rowSums(sce@assays$RNA@counts>0)>3,]
dim(sce)
sce <- NormalizeData(sce,normalization.method = "LogNormalize",scale.factor = 10000)#默认参数
###细胞周期评分，判断是否影响表达
#S.score 较高为S期，G2M.Score较高的为G2M期，都比较低的为G1期
s_feature <- cc.genes.updated.2019$s.genes
g2m_feature <- cc.genes$g2m.genes
sce <- CellCycleScoring(sce,
s.features = s_feature,
g2m.features = g2m_feature,
set.ident = T)
VlnPlot(sce,features = c('S.Score','G2M.Score'),group.by = 'orig.ident',pt.size = 0)
sce <- FindVariableFeatures(sce,selection.method = "vst",nfeatures = 2000)#默认参数
Top10 <- head(VariableFeatures(sce),10)
plot1 <- VariableFeaturePlot(sce)
plot2 <- LabelPoints(plot1,points = Top10,repel = T)
plot2
save_plot(plot2,filename = "SingleCell/QC/HVG_villion_plot",style = "g")
dim(sce)
allGenes <- rownames(sce)
sce <- ScaleData(sce,features = allGenes)
sce <- RunPCA(sce,features = VariableFeatures(object = sce))
p<-VizDimLoadings(sce, dims = 1:2)
p
save_plot(p,filename = "SingleCell/QC/top2PCA_gene",style = "g",width = 10)
p <-DimPlot(sce,reduction = "pca",group.by = "group")
p
save_plot(p,filename = "SingleCell/QC/PCA_DimPlot",style = "g")
DimHeatmap(sce, dims = 1:2)
sce <- JackStraw(sce, num.replicate = 200,dims = 40)#JackStraw方法
sce <- ScoreJackStraw(sce, dims = 1:40)
p <-JackStrawPlot(sce, dims = 1:40)
p
save_plot(p,filename = "SingleCell/QC/JackStrawPlot",style = "g")
p <-ElbowPlot(sce,ndims = 40)#看拐点
p
save_plot(p,filename = "SingleCell/QC/ElbowPlot",style = "g")
#前20个主成分
sce <- RunUMAP(sce,reduction = 'pca',dims = 20)
sce <- RunTSNE(sce,reduction = "pca",dims = 20)
p <-DimPlot(sce,reduction = "umap",group.by = "group")
p
save_plot(p,filename = "SingleCell/QC/pca20_umapPlot",style = "g")
p
p <-DimPlot(sce,reduction = "tsne",group.by = "group")
p
save_plot(p,filename = "SingleCell/QC/pca20_tsnePlot",style = "g")
sce <- FindNeighbors(sce,dims = 20)

sce_res <- sce
for (i in c(0.01, 0.05, 0.1, 0.15, 0.2, 0.3,0.4, 0.5,0.8,1)){
  sce_res <- FindClusters(sce_res,resolution = i)
}
p<-clustree(sce_res,prefix = 'RNA_snn_res.')
p
#resolution=0.8
sce <- FindClusters(sce,resolution = 0.5)
p<-DimPlot(sce,reduction = 'tsne',group.by = 'seurat_clusters',label = T)
p
p<-clustree(sce_res,prefix = 'RNA_snn_res.')
save_plot(p,filename = "SingleCell/QC/resolution_clustree_plot",style = "g",height = 10)
p<-DimPlot(sce,reduction = 'tsne',group.by = 'seurat_clusters',label = T)
p
save_plot(p,filename = "SingleCell/QC/Cluster_tsnePlot",style = "g")
p<-DimPlot(sce,reduction = 'umap',group.by = 'seurat_clusters',label = T)
p
save_plot(p,filename = "SingleCell/QC/Cluster_umapPlot",style = "g")
xx<-table(sce$singleR_label) %>% as.data.frame()
colnames(xx)<- c("Cell","Number")
fwrite(xx,file = "返修/细胞数量.csv")

mark2 <- c("CD3E","TRAC",#T-Cells
           "CD4","IL7R",#CD-4
           "CD8A","CD8B",#CD8
           "FOXP3","CTLA4",#Treg
           "GNLY","PRF1","XCL1",#NK-clees
           "LYZ",#myeloid lineage cells
           "FCGR3A","S100A8",#Mono
           "CD14","PLTP","MRCL","NLRP3","IL1B",#Macrophage cells
           "XCR1",#mDC1
           "FCER1A","CD1C",#mDC2
           "TNFRSF21",#pDC
           "CD79A",#B-Cell
           "EPCAM","CDH1","KRT18","KRT19")#tumour cells
gene_check2 <- unique(intersect(rownames(sce),makr2))
gene_check2 <- unique(intersect(rownames(sce),mark2))
p <- DotPlot(sce,features = gene_check2,assay = "RNA")+coord_flip()
p2 <- VlnPlot(sce,features = gene_check2[1:9])
p3 <- VlnPlot(sce,features = gene_check2[10:18])
p4 <- VlnPlot(sce,features = gene_check2[19:27])

p5 <- FeaturePlot(sce,features = gene_check2[1:9])
p6 <- FeaturePlot(sce,features = gene_check2[10:18])
p7 <- FeaturePlot(sce,features = gene_check2[19:27])

p <- DotPlot(sce,features = gene_check2,assay = "RNA")+coord_flip()
save_plot(p,filename = "SingleCell/annotation/markGene_dotplot",style = "g")
p <- FeaturePlot(sce,features = c( "GNLY","PRF1","XCL1"))
p
save_plot(p,filename = "SingleCell/annotation/NK_featurePlot",style = "g")

marker <- data.frame(cluster=0:20,cell=0:20)
marker[marker$cluster %in% c(-4),2] <- 'Other cells'
marker[marker$cluster != 4,2] <- 'Other cells'
sce@meta.data$cell_tyle <- sapply(sce@meta.data$seurat_clusters,function(x){marker[x,2]})

DimPlot(sce,reduction = "umap",group.by = "cell_type",label = T)
DimPlot(sce,reduction = "umap",group.by = "cell_tyle",label = T)
DimPlot(sce,reduction = "tsne",group.by = "cell_tyle",label = T)
p <-DimPlot(sce,reduction = "umap",group.by = "cell_tyle",label = T)
p <-DimPlot(sce,reduction = "umap",group.by = "cell_tyle",label = T)+theme_prism()
p
save_plot(p,filename = "SingleCell/annotation/NK_cell_umapPlot",style = "g")
p <-DimPlot(sce,reduction = "tsne",group.by = "cell_tyle",label = T)+theme_prism()
save_plot(p,filename = "SingleCell/annotation/NK_cell_tsnePlot",style = "g")


library(SingleR)
library(celldex)
load("../../rawdata/SingleR/HumanAnno_SingleR.Rdata")
anno <- SingleR(sce@assays$RNA@data,
                              ref = hpca.se,
                              labels = hpca.se$label.main,
                              clusters = sce@meta.data$seurat_clusters,
                              assay.type.test = "logcounts",
                              assay.type.ref = "logcounts",
                              de.method = "wilcox")
anno$labels[4] <- "NK cells"
sce@meta.data$singleR_label <- unlist(lapply(sce@meta.data$seurat_clusters, function(x){anno$labels[x]}))
p <-DimPlot(sce,reduction = "umap",group.by = "singleR_label",label = T)+theme_prism()
p
save_plot(p,filename = "SingleCell/annotation/NK_cell_umapPlot",style = "g")
save_plot(p,filename = "SingleCell/annotation/NK_cell_umapPlot",style = "g")
p
p <-DimPlot(sce,reduction = "tsne",group.by = "singleR_label",label = T)+theme_prism()
save_plot(p,filename = "SingleCell/annotation/NK_cell_tsnePlot",style = "g")
saveRDS(sce,"BC_cellAnno.rds")

diffmarks <- FindAllMarkers(sce,only.pos = T,min.pct = 0.25,logfc.threshold = 1)
View(diffmarks)
write.table(diffmarks,"SingleCell/FindAllmarks_result.txt",row.names = T,sep = "\t",quote = F)
NKdiff <- diffmarks %>% filter(cluster==4&p_val_adj < 0.05)
write.table(NKdiff,"SingleCell/NK_marks_result.txt",row.names = T,sep = "\t",quote = F)


###细胞通讯
rm(list = ls())
library(CellChat)
library(patchwork)
library(Seurat)
options(stringsAsFactors = FALSE)
setwd("D:/WORK/Project/乳腺癌NC细胞/")

sce.all.int<- readRDS("乳腺癌NC细胞/SingleCell/BC_cellAnno.rds")

data.input <- sce.all.int@assays$RNA@data
meta.data <- sce.all.int@meta.data
table(sce.all.int$singleR_label)
cellchat <- createCellChat(object=data.input,
                           meta = meta.data,
                           group.by='singleR_label')

cellchatDB <- CellChatDB.human
cellchat@DB <- cellchatDB
rm(sce.all.int)
gc()
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

cellchat <- updateCellChat(cellchat)

cellchat <- computeCommunProb(cellchat,population.size = F)
cellchat <- filterCommunication(cellchat,min.cells = 10)

df.net <- subsetCommunication(cellchat)
df.pathway <- subsetCommunication(cellchat,slot.name = 'netP')

cellchat <- computeCommunProbPathway(cellchat)

cellchat <- aggregateNet(cellchat)
dir.create("返修2/NK细胞通讯")
pdf("GSE161529/Cellchat/Interaction_Number_wight.pdf",onefile = F,width = 10,height = 6)
groupSize <- table(cellchat@idents) %>% as.numeric()
par(mfrow = c(1, 2), xpd = TRUE)
netVisual_circle(cellchat@net$count,vertex.weight = groupSize,
                 weight.scale = T,label.edge = F,
                 title.name = 'number of Interaction')
pdf("Cellchat/Interaction_Weight.pdf",width = 8,height = 6)
netVisual_circle(cellchat@net$weight,vertex.weight = groupSize,
                 weight.scale = T,label.edge = F,
                 title.name = 'Interaction Weight')
dev.copy2pdf()
file.rename("Rplot.pdf","返修2/NK细胞通讯//Interaction_number_Weight.pdf")

p1 <- netVisual_heatmap(cellchat,color.heatmap = "Reds",font.size = 15,font.size.title = 15)
p1
p2 <- netVisual_heatmap(cellchat,color.heatmap = "Reds",measure = "weight",font.size = 15,font.size.title = 15)
pdf("返修2/NK细胞通讯/Interaction_Number_Weight_heatmap.pdf",width = 15,height = 6)
p1+p2
dev.off()

mat <- cellchat@net$weight
par(mfrow = c(3,3),xpd=T)
pdf("返修2/NK细胞通讯//CellNetwork.pdf",width = 8,height = 6)
for (i in 1:nrow(mat)){
  mat2 <- matrix(0,nrow = nrow(mat),ncol = ncol(mat),
                 dimnames = dimnames(mat))
  mat2[i,] <- mat[i,]
  netVisual_circle(mat2,vertex.weight = groupSize,
                   weight.scale = T,edge.weight.max = max(mat),
                   title.name = rownames(mat)[i])
}
dev.off()


###NK基因功能富集
install.packages("genekitr")
library(clusterProfiler)
library(data.table)
library(CyDataPro)
library(org.Hs.eg.db)
library(tidyverse)
gene<- fread("返修/SingleCell/NK_Gene_FindAllmarker_result.csv",data.table = F)
dir.create("返修2/NK基因富集")
Erich_GO_KEGG(gene$gene,out = "返修2/NK基因富集/NK_Gene",
              top=10,method = "GO",plotStyle = "bar",width = 8,height = 6,color = c("#E69056","#325CAC"))
Erich_GO_KEGG(gene$gene,out = "返修2/NK基因富集/NK_Gene",
              top=10,method = "KEGG",plotStyle = "bar",width = 8,height = 6,color = c("#1F78B4","#E31A1C"))

edit(plotEnrich)
#####TCGA-BRAC#####
rm(list = ls())

library(data.table)
library(tidyverse)
library(tidymodels)
library(ggprism)
library(discrim)
library(ggDCA)
library(survival)
library(survminer)
setwd("D:/WORK/Project/P2/")

expre <- fread("../../rawdata/TCGA/TCGA-BRCA/exprmatx_without_normal.txt",data.table = F)
rownames(expre) <- expre[,1]
expre <- expre[,-1]
NK_gene <- read.delim("SingleCell/NK_marks_result.txt")
NK_expre <- expre[intersect(rownames(expre),NK_gene$gene),]
clin <- read.delim("../../rawdata/TCGA/TCGA-BRCA/Clin_process.txt",row.names = 1,check.names = F)
clin <- clin[clin$OS>0,]
clin$status <- ifelse(clin$status=="Alive",1,0)
NK_expre <- NK_expre[,intersect(colnames(NK_expre),rownames(clin))]
###批量单因素

cox <- apply(
  NK_expre,1,function(x){
    clin$genes <- ifelse(x > median(x),"High","Low")
    cox_genes <- coxph(Surv(OS,status)~genes,data = clin)
    beta <- coef(cox_genes)
    se <- sqrt(diag(vcov(cox_genes)))
    HR <- exp(beta)
    HRse <- HR * se
    cox_need <- round(cbind(coef=beta,
                            se=se,
                            z=beta/se,
                            p = 1-pchisq((beta/se)^2,1),
                            HR = HR, 
                            HRse = HRse,
                            HRz = (HR-1)/HRse,
                            HRp = 1-pchisq(((HR-1)/HRse)^2,1),
                            over_95 = exp(beta-qnorm(.975,0,1)*se),
                            upper_95 = exp(beta+qnorm(.975,0,1)*se)),3)
    return(cox_need["genesLow",])
  }
)
genes_cox <- t(cox)
write.table(genes_cox,"ML/NK_cox_result.txt",sep = "\t",quote = F,row.names = T)
gene_cox <- read.delim("ML/NK_cox_result.txt")
cox_gene <- genes_cox[genes_cox[,4] < 0.05,]

NK_expre <- t(NK_expre)
logd <- dectect_log(NK_expre)
identical(rownames(NK_expre),rownames(clin))
expre <- NK_expre[rownames(clin),]
# expre <- expre[,rownames(cox_gene)]
identical(rownames(expre),rownames(clin))
data <- cbind(clin[,2:3],expre)
data <- data[data$OS>0,]

data <- data[,-2]
data$status <- factor(data$status)
data$status <- ifelse(data$status==1,"Alive","Dead")
library(rstatix)
result <- melt(data) %>% group_by(variable) %>% 
  wilcox_test(value~status) %>% 
  adjust_pvalue(method = "fdr") %>% 
  add_significance("p.adj")

p_gene <- result %>% filter(p < 0.05) %>% select(variable)
data1 <- data.frame(status=data[,1],data[,p_gene$variable])
data1$status <- ifelse(data1$status=="Alive",1,0)
data1$status <- factor(data1$status)

data$status <- ifelse(data$status=="Alive",1,0)
data$status <- factor(data$status)
str(data)
set.seed(12345678)
split <- initial_split(data,0.8,strata = status)
train_data <- training(split)
test_data <- testing(split)




rec <- recipe(status~.,train_data)
  step_corr(all_predictors()) %>%
  step_center(all_predictors(), -all_outcomes()) %>%
  step_scale(all_predictors(), -all_outcomes()) %>%
  prep()


# XGBoost
xgb_mod <- boost_tree() %>%           
  set_engine("xgboost") %>%           
  set_mode("classification")   
# 决策树
dt_mod <- decision_tree() %>%           
  set_engine("rpart") %>%           
  set_mode("classification")
# 逻辑回归
logistic_mod <-          
  logistic_reg() %>%          
  set_engine('glm') %>% 
  set_mode("classification")

# nnet
nnet_mod <-          
  mlp() %>%          
  set_engine('nnet') %>%          
  set_mode('classification')
# 朴素贝叶斯        
naivebayes_mod <-          
  naive_Bayes() %>%          
  set_engine('naivebayes') %>% 
  set_mode("classification")
#KNN
kknn_mod <-          
  nearest_neighbor() %>%          
  set_engine('kknn') %>%          
  set_mode('classification')
# 随机森林
rf_mod <-          
  rand_forest() %>%          
  set_engine('randomForest') %>%          
  set_mode('classification')
# SVM        
svm_mod <-          
  svm_rbf() %>%          
  set_engine('kernlab') %>%          
  set_mode('classification')

# 设置重采样
set.seed(12345678)
folds <- bootstraps(train_data,20)
folds <- vfold_cv(train_data,v=10)
# 控制条件，保存预测值      
ctr <- control_resamples(save_pred = TRUE,verbose = T)


wf <- workflow_set(preproc=list(rec),          
                   models=list(xgb=xgb_mod,          
                               dt=dt_mod,          
                               log= logistic_mod,          
                               nb=naivebayes_mod,          
                               nnet=nnet_mod,          
                               knn=kknn_mod,          
                               rf=rf_mod,          
                               svm=svm_mod))
wf


# 模型拟合 
wf_res <- wf %>% 
  workflow_map("fit_resamples", #重采样拟合
               resamples=folds, # 重采样 
               control=ctr)
wf_res

a <- rank_results(wf_res,rank_metric = "roc_auc") %>% 
  filter(.metric=="roc_auc") %>% 
  select(model,mean)
a
write.table(a,file = "ML/model_meanAUC.txt",row.names = F,quote = F,sep = "\t")
p <-autoplot(wf_res,metric = "accuracy")+
  theme_bw()+mytheme
save_plot(p,filename = "ML/Multil_ML_ACC_boxplot",style = "g")

b<-collect_predictions(wf_res) %>%           
  group_by(model) %>%          
  roc_curve(status,.pred_0) %>%          
  ggplot(aes(x=1-specificity,y=sensitivity,color=model))+          
  geom_line(lwd=1)+
  geom_segment(x=0,y=0,xend = 1,yend = 1,color="red",size=2,linetype=6)+
  theme_bw()+mytheme+
  theme(legend.position = c(.98, .55),          
        legend.justification = c("right", "top"),
        legend.key.width = unit(2,"cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
save_plot(b,filename = "ML/Multil_ML_ROC",style = "g")


##最好模型运用于测试集


# svm_spce <- svm_rbf(cost = tune(),margin = tune(),rbf_sigma = tune()) %>% 
#   set_engine("kernlab") %>% 
#   set_mode("classification") 
# svm_grid <- grid_max_entropy(extract_parameter_set_dials(svm_spce),size = 20)
# svm_wf <- workflow() %>% 
#   add_model(svm_spce) %>% 
#   add_formula(status~.)
# set.seed(12345678)
# svm_folds <- vfold_cv(train_data,strata = status,v=10)
# set.seed(12345678)
# svm_res <- tune_grid(svm_wf,
#                      resamples = svm_folds,
#                      grid = svm_grid,
#                      control = control_grid(save_pred = T,verbose = T))

rm_space <- rand_forest(mtry = tune(),trees = tune()) %>% 
  set_engine('randomForest') %>%          
  set_mode('classification')
rm_wf <- workflow() %>% add_recipe(rec) %>% add_model(rm_space)
rm_folds <- vfold_cv(train_data,strata = status,v=10)
rm_grid <- grid_regular(mtry(range = c(5L,10L)),trees(range = c(500L,1000L)),levels = c(mtry=5,trees=6))
rm_res <- tune_grid(rm_wf,resamples = rm_folds,grid = rm_grid)
collect_metrics(rm_res)
autoplot(rm_res,metric = "roc_auc")+
  theme_classic()
best_auc <- select_best(rm_res,"roc_auc")
fit_rm <- finalize_workflow(rm_wf,best_auc) %>% fit(data = train_data)

pre <- predict(fit_rm,new_data = test_data,type = "prob")
pre %>% bind_cols(t=test_data$status) %>% roc_auc(t,.pred_0)

pred = predict(fit_rm,test_data) %>% 
  bind_cols(select(test_data,status)) %>% 
  bind_cols(predict(fit_rm,test_data,type = "prob"))

library(pROC)
roc1 <- roc(pred$status,pred$.pred_1)
p <-ggroc(roc1,legacy.axes = T,size=2,color="red")+
  geom_segment(x=0,y=0,xend = 1,yend = 1,size=2,linetype=4)+
  annotate("text",x=0.75,y=0.125,label=paste("AUC = ", round(roc1$auc,3)),size=6)+
  ggtitle("Rand_forest-ROC")+
  theme_bw()+mytheme
save_plot(p,filename = "ML/RF_ROC_test",style = "g")
###变量重要性
a <-fit_rm %>% extract_fit_parsnip() %>%
  vip(num_features=20,geom="col",mapping=aes(fill=Importance))+
  scale_fill_distiller(palette = "Spectral")+
  theme_classic()+mytheme
save_plot(a,filename = "ML/rf_top20_importance_plot",style = "g")

b <-fit_rm %>% extract_fit_parsnip() %>%
  vip(num_features=100)
improtance <- b$data

write.table(improtance,file = "ML/rf_Importance.txt",sep = "\t",quote = F,row.names = F)
                               
#RiskScore
rm(list = ls())
expre <- fread("../../rawdata/TCGA/TCGA-BRCA/exprmatx_without_normal.txt",data.table = F)
rownames(expre) <- expre[,1]
expre <- expre[,-1]
improtance <- read.delim("ML/rf_Importance.txt")
T20 <- improtance[1:20,]
exp1 <- expre[T20$Variable,]
clin <- read.delim("../../rawdata/TCGA/TCGA-BRCA/Clin_process.txt",row.names = 1,check.names = F)
clin <- clin[clin$OS>0,]
clin$status <- ifelse(clin$status=="Alive",1,0)
identical(colnames(exp1),rownames(clin))
exp2 <- exp1[intersect(colnames(exp1),rownames(clin))] %>% t()
clin2 <- clin[intersect(colnames(exp1),rownames(clin)),]
meta <- cbind(clin2[,2:3],exp2)
meta$status <- as.numeric(meta$status)
multiCox <- coxph(Surv(OS,status)~.,data = meta)
RS<- predict(multiCox,type = "risk",newdata = meta) %>% as.data.frame()
RS<- RS %>% rownames_to_column()
colnames(RS) <- c("sample","risk_score")
clin <- clin %>% rownames_to_column("sample")
RS_meta <-RS %>% inner_join(clin)
RS_meta$Group <- ifelse(RS_meta$risk_score > median(RS_meta$risk_score),"High","Low")
table(RS_meta$Group)
write.table(RS_meta,file = "ML/Clin_with_Risk.txt",row.names = F,quote = F,sep = "\t")
a <- read.delim("ML/Clin_with_Risk.txt")
median(a$risk_score)
fit_tcga <- survfit(Surv(OS,status)~Group,data=RS_meta)
p <- ggsurvplot(fit_tcga,data = RS_meta,
                censor.shape="|",censor.size=4,
                # conf.int = T,
                # con.int.style ="ribnbon",
                # con.in.alpha = 0.2,
                pval = TRUE,
                #palette = "lancet",
                surv.median.line = "hv",
                ggtheme = theme_prism(),
                legend ="top",
                legend.abs = c("c2","c1"),
                xlab = "OS_time(days)",
                ylab = "Survival probablity",
                title = "Survival Curves by Risk Group",
                break.x.by =2000,
                break.y.by = 0.2,
                risk.table = F,
                # risk.table.col = "strata",
                # risk.table.height = 0.2,
                risk.table.y.text = FALSE)
save_plot(p,filename = "ML/RS_Group_Surv_plot",style = "xx",width = 8,height = 6)

##GEO验证####
geo <- read.delim("../../rawdata/GEO/GSE20685_exp.txt")
geo_pha <- read.delim("../../rawdata/GEO/GSE20685_clin.txt")
geo_pha$Status<- ifelse(geo_pha$Status==1,0,1)
geo_pha <- geo_pha %>% filter(OS > 0)
geo <- t(geo)
geo_pha <- geo_pha[rownames(geo),]
identical(rownames(geo_pha),row.names(geo))
intersect(T20$Variable,colnames(geo))
geo_gene <- geo[,intersect(T20$Variable,colnames(geo))] %>% as.data.frame()
exp3 <- exp2[,intersect(T20$Variable,colnames(geo))]
clin2 <- clin[intersect(colnames(exp1),rownames(clin)),]
meta <- cbind(clin2[,2:3],exp3)
meta$status <- as.numeric(meta$status)
multiCox <- coxph(Surv(OS,status)~.,data = meta)
RS<- predict(multiCox,type = "risk",newdata = meta) %>% as.data.frame()
RS<- RS %>% rownames_to_column()
colnames(RS) <- c("sample","risk_score")
mevalue <- median(RS$risk_score)
geo_preditc<-predict(multiCox,type = "risk",newdata = geo_gene) %>% as.data.frame()
colnames(geo_preditc) <- "risk_score"
geo_preditc <- geo_preditc %>% rownames_to_column("sample")

geo_pha <- geo_pha %>% rownames_to_column("sample")
geo_meta <-geo_preditc %>% inner_join(geo_pha)
geo_meta$Group <- ifelse(geo_meta$risk_score > mevalue,"High","Low")
table(geo_meta$Group)
write.table(geo_meta,file = "GEO_exam/geo_RS_clin.txt",row.names = F,sep = "\t",quote = F)
##生存分析####
geo_meta$Os_time <- geo_meta$OS*365
test <- geo_meta %>% filter(Os_time >0)
test$Status <- ifelse(test$Status==1,0,1)
fit <- survfit(Surv(Os_time,Status)~Group,data = test)

p <- ggsurvplot(fit,data = test,
                censor.shape="|",censor.size=4,
                # conf.int = T,
                # con.int.style ="ribnbon",
                # con.in.alpha = 0.2,
                pval = TRUE,
                #palette = "lancet",
                surv.median.line = "hv",
                ggtheme = theme_prism(),
                legend ="top",
                legend.abs = c("c2","c1"),
                xlab = "OS_time(Years)",
                ylab = "Survival probablity",
                title = "Survival Curves by Risk Group",
                break.x.by =500,
                break.y.by = 0.2,
                risk.table = F,
                # risk.table.col = "strata",
                # risk.table.height = 0.2,
                risk.table.y.text = FALSE)
p
dir.create("GEO_exam")
save_plot(p,filename = "GEO_exam/GEO_exam_Surv_plot",style = "xx",width = 8,height = 6)

library(timeROC)
time_roc_res <- timeROC(T=TCGA$OS,delta = TCGA$status,
                        marker = TCGA$risk_score,
                        cause = 1,weighting = "marginal",
                        times = c(1 * 365, 3 * 365, 5 * 365),
                        ROC = T,iid = T)
time_roc_res$AUC
dat = data.frame(fpr = as.numeric(time_roc_res$FP),
                 tpr = as.numeric(time_roc_res$TP),
                 time = rep(as.factor(c(365*1,365*3,365*5)),each = nrow(time_roc_res$TP)))
p <-ggplot() +
  geom_line(data = dat,aes(x = fpr, y = tpr,color = time),size = 2) +
  scale_color_manual(name = NULL,values = c("#92C5DE", "#F4A582", "#66C2A5"),
                     labels = paste0("AUC of ",c(1,3,5),"-y survival: ",
                                     format(round(time_roc_res$AUC,2),nsmall = 2)))+
  geom_line(aes(x=c(0,1),y=c(0,1)),color = "grey")+
  theme_prism()+
  theme(plot.title = element_text(hjust = 0.5,size = 22),
        panel.grid = element_blank(),
        legend.background = element_rect(linetype = 1, size = 0.2, colour = "black"),
        legend.position = c(0.765,0.125),
        axis.text.x=element_text(hjust = 0.5,colour="black",size=16),
        axis.text.y=element_text(hjust = 0.5,colour="black",size=16),
        axis.title.x=element_text(size = 16),
        axis.title.y=element_text(size = 16),
        panel.border = element_blank(),axis.line = element_line(colour = "black",size=1),
        legend.text=element_text(colour="black",size=14),
        legend.title=element_text(colour="black", size=15))+
  scale_x_continuous(expand = c(0.005,0.005))+
  scale_y_continuous(expand = c(0.005,0.005))+
  labs(x = "1 - Specificity",
       y = "Sensitivity",title = "TCGA")+
  coord_fixed()
save_plot(p,filename = "ML/RS_ROC_plot",style = "g")
library(survivalROC)
sRocFuction=function(td=null,gene=null){
  par(mar= c(5,5,1,1),cex.lab=1.2,cex.axis= 1.2)
  sROC=survivalROC(Stime=td$OS, status=td$Status, marker = gene, predict.time =1, method="KM")
  plot(sROC$FP, sROC$TP, type="l", xlim=c(0,1), ylim=c(0,1),col="red", 
       xlab="False positive rate", ylab="True positive rate",
       lwd = 2, cex.main=1.3, cex.lab=1.5, cex.axis=1.2, font=1.2)
  abline(0,1)
  aucText=paste0("1 years"," (AUC=",sprintf("%.3f",sROC$AUC),")")
  sROC3=survivalROC(Stime=td$OS, status=td$Status, marker = gene, predict.time =3, method="KM")
  lines(sROC3$FP, sROC3$TP, type="l", xlim=c(0,1), ylim=c(0,1),col="green",lwd = 2)
  aucText3=paste0("3 years"," (AUC=",sprintf("%.3f",sROC3$AUC),")") 
  sROC5=survivalROC(Stime=td$OS, status=td$Status, marker = gene, predict.time =5, method="KM")
  lines(sROC5$FP, sROC5$TP, type="l", xlim=c(0,1), ylim=c(0,1),col="blue",lwd = 2)
  aucText5=paste0("5 years"," (AUC=",sprintf("%.3f",sROC5$AUC),")") 
  legend("bottomright", c(aucText,aucText3,aucText5),
         lwd=2,bty="n",col=c("red","green","blue"))
}
pdf(file = paste0("GEO_exam/GEO_exam_ROC_plot",".pdf"),width = 8,height = 6)
sRocFuction(td=test,gene = test$risk_score)
dev.off()
png(filename = paste0("GEO_exam/GEO_exam_ROC_plot",".png"),width = 8,height = 6,units = "in",res = 300)
sRocFuction(td=test,gene = test$risk_score)
dev.off()
p <-sRocFuction(td=test,gene = test$risk_score)
save_plot(p,filename = "GEO_exam/GEO_exam_ROC_plot",style = "xx")

TCGA <- read.delim("ML/Clin_with_Risk.txt")
TCGA$OS_time <- TCGA$OS/365
tcga_ROC_func <- edit(sRocFuction)
pdf(file = paste0("ML/RS_ROC_plot",".pdf"),width = 8,height = 6)
sRocFuction(td=test,gene = test$risk_score)
dev.off()
png(filename = paste0("ML/RS_ROC_plot",".png"),width = 8,height = 6,units = "in",res = 300)
sRocFuction(td=test,gene = test$risk_score)
tcga_ROC_func(td=TCGA,gene = TCGA$risk_score)

#######病理参数比较#####
rm(list = ls())

library(gmodels)
library(ggpie)
library(aplot)
install.packages("ggpie")
dia2 <- clin %>% select(status,Group) %>% as.factor()
dia2$status <- factor(dia2$status,levels = c("1","0"))
clin <- read.delim("ML/Clin_with_Risk.txt",row.names = 1)
data <- clin %>% mutate(Status = ifelse(status==1 ,"Alive","Dead"),)

data <- data %>% select(Group,Status,Stage,Gender,Age)

ChisqTset <-function(data){
  info <- c("Status","Stage","Gender","Age")
  result <- tibble(variable = character(), Pvalue = character(),method = character())
  for (i in 2:ncol(data)) {
    x = info[i-1]
    a = CrossTable(data[,1],data[,i],chisq = T,fisher = T)
    pvalue =a$chisq$p.value
    pcor = a$chisq.corr$p.value
    result = rbind(result,tibble(variable=x,Pvalue=pvalue,method="Chisq"))
  }
  return(result)
}

Chisq_result <- ChisqTset(data)
write.table(Chisq_result,"ML/RS_Clin_chisqTest_result.txt",row.names = F,quote = F,sep = "\t")
High <- data %>% filter(Group=="High")
Low <- data %>% filter(Group=="Low")
p1 <-ggdonut(High,group_key = "Status",count_type = "full",
             r0 = 0.5,r1=1,
             border_color = "white",label_type = "none",
             label_pos = "out",label_size = 8)+
  scale_fill_manual(values = c("#999999","#4D4D4D"))+labs(fill="Status")+
  theme(text = element_text(size = 15,colour = "black",face = "bold"),
        legend.position = "none")
p2 <- ggdonut(High,group_key = "Stage",count_type = "full",
              r0 = 0.5,r1=1,
              border_color = "white",label_type = "none",
              label_pos = "out",label_size = 8)+
  scale_fill_manual(values = c("#9ECAE1","#C6DBEF","#2171B5","#4292C6","#6BAED6"))+labs(fill="Stage")+
  theme(text = element_text(size = 15,colour = "black",face = "bold"),
        legend.position = "none")
p3 <- ggdonut(High,group_key = "Gender",count_type = "full",
              r0 = 0.5,r1=1,
              border_color = "white",label_type = "none",
              label_pos = "out",label_size = 8)+
  scale_fill_manual(values = c("#FF7F00","#FDBF6F"))+labs(fill="Gender")+
  theme(text = element_text(size = 15,colour = "black",face = "bold"),
        legend.position = "none")
p4 <- ggdonut(High,group_key = "Age",count_type = "full",
              r0 = 0.5,r1=1,
              border_color = "white",label_type = "none",
              label_pos = "out",label_size = 8)+
  scale_fill_manual(values = c("#33A02C","#B2DF8A"))+labs(fill="Age")+
  theme(text = element_text(size = 15,colour = "black",face = "bold"),
        legend.position = "none")

p5 <-ggdonut(Low,group_key = "Status",count_type = "full",
             r0 = 0.5,r1=1,
             border_color = "white",label_type = "none",
             label_pos = "out",label_size = 8)+
  scale_fill_manual(values = c("#999999","#4D4D4D"))+labs(fill="Status")+
  theme(text = element_text(size = 15,colour = "black",face = "bold"),
        legend.position = "bottom")
p6 <- ggdonut(Low,group_key = "Stage",count_type = "full",
              r0 = 0.5,r1=1,
              border_color = "white",label_type = "none",
              label_pos = "out",label_size = 8)+
  scale_fill_manual(values =  c("#9ECAE1","#C6DBEF","#4292C6","#2171B5","#6BAED6"))+labs(fill="Stage")+
  theme(text = element_text(size = 15,colour = "black",face = "bold"),
        legend.position = "bottom")
p7 <- ggdonut(Low,group_key = "Gender",count_type = "full",
              r0 = 0.5,r1=1,
              border_color = "white",label_type = "none",
              label_pos = "out",label_size = 8)+
  scale_fill_manual(values = c("#FF7F00","#FDBF6F"))+labs(fill="Gender")+
  theme(text = element_text(size = 15,colour = "black",face = "bold"),
        legend.position = "bottom")
p8 <- ggdonut(Low,group_key = "Age",count_type = "full",
              r0 = 0.5,r1=1,
              border_color = "white",label_type = "none",
              label_pos = "out",label_size = 8)+labs(fill="Age")+
  theme(text = element_text(size = 15,colour = "black",face = "bold"),
        legend.position = "bottom")+
  scale_fill_manual(values = c("#33A02C","#B2DF8A"))
apH <- p1 %>% insert_right(p2) %>% insert_right(p3) %>% insert_right(p4)
apL <- p5 %>% insert_right(p6) %>% insert_right(p7) %>% insert_right(p8)
ap <-apL+theme(legend.position = "bottom")

all <-ggarrange(p1,p2,p3,p4,p5,p6,p7,p8,ncol = 4,nrow = 2)
save_plot(xx,filename = "ML/RS_Clin_chisq_PiePlot",style = "g",width = 15,height = 8)
all
xx <-annotate_figure(all,top = text_grob("High",color = "black",face = "bold",size = 14,x=0,y=0.6,hjust = 1),
                     left = text_grob("Low",color = "black",face = "bold",size = 14,y=0.5))




####表达量比较###
rm(list = ls())
load("../../rawdata/DEG.Rdata")
exp <- read.delim("../../rawdata/TCGA/TCGA-BRCA/exprmatx_without_normal.txt",check.names = F)
Top20 <- read.delim("ML/rf_Importance.txt")[1:20,]

gene_exp <- exp[Top20$Variable,]
a <-sample_group %>% as.data.frame()
clin <- read.delim("ML/Clin_with_Risk.txt",row.names = 1)
group <- data.frame(Sample = rownames(clin),Group=clin$Group)
write.table(group,file = "../../rawdata/TCGA/TCGA-BRCA/Sample_Group.txt",row.names = F,sep = "\t",quote = F)
gene_exp <- gene_exp %>% rownames_to_column("Gene")
data<-reshape2::melt(gene_exp,id.vars = "Gene")
names(data) <- c("Gene","Sample","values")
data1 <- inner_join(data,group)
p <-ggplot(data1,aes(Gene,values,fill=Group))+
  geom_boxplot(color="black",alpha=0.6,outlier.shape = 21,outlier.size = 1.2)+
  scale_fill_manual(values = c("#d9352a","#4979b6"))+
  labs(x="Gene",y="Expression")+
  stat_compare_means(aes(group= Group),label = "p.signif",size=5)+theme_prism()+
  theme(axis.text.x = element_text(angle = 45))
  
dir.create("TCGA")
save_plot(p,filename = "TCGA/Top20Gene_expre_RsGroup_boxplot",style = "g",width = 12,height = 6)


edit(UCSC_TCGAdownload)
group <- group %>% column_to_rownames("Sample")
# library(ComplexHeatmap)
# anncol <- HeatmapAnnotation(Group = group$Group,
#                             col = list(RsGroup = c("Tumor"= "#4979b6", "normal"="#d9352a")))
# plot_data <- t(gene_exp)
# plot_data <- plot_data[,rownames(group)]
# identical(colnames(plot_data),rownames(group))
# data <- apply(plot_data,1,scale)
# rownames(data) <- colnames(plot_data)
# Heatmap(plot_data,na_col = "white",
#             show_column_names = F,
#             show_row_names = T,name = "GSVA Score",
#             column_order = c(colnames(plot_data)[c(grep("Tumor",group$Group),grep("normal",group$Group))]),
#             column_split = group$Group,
#             cluster_columns = F,column_title = NULL,
#             top_annotation = anncol)

##体细胞突变###
library(devtools)
library(data.table)
library(maftools)
library(tidyverse)
rm(list=ls())
maf <- fread("../../rawdata/TCGA/TCGA-BRCA/TCGA.BRCA.mutect.995c0111-d90b-4140-bee7-3845436c3b42.DR-10.0.somatic (1).maf.gz")
maf$Tumor_Sample_Barcode <- substr(maf$Tumor_Sample_Barcode,1,16)
clin <- read.delim("ML/Clin_with_Risk.txt",check.names = F,row.names = 1)
clin <- clin %>% rownames_to_column("Tumor_Sample_Barcode")
mafd <- read.maf(maf,clin)
vc_cols = c("#1F6CB6","red3","#70BC6B","#F4B9C5","#784194","#B81A7B","#A65628","#9E1F63","#005A32")  # RdYlBu
names(vc_cols) = c(
  'Missense_Mutation',
  'Multi_Hit',
  'Frame_Shift_Del',
  'Nonsense_Mutation',
  'Frame_Shift_Ins',
  'In_Frame_Ins',
  'Splice_Site',
  'In_Frame_Del',
  'Nonstop_Mutation'
)
gene <- read.delim("TCGA/Surv_key_gene.txt")

pdf("TCGA/Key_gene_mutation_oncoplot.pdf",width = 8,height = 6)
oncoplot(mafd,clinicalFeatures = "Group",sortByAnnotation = T,genes = Top20$Variable,
         showTumorSampleBarcodes = F,
         bgCol = "gray",colors = vc_cols,fontSize = .7)
dev.off()
pdf("TCGA/Top30_gene_mutation_oncoplot.pdf",width = 8,height = 6)
oncoplot(mafd,clinicalFeatures = "Group",sortByAnnotation = T,top = 30,
         showTumorSampleBarcodes = F,removeNonMutated = T,borderCol = NA,colors = vc_cols,fontSize = .7,sepwd_samples = 0)

dev.off()

###TMB###
tmb <- tmb(maf = mafd)
group <- clin %>% select("Tumor_Sample_Barcode","Group")
tmbGroping<- inner_join(tmb,group) %>% filter(total != 0)
data_p <- tmbGroping %>% select(total_perMB,Group) %>% melt()
p <-ggplot(data_p,aes(Group,value,fill=Group))+
  geom_violin(trim = T,aes(fill=Group),color="white",cex=1,alpha=0.4)+
  geom_boxplot(width=0.15,position = position_dodge(0.9),color="white")+
  # geom_signif(comparisons = list(c("High","Low")),step_increase = 0.1,size=1,test = "wilcox.test",textsize = 5)+
  #  geom_jitter()+
    stat_compare_means(aes(group=Group),size=5)+
  ylab("TMB")+xlab("Risk Group")+ylim(0.8,1.05)+
  theme_prism()+
  scale_fill_manual(values = c("#CA0020","#7FBC41"))
dir.create("返修/TCGA/TMB")
save_plot(p,filename = "返修/TCGA/TMB/RsGroup_TMB_boxplot",style = "g")

##预后
library(survival)
library(survminer)
group <- clin %>% dplyr::select("Tumor_Sample_Barcode","status","OS","Group")
survtmb <-inner_join(tmb,group) %>% filter(total != 0)
survtmb$TMB_Group <-if_else(survtmb$total_perMB >= median(survtmb$total_perMB), "High", "Low")

fit <- survfit(Surv(OS,status)~TMB_Group,data = survtmb)

p <- ggsurvplot(fit,data = survtmb,
                censor.shape="|",censor.size=4,
                # conf.int = T,
                # con.int.style ="ribnbon",
                # con.in.alpha = 0.2,
                pval = TRUE,
                # pval.coord = c(800,0.1),
                # pval.method = T,
                # pval.method.coord = c(1,0.1),
                #palette = "lancet",
                surv.median.line = "hv",
                ggtheme = theme_classic()+mytheme,
                legend =c(0.8,0.5),
                legend.abs = c("High","Low"),
                xlab = "OS_time(Days)",
                ylab = "Survival probablity",
                title = "Survival Curves by TMB Group",
                break.x.by =1000,
                break.y.by = 0.2,
                risk.table = F,
                # risk.table.col = "strata",
                # risk.table.height = 0.2,
                risk.table.y.text = FALSE)
p
save_plot(p ,filename = "返修/TCGA/TMB/Surv_TMBGroup_kvPlot",style = "xx")

survtmb$ConGroup <- ifelse(survtmb$Group =="High"& survtmb$total_perMB >= median(survtmb$total_perMB),"riskHigh-TMBHigh",
                           ifelse(survtmb$Group=="Low"&survtmb$total_perMB >= median(survtmb$total_perMB),"riskLow-TMBHigh",
                                  ifelse(survtmb$Group== "High"&survtmb$total_perMB < median(survtmb$total_perMB),"riskHigh-TMBLow","riskLow-TMBLow")))
table(survtmb$ConGroup)
fit <- survfit(Surv(OS,status)~ConGroup,data = survtmb)

p <- ggsurvplot(fit,data = survtmb,
                censor.shape="|",censor.size=4,
                # conf.int = T,
                # con.int.style ="ribnbon",
                # con.in.alpha = 0.2,
                pval = TRUE,
                # pval.coord = c(800,0.1),
                # pval.method = T,
                # pval.method.coord = c(1,0.1),
                #palette = "lancet",
                surv.median.line = "hv",
                ggtheme = theme_classic()+mytheme,
                legend =c(0.8,0.5),
                legend.abs = c("riskHigh-TMBHigh","riskHigh-TMBLow","riskLow-TMBHigh","riskLow-TMBLow"),
                xlab = "OS_time(Days)",
                ylab = "Survival probablity",
                title = "Survival Curves by riskScore-TMB Group",
                break.x.by =1000,
                break.y.by = 0.2,
                risk.table = F,
                # risk.table.col = "strata",
                # risk.table.height = 0.2,
                risk.table.y.text = FALSE)
p
save_plot(p,filename = "返修/TCGA/TMB/Surv_risk-TMB_plot",style = "xx")
write.table(tmbGroping,file = "返修/TCGA/TMB/TMB_result.txt",sep = "\t",quote = F,row.names = F)
##拷贝数变异
library(ChIPseeker)
library(GenomicRanges)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library(bigmemory)
library(tidyverse)
install.packages("bigmemory")
BiocManager::install("TxDb.Hsapiens.UCSC.hg38.knownGene")
setwd("D:/WORK/Project/P2")
rm(list = ls())
data <- read.delim(gzfile("../../rawdata/TCGA/TCGA-BRCA/TCGA-BRCA.ggene_cnv.tsv.gz"),check.names = F,row.names = 1)

# Key_gene <- read.delim("Ml/rf_Importance.txt")[1:20,] 
# cnv_dat <- data[data$`Gene Symbol` %in% Key_gene$Variable,]
# 
# 
# p <- data
# p$Chrom <- paste0("chr",p$Chrom)
# peak <- GRanges(sample = p[,1],
#                 Segment_Mean = p[,5],
#                 seqnames = Rle(p[,2]),
#                 ranges = IRanges(p[,3],p[,4]),
#                 strand = rep(c("*"),nrow(p)))
# txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
# gc()
# memory.limit()
# rm(data)
# rm(p)
# object.size(peak)
# options(future.globals.maxSize = 1000 * 1024^3)
# peakAnno <- annotatePeak(peak,tssRegion = c(-3000,3000),
#                          TxDb =TxDb.Hsapiens.UCSC.hg38.knownGene ,annoDb = "org.Hs.eg.db")

id <- read.delim("../../rawdata/TCGA/TCGA-BRCA/gencode.v22.annotation.gene.probeMap",row.names = 1)
geneinfo <- na.omit(id[,1:2])
int <- intersect(rownames(geneinfo),rownames(data))
geneinfo <- geneinfo[int,] %>% na.omit() #基因信息行数与表达矩阵匹配
geneinfo <- geneinfo[!duplicated(geneinfo$gene),] #去除重复
data1 <- data[rownames(geneinfo),]
which(is.na(geneinfo$gene))
tail(geneinfo)
class(geneinfo)
class(expr)
rownames(data1) <- geneinfo$gene 

gene_data <- read.delim("ML/rf_Importance.txt")[1:20,]
clin <- read.delim("ML/Clin_with_Risk.txt")
High <- clin %>% filter(Group == "High")
Low <- clin %>% filter(Group == "Low")
int <- intersect(gene_data$Variable,rownames(data1))
cnv_dat <- data1[int,]
cnv_dat[cnv_dat > 0.3] =0.4
cnv_dat[cnv_dat < -0.3] = -0.4
cnv_dat[cnv_dat < 0.3 & cnv_dat > -0.3 ]=0

cnv_dat[cnv_dat == 0.4] <- "Amplifications"
cnv_dat[cnv_dat == -0.4 ] <- "Deletions"
cnv_dat[cnv_dat == 0 ] <- "neutral"

cnv <- as.data.frame(t(cnv_dat))
cnv <- as.data.frame(t(apply(cnv,2,table)))
dat <- cnv[,1:2]
dat$Gene <- rownames(dat)
df <- reshape2::melt(dat)
library(ggalt)

p <-ggplot(dat,aes(y=Gene))+
  geom_point(data = df,aes(x=value,color=variable),size=3)+
  geom_dumbbell(aes(x=Deletions,xend =Amplifications),size=3,color = "#e3e2e1",
                colour_x = "green",colour_xend = "red",
                dot_guide = T,dot_guide_size = 0.25)+
  theme_bw()+coord_flip()+
  scale_color_manual(name = "",values = c("green","red"))+mytheme+
  theme(axis.text.x = element_text(angle = 45,hjust = 1))+
  xlab("CNV")
save_plot(p ,filename = "TCGA/ALL_CNV_plot",style = "g")
write.table(dat,file = "TCGA/All_CNV_result.txt",row.names = T,quote = F,sep = "\t")

cnv_High <- cnv_dat[,intersect(colnames(cnv_dat),High$sample)]
cnv_Low <- cnv_dat[,intersect(colnames(cnv_dat),Low$sample)]
dat_H <- as.data.frame(t(cnv_High))
cnv_H <- as.data.frame(t(apply(dat_H,2,table)))
cnv_H <- cnv_H[,1:2]
cnv_H$Gene <- rownames(cnv_H)
cnv_HD <- reshape2::melt(cnv_H)
cnv_HD$Group <- "High"

dat_L <- as.data.frame(t(cnv_Low))
cnv_L <- as.data.frame(t(apply(dat_L,2,table)))
cnv_L <- cnv_L[,1:2]
cnv_L$Gene <- rownames(cnv_L)
cnv_LD <- reshape2::melt(cnv_L)
cnv_LD$Group <- "High"

All <- rbind(cnv_H,cnv_L)
p <-ggplot(cnv_L,aes(y=Gene))+
  geom_point(data = cnv_LD,aes(x=value,color=variable),size=3)+
  geom_dumbbell(aes(x=Deletions,xend =Amplifications),size=3,color = "#e3e2e1",
                colour_x = "green",colour_xend = "red",
                dot_guide = T,dot_guide_size = 0.25)+
  theme_bw()+coord_flip()+
  scale_color_manual(name = "",values = c("green","red"))+mytheme+
  theme(axis.text.x = element_text(angle = 45,hjust = 1))+
  xlab("CNV")+
  ggtitle("Low Risk")
  
p <- save_plot(p ,filename = "TCGA/LowRisk_CNV_plot",style = "g")
write.table(cnv_L,file = "TCGA/LowRisk_CNV_result.txt",row.names = T,quote = F,sep = "\t")

p <-ggplot(cnv_H,aes(y=Gene))+
  geom_point(data = cnv_HD,aes(x=value,color=variable),size=3)+
  geom_dumbbell(aes(x=Deletions,xend =Amplifications),size=3,color = "#e3e2e1",
                colour_x = "green",colour_xend = "red",
                dot_guide = T,dot_guide_size = 0.25)+
  theme_bw()+coord_flip()+
  scale_color_manual(name = "",values = c("green","red"))+mytheme+
  theme(axis.text.x = element_text(angle = 45,hjust = 1))+
  xlab("CNV")+
  ggtitle("High Risk")
save_plot(p ,filename = "TCGA/HighRisk_CNV_plot",style = "g")
write.table(cnv_L,file = "TCGA/HighRisk_CNV_result.txt",row.names = T,quote = F,sep = "\t")

####甲基化##########

rm(list = ls())
setwd("D:/WORK/Project/P2/")
library(data.table)
library(tidyverse)
library(ggsignif)
library(ggprism)
library(ggpubr)
get_miss_rate <- function(x){sum(is.na(x))/length(x)*100}
data <- fread("../../rawdata/TCGA/TCGA-BRCA/TCGA-BRCA.methylation450.tsv.gz",data.table = F)
idmap <- fread("../../rawdata/TCGA/TCGA-BRCA/illuminaMethyl450_hg38_GDC",data.table = F)
key_gene <- read.delim("ML/rf_Importance.txt")[1:20,]

rownames(data)<- data[,1]
data <- data[,-1]
data <- data %>% rownames_to_column("id")
miss_prob <- data %>% apply(1, get_miss_rate) %>% as.data.frame() %>% rename(x=".")
prob <- miss_prob %>% filter(x<10) %>% rownames()
data_clean <- data[prob,] %>% rownames_to_column("id")

Pro_ID <- select(idmap,c("#id","gene")) 
Pro_ID <-Pro_ID[Pro_ID$gene!='.',] #去掉空值

Pro_ID$gene<- sapply(Pro_ID$gene, function(x){
  unlist(str_split(x,","))[1]
})
colnames(Pro_ID) <- c("id","Geneid")
exp_symbol <- merge(Pro_ID,data,by = "id") %>%
  select(-id) %>% group_by(Geneid) %>% 
  summarise_all(max) %>% filter(!Geneid =='')
Methy_data <- exp_symbol %>% column_to_rownames("Geneid")
s <- unique(exp_symbol$Geneid[duplicated((exp_symbol$Geneid))])
inner<- intersect(key_gene$Variable,rownames(Methy_data))
Methy_key <- Methy_data[inner,]
miss_prob <- Methy_key %>% apply(1, get_miss_rate) %>% as.data.frame() %>% rename(x=".")
prob <- miss_prob %>% filter(x<10) %>% rownames()
Methy_key <- Methy_key[prob,] 

cluster <- read.delim("ML/Clin_with_Risk.txt",check.names = F,row.names = 1)
inner <- intersect(rownames(cluster),colnames(Methy_key)) 
data_f <- Methy_key[,inner] %>% rownames_to_column("Geneid")
cluster1 <- cluster[inner,] %>% rownames_to_column("variable")
cluster1 <- cluster1 %>% select(variable,Group)
data1 <- melt(data_f)
data_p <- inner_join(data1,cluster1) %>% na.omit()
which(is.na(data_p$value))

p <-ggplot(data_p,aes(Geneid,value,fill=Group))+
  geom_boxplot(color="black",alpha=0.6,outlier.shape = 21,outlier.size = 1.2)+
  scale_fill_manual(values = c("#E41A1C","#4DAF4A"))+
  labs(x="Gene",y="Methylation Leve")+
  stat_compare_means(aes(group= Group),label = "p.signif",size=5)+theme_prism()+
  theme(axis.text.x = element_text(angle = 45))
save_plot(p,filename = "TCGA/Gene_Methy_boxplot",style = "g",width = 10)
p <-ggplot(data_p,aes(Group,value,fill=Group))+
  geom_violin(trim = F,aes(file=Group),color="white",cex=1,alpha=0.5)+
  geom_boxplot(width=0.1,position = position_dodge(0.9),color="white")+
  geom_signif(comparisons = list(c("High","Low")),step_increase = 0.1,size=1,test = "wilcox.test",textsize = 5)+
  #  geom_jitter()+
  #  stat_compare_means(aes(group=Cluster),label.x = 0.6,size=5,)+
  ylab("Methylation Level")+xlab("Risk Subtype")+ylim(0.8,1.05)+
  theme_prism()+
  scale_fill_manual(values = c("#B3CDE3","#DECBE4"))
save_plot(p,filename = "TCGA/Group_keyGene_Methy_boxplot",style = "g")  
write.table(data_p,file = "TCGA/KeyGene_MethyLevel.txt",row.names = F,quote = F,sep = "\t")

###独立预后###
rm(list = ls())
library(survival)
library(survminer)
library(plyr)
library(ezcox)
library(tidyverse)
install.packages("ezcox")
setwd("D:/WORK/Project/P2/")
clin <- read.delim("ML/Clin_with_Risk.txt",row.names = 1)
data <- cbind(clin[,3:ncol(clin)],riskScore=clin[,1])
data <- data[,-9]
data$Stage <-factor(data$Stage)
data <- data %>%dplyr::rename("riskGroup" = "Group")
variable.names<- colnames(data)[c(3:ncol(data))]
variable.names <- c("Gender","Age","riskGroup","riskScore")
result  <- ezcox(data,time = "OS",status = "status",covariates = variable.names,return_models = T)
dir.create("TCGA/IndependentFactor")
mod <- get_models(result)
a <- result$res %>% as.data.frame()
write.table(a,file = "TCGA/IndependentFactor/Single_cox_result.txt",row.names = F,sep = "\t",quote = F)
pdf("TCGA/IndependentFactor/Singlle1.pdf",width = 8,height = 6)
show_models(mod,model_names = paste0("Model ", 1:4))
dev.off()
png("TCGA/IndependentFactor/Singlle.png",width = 8,height = 6,units = "in",res = 300)
show_models(mod,model_names = paste0("Model ", 1:4))
dev.off()

td <- cbind(data[,1:2],data[,variable.names])
tdmultiCox=coxph(Surv(OS, status) ~ ., data = td)
x <-summary(tdmultiCox) 
pvalue=signif(as.matrix(x$coefficients)[,5],2)
HR=signif(as.matrix(x$coefficients)[,2],2)
low=signif(x$conf.int[,3],2)
high=signif(x$conf.int[,4],2)
multi_res=data.frame(p.value=pvalue,
                     HR=paste(HR," (",low,"-",high,")",sep=""),
                     stringsAsFactors = F
)
write.table(multi_res,file = "multi_cox_result.txt",row.names = F,quote = F,sep = "\t")
p <-ggforest(tdmultiCox,
         main = "Hazard ratio",
         cpositions = c(0.02,0.15, 0.35), #前三列的位置，第二列是样品数，设了个负值，相当于隐藏了
         fontsize = 1.3, #字体大小
         refLabel = "reference", 
         noDigits = 2)

save_plot(p,filename = "TCGA/IndependentFactor/Multi",style = "g",width = 12)

variable.names <- c("Gender","Age","riskGroup","riskScore","Stage")
td <- cbind(data[,1:2],data[,variable.names])
tdmultiCox=coxph(Surv(OS, status) ~ ., data = td)

Age <- coxph(Surv(OS, status) ~ Age, data = td)
Gender <- coxph(Surv(OS, status) ~ Gender, data = td)
Group <- coxph(Surv(OS, status) ~ riskGroup, data = td)
Score <- coxph(Surv(OS, status) ~ riskScore, data = td)
Stage <- coxph(Surv(OS, status) ~ Stage, data = td)
library(ggprism)
library(ggDCA)
library(regplot)

dca_cox <- dca(tdmultiCox,Age,Gender,Group,Score,Stage,model.names = c("Age+Gender+riskGroup+riskScore","Age","Gender","riskGroup","riskScore","Stage"))
p <-ggplot(dca_cox,
           linetype = F,lwd=1.2)+
  theme_classic()+  
  theme_prism(base_size =17)+
  theme(legend.position="top")+
  scale_x_continuous(
    limits = c(0.25, 1),
    guide = "prism_minor") +
  scale_y_continuous(
    limits = c(-0.01, 0.4),
    guide = "prism_minor")+
  scale_colour_prism(         
    palette = "candy_bright",
    name = "Cylinders")
save_plot(p,filename = "TCGA/IndependentFactor/DCA_plot",style = "g",width = 10,height = 6)



pdf(file = "TCGA/IndependentFactor/Nomogram_plot",width = 8,height = 8,onefile = FALSE)
regplot(tdmultiCox,observation = T,
        failtime = c(365,365*2),
        prfail = T,showP = T,
        droplines=F,
        rank = "sd",
        interval = "confidence",points = T,title = "Nomogram")
dev.off()
save_plot(p,filename = "TCGA/Lasso_cox/Nomogram_plot",style = "x",width = 8,height = 6)


library(nomogramFormula)

library(rms)
library(survivalROC)
td$Gender <- ifelse(td$Gender=="female",1,0)
td$Age <- ifelse(td$Age == ">45",1,0)
td$riskGroup <- ifelse(td$riskGroup=="High",1,0)
td$Stage <- ifelse(td$Stage == "I",1,ifelse(td$Stage == "II",2,
                                            ifelse(td$Stage == "III",3,
                                                   ifelse(td$Stage == "IV",4,0))))
td$OS_t <- td$OS/365
dd <- datadist(td)
options(datadist="dd")
str(td)
f <- cph(formula(Surv(OS_t, status) ~Age + riskGroup + Stage + riskScore), data = td, x = TRUE, y = TRUE,
         surv = TRUE, time.inc = 3)

surv <- Surv(f)
nomo <- nomogram(f)
result <- formula_rd(nomo)
result$formula
td$points <- points_cal(formula = result$formula,rd=td)
pdf("TCGA/IndependentFactor/NomoPoint_ROC.pdf",width = 8,height = 6)
sRocFuction(OS=td$OS_t,status = td$status,marker = td$points)
dev.off()
png("TCGA/IndependentFactor/NomoPoint_ROC.png",width = 8,height = 6,units = "in",res = 300)
sRocFuction(OS=td$OS_t,status = td$status,marker = td$points)
dev.off()
write.table(td,file = "TCGA/IndependentFactor/Nomopoint_result.txt",quote = F,sep = "\t",row.names = T)


###GSEA功能差异#######
rm(list = ls())
library(GSVA)
library(GSEABase)
library(rstatix)
library(reshape2)
exp <- read.delim("../../rawdata/TCGA/TCGA-BRCA/exprmatx_without_normal.txt",check.names = F)
Group <- read.delim("ML/Clin_with_Risk.txt",check.names = F,row.names = 1)
exp <- exp[,Group$sample]
gmt <- getGmt("../../rawdata/GSEA/c5.go.bp.v2023.1.Hs.symbols.gmt")
gsva_result <- gsva(as.matrix(exp),gmt,method="gsva",
                    kcdf = "Poisson",parallel.sz=15)

gsva <- t(gsva_result) %>% as.data.frame()
dir.create("TCGA/GSVA")
write.table(gsva,"TCGA/GSVA/gsva_BP_result.txt",sep = "\t",quote = F,row.names = T)
inner <- intersect(rownames(Group),rownames(gsva))
g_result <- gsva[inner,] 
clin_data <- Group[inner,]
identical(rownames(g_result),rownames(clin_data))
g_result$Group <- clin_data$Group
result <- melt(g_result) %>% group_by(variable) %>% 
  wilcox_test(value~Group) %>% 
  adjust_pvalue(method = "fdr") %>% 
  add_significance("p.adj")
write.table(result,"TCGA/GSVA/gsva_wilcox_test.txt",sep = "\t",quote = F)
sig_gsva <- result %>% filter(p.adj.signif != "ns")
top30<- arrange(sig_gsva,p.adj)[1:30,]
write.table(top30,"TCGA/ConsensusCluster/GSVA/top30_gsva_result.txt",quote = F,sep = "\t",row.names = F)
top30<- read.delim("TCGA/ConsensusCluster/GSVA/top30_gsva_result.txt")
top_gsva <- g_result[,top30$variable]

library(ComplexHeatmap)
anncol <- HeatmapAnnotation(Group = Group$Group,
                            col = list(Group = c("Low"= "#1A9641", "High"="#D7191C")),
                            annotation_legend_param = list(Group = list(nrow=1)))
plot_data <- t(top_gsva)
rownames(plot_data) <- str_remove(rownames(plot_data),pattern = "GOBP_")
p <-Heatmap(plot_data,na_col = "white",
            show_column_names = F,
            show_row_names = T,name = "GSVA Score",row_names_gp = grid::gpar(fontsize=11),
            column_order = c(colnames(plot_data)[c(grep("High",Group$Group),grep("Low",Group$Group))]),
            column_split = Group$Group,
            cluster_columns = F,column_title = NULL,
            top_annotation  = anncol,
            heatmap_legend_param = list(direction="horizontal"),
            column_names_rot = 90,row_names_max_width = unit(20,"cm"),
)
pdf("TCGA/GSVA/gsva_top30_heatmap.pdf",width = 15,height = 10)
draw(p,heatmap_legend_side="bottom",
     annotation_legend_side="bottom",merge_legend=T)
dev.off()
png("TCGA/GSVA/gsva_top30_heatmap.png",width = 15,height = 10,units = "in",res = 300)
draw(p,heatmap_legend_side="bottom",
     annotation_legend_side="bottom",merge_legend=T)
dev.off()

###免疫检查点#########
rm(list = ls())
IGG <- readRDS("../../../rawdata/IGG.Rds")
exp <- read.delim("../mRNA/TCGA/exprmatx_without_normal.txt",row.names = 1,check.names = F)
Rs_g <- read.delim("TCGA/Clin_process.txt",row.names = 1,check.names = F)
inner <- intersect(rownames(Rs_g),colnames(exp))
exp1 <- exp[,inner]
Rs_g1 <- Rs_g[inner,]
plot_data <- exp1[intersect(IGG,rownames(exp1)),]

data <- t(plot_data) %>% as.data.frame()
identical(rownames(data),rownames(Rs_g1))
data$Cluster<- Rs_g1$Cluster
Cluster_data <- melt(data)
result <- melt(data) %>% group_by(variable) %>% 
  wilcox_test(value~Cluster) %>% 
  adjust_pvalue(method = "fdr") %>% 
  add_significance("p.adj")
dir.create("TCGA/ConsensusCluster/IGGs")
write.table(result,file = "TCGA/ConsensusCluster/IGGs/Iggs_wilcox_test.txt",row.names = F,sep = "\t",quote = F)
sig_ICC <- result %>% filter(p.adj.signif != "ns")
plot_data <- Cluster_data %>% filter(variable %in% sig_ICC$variable)
p <- ggplot(plot_data,aes(variable,value,fill=Group))+
  geom_boxplot(color="black",alpha=0.6,outlier.shape = 21,outlier.size = 1.2)+
  scale_fill_manual(values = c("#D7191C","#1A9641"))+
  labs(x="IGGs",y="Expression level")+
  stat_compare_means(aes(group= Group),label = "p.signif",size=5)+
  theme_prism()+theme(axis.text.x = element_text(angle = 45))
p
dir.create("TCGA/Lasso_cox/IGG")

save_plot(p,filename = "TCGA/IGGs/Cluster_group_IGG_plot",width =10 ,height = 6,style = "g")

###TME ssGSEA###########
library(ggpubr)
library(lemon)
rm(list = ls())
tme <- read.delim("../P1/mRNA/TCGA/ConsensusCluster/TME_Score/Estamite_result.txt",row.names = 1)
RS_group <- read.delim("Ml/Clin_with_Risk.txt",row.names = 1)
tme1 <- tme[rownames(RS_group),]
tme1 <- tme1[,-5]

identical(rownames(tme1),rownames(RS_group))
tme1$RS_group <- RS_group$Group
plot_data <- melt(tme1)
dir.create("TCGA/TME")
p <-ggplot(plot_data,aes(RS_group,value,fill=RS_group))+
  geom_violin(trim = F,aes(fill=RS_group),color="white",cex=1,alpha=0.5)+
  geom_boxplot(width=0.1,position = position_dodge(0.9),color="gray")+
  # geom_signif(comparisons = list(c("High","Low")),step_increase = 0.1,size=1,test = "wilcox.test",textsize = 5)+
  facet_rep_wrap(.~variable,scales = 'free',repeat.tick.labels = 'left',ncol = 4)+
  stat_compare_means(comparisons = combn(unique(plot_data$RS_group), 2, simplify =FALSE),
                     method = 't.test',size=6)+
  theme_bw()+
  scale_fill_manual(values = c("#d9352a","#4979b6"))+ylab("TME Score")+xlab("")+
  theme(panel.grid=element_blank())+mytheme
  
write.table(tme1,"TCGA/TME/RS_TME_result.txt",quote = F,sep = "\t",row.names = T)
save_plot(p,filename = "TCGA/TME/RS_group_estimate_boxplot",style = "g",width = 12,height = 8)
##ssGSEA###
g_result <- read.delim("../P1/mRNA/TCGA/ConsensusCluster/GSVA/ssGSEA_result,txt",check.names = F,row.names = 1)
g_result <- select(g_result,-RSGroup,-Cluster)
clin_data <- read.delim("ML/Clin_with_Risk.txt",check.names = F,row.names = 1)
inner <- intersect(rownames(clin_data),rownames(g_result))
g_result <- g_result[inner,] 
clin_data <- clin_data[inner,]
identical(rownames(g_result),rownames(clin_data))
g_result$RS_group <- clin_data$Group
result <- melt(g_result) %>% group_by(variable) %>% 
  wilcox_test(value~RS_group) %>% 
  adjust_pvalue(method = "fdr") %>% 
  add_significance("p.adj")
write.table(result,"TCGA/TME/ssgsea_wilcoxTest.txt",row.names = F,sep = "\t",quote = F)
write.table(g_result,"TCGA/TME/ssGSEA_result.txt",row.names = T,quote = F,sep = "\t")
plot_data <- melt(g_result)
data_arrange <- plot_data %>% group_by(variable) %>% 
  summarise(de = median(value)) %>% 
  arrange(desc(de)) %>% pull(variable)
plot_data$variable <- factor(plot_data$variable,levels = unique(data_arrange))
p <- ggplot(plot_data,aes(variable,value,fill=RS_group))+
  geom_boxplot(color="black",alpha=0.6,outlier.shape = 21,outlier.size = 1.2)+
  scale_fill_manual(values = c("#E7298A","#1B9E77"))+
  labs(x="Cell type",y="ssGSEA Score")+
  stat_compare_means(aes(group= RS_group),label = "p.signif",size=5)+
  theme_prism()+
  theme(axis.text.x = element_text(angle = 45))
save_plot(p ,filename = "TCGA/TME/ssgsea_RSGroup_boxplot",style = "g",width = 10,height = 6)

###药物敏感性####
rm(list = ls())
Drug <- read.delim("../../rawdata/TCGA/TCGA-BRCA/BRCA_DrugPredictions.csv",sep = ",",row.names = 1)
Cluster <- read.delim("ML/Clin_with_Risk.txt",row.names = 1,check.names = F)
inner <- intersect(rownames(Cluster),rownames(Drug))
Drug1 <-Drug[inner,] 
identical(rownames(Drug1),rownames(Cluster))

data <- Drug1
identical(rownames(data),rownames(Cluster))
data$Group <- Cluster$Group
group_data <- melt(data)
result <- melt(data) %>% group_by(variable) %>% 
  wilcox_test(value~Group) %>% 
  adjust_pvalue(method = "fdr") %>% 
  add_significance("p.adj")
dir.create("TCGA/Drug")
write.table(result,"TCGA/Drug/Drug_wilcox_test.txt",row.names = F,sep = "\t",quote = F)
sig_Drug <- result %>% filter(p.adj.signif != "ns")
plot_data <- group_data %>% filter(variable %in% sig_Drug$variable)
p <- ggplot(plot_data,aes(variable,value,fill=Group))+
  geom_boxplot(width=1,position = position_dodge(0.9),color="black",outlier.shape = 21)+
  scale_fill_manual(values = c("#E7298A","#CAB2D6"))+
  labs(x="",y="Drug response level")+
  stat_compare_means(aes(group= Group),label = "p.signif",size=5)+
  theme_prism()+theme(axis.text.x = element_text(angle = 45))+
  geom_hline(yintercept = 0,lty=4,col="black",lwd=0.8)

p
save_plot(p ,filename = "TCGA//Drug/Drug_Sign_boxplot",style = "g",width = 20,height = 6)

plot_data <- data[,sig_Drug$variable] %>% t()
exp <- apply(plot_data, 1, scale)
rownames(exp) <- colnames(plot_data)
exp <- t(exp)
anncol <- HeatmapAnnotation(RsGroup = Cluster$Group,
                            col = list(RsGroup = c("Low"= "#1A9641", "High"="#D7191C")))
p <-Heatmap(exp,na_col = "white",
        show_column_names = F,
        # width = ncol(mat)*unit(5, "mm"), 
        height = nrow(exp)*unit(4, "mm"),
        show_row_names = T,name = "Drug response level",row_names_gp = grid::gpar(fontsize=11),
        column_order = c(colnames(exp)[c(grep("High",Cluster$Group),grep("Low",Cluster$Group))]),
        column_split = Cluster$Group,
        cluster_columns = F,column_title = NULL,
        top_annotation = anncol,)
pdf("TCGA/Drug/SigDrug_heatmap.pdf",width = 10,height = 15)
draw(p,heatmap_legend_side="right",
     annotation_legend_side="right",merge_legend=T)
dev.off()
png("TCGA/Drug/SigDrug_heatmap.png",width = 10,height = 15,units = "in",res = 300)
draw(p,heatmap_legend_side="right",
     annotation_legend_side="right",merge_legend=T)
dev.off()

data_arrange <- sig_Drug %>% group_by(variable) %>% 
  arrange(desc(p.adj)) %>% pull(variable)
sig_Drug <- sig_Drug %>% arrange(p.adj)
Top10 <- sig_Drug[1:10,]
plot_data <- group_data %>% filter(variable %in% Top10$variable)
p <- ggplot(plot_data,aes(variable,value,fill=Group))+
  geom_violin(trim = F,aes(fill=Group),color="white",alpha=0.5)+
  geom_boxplot(width=0.2,position = position_dodge(0.9),color="gray",outlier.shape = 21)+
  scale_fill_manual(values = c("#E7298A","#1A9641"))+
  labs(x="",y="Drug response level")+
  stat_compare_means(aes(group= Group),label = "p.signif",size=5)+
  theme_prism()+theme(axis.text.x = element_text(angle = 45,hjust  = 0.9,vjust = 0.9))+
  geom_hline(yintercept = 0,lty=4,col="black",lwd=0.8)
p
save_plot(p ,filename = "TCGA//Drug/Drug_Sign10_boxplot",style = "g",width = 8,height = 6)

###HPA下载#########
library(BiocStyle)
install.packages("HPAanalyze")
BiocManager::install("HPAanalyze")
library(HPAanalyze)
library(dplyr)
library(tibble)
library(readr)
library(tidyr)
dir.create("TCGA/HPA")
genes <- read.delim("ML/rf_Importance.txt")[1:20,]
genes <- genes$Variable
tissue="Breast"
filgene <- genes[18:20]
gene <-"AREG"
for (gene in filgene) {
  #获得HPA网站中该基因的xml文件
  hpa_target_gene<-try(hpaXmlGet(gene))
  if(class(hpa_target_gene)[1] == "try-error"){
    print(paste0(gene,"不在数据库！"))
    next
  }else{
    #将xml中组织染色的信息提取出来
    
    hpa_target_gene_fig_url<-hpaXmlTissueExpr(hpa_target_gene)
    hpa_target_gene_fig_url_1<-try(as.data.frame(hpa_target_gene_fig_url[[1]]))
    if(class(hpa_target_gene_fig_url_1) == "try-error"){
      print(paste0(gene," 下载不了"))
      next
    }else{
      tmp <- try(hpa_target_gene_fig_url_1[1:6,1:18])
      if(class(tmp) == "try-error"){
        print(paste0(gene," 下载不了"))
        next
      }else{
        #选择自己感兴趣的组织
        hpa_target_gene_fig_url_tissue<-hpa_target_gene_fig_url_1[hpa_target_gene_fig_url_1$tissueDescription2==tissue,]
        # hpa_target_gene_fig_url_tissue<-hpa_target_gene_fig_url_2[hpa_target_gene_fig_url_2$tissueDescription2==tissue,]
        
        # hpa_target_gene_fig_url_2<-as.data.frame(hpa_target_gene_fig_url[[2]])
        # hpa_target_gene_fig_url_2[1:6,1:18]#为该组织该基因单独建个文件夹储存
        picDir <- paste('TCGA/HPA/',gene, tissue,"IHC-2/", sep = "_")
        if (!dir.exists(picDir)) {
          dir.create(picDir)
        }
        
        
        for (i in 1:nrow(hpa_target_gene_fig_url_tissue)) {
          file_url<-hpa_target_gene_fig_url_tissue$imageUrl[i]
          file_dir<-paste(picDir,gene,tissue,hpa_target_gene_fig_url_tissue$patientId[i],hpa_target_gene_fig_url_tissue$tissueDescription1[i],hpa_target_gene_fig_url_tissue$tissueDescription2[i],".tiff",sep = "_")
          print(paste0("正在下载：",gene))
          download.file(url = file_url,destfile = file_dir,mode = "wb")
        }
      }
    }
  }
 
}

rm(sce)
#################返修############
setwd("D:/WORK/Project/乳腺癌NC细胞/")
gc()
library(Seurat)
sce<- readRDS("乳腺癌NC细胞/SingleCell/BC_cellAnno.rds")
table(sce$singleR_label)
diffmarks <- FindAllMarkers(sce,only.pos = T,min.pct = 0.25,logfc.threshold = 1,
                             group.by="singleR_label",test.use = "t",min.diff.pct = 0.1)
mark <- diffmarks %>% filter(cluster=="NK cells")
TT <- diffmarks %>% filter(cluster !="NK cells")
inn<-intersect(mark$gene,TT$gene)
mark<- mark %>% filter(!gene %in% inn)

mark<- fread("返修/NK_Gene_FindAllmarker_result.csv")


gene_cell_exp <- AverageExpression(sce,
                                   features = mark$gene,
                                   group.by = 'singleR_label',
                                   slot = 'data') 
gene_cell_exp <- as.data.frame(gene_cell_exp$RNA)
library(ComplexHeatmap)
df <- data.frame(colnames(gene_cell_exp))
colnames(df) <- 'Cell_Type'
top_anno = HeatmapAnnotation(df = df,#细胞名/cluster
                             border = T,
                             show_annotation_name = F,
                             col = list(Cell_Type = c('B_cell'="#9ECABE",
                                                      'CMP'="#F6F5B4",
                                                      'Endothelial_cells'="#2F528F",
                                                      "Epithelial_cells"="#E3AD68",
                                                      "Fibroblasts"="#ACD45E",
                                                      "Macrophage"="#D95F02",
                                                      "NK cells"="#CAB2D6",
                                                      "Smooth_muscle_cells" ="#E64B35FF",
                                                      "T_cells"="#00A087FF",
                                                      "Monocyte"="#8491B4FF")),
                             gp = gpar(col = 'black')
)#颜色设置


#数据标准化缩放一下
dir.create("返修")
marker_exp <- t(scale(t(gene_cell_exp),scale = T,center = T)) %>% na.omit()
pdf("./返修/SingleCell/NK_Gene_agvExpression_heatmap.pdf",width = 8,height = 10)
Heatmap(marker_exp,
        cluster_rows = T,
        cluster_columns = T,
        show_column_names = F,
        show_row_names = T,
        column_title = NULL,
        heatmap_legend_param = list(
          title=' '),
        col = colorRampPalette(c("#33A02C","white","#E31A1C"))(100),
        border = 'gray',
        rect_gp = gpar(col = "gray", lwd = 1),
        row_names_gp = gpar(fontsize = 15),
        column_names_gp = gpar(fontsize = 20),
        top_annotation = top_anno)
dev.off()
re<- diffmarks %>% filter(gene %in% mark$gene)
fwrite(mark,file = "返修/NK_Gene_FindAllmarker_result.csv")

Idents(sce) <- sce$singleR_label
i=mark$gene[1]
dir.create("返修/NK_Geneplot/UMAP")
my.colors<- colorRampPalette(c("lightblue","white","darkred"))(100)
for(i in mark$gene){
  p1<-FeaturePlot(sce, features = i,min.cutoff = "q10",reduction = "tsne",label = T,label.size = 4,col=my.colors)+mytheme
  p1
  CyDataPro::save_plot(p1, filename = paste0("返修/NK_Geneplot/Tsne/",i),style = "g")
  p2<- FeaturePlot(sce, features = i, min.cutoff = "q10",reduction = "umap",label = T,label.size = 4,col=my.colors)+mytheme

  CyDataPro::save_plot(p2, filename = paste0("返修/NK_Geneplot/UMAP//",i),style = "g")
}


#####TCGA-BRAC#####
rm(list = ls())

library(data.table)
library(tidyverse)
library(tidymodels)
library(ggprism)
library(discrim)
library(ggDCA)
library(survival)
library(survminer)
setwd("D:/WORK/Project/乳腺癌NC细胞/")

expre <- fread("../../rawdata/exprmatx_without_normal.txt",data.table = F)
rownames(expre) <- expre[,1]
expre <- expre[,-1]
NK_gene <- fread("返修/SingleCell/NK_Gene_FindAllmarker_result.csv")
NK_expre <- expre[intersect(rownames(expre),NK_gene$gene),]
clin <- read.delim("../../rawdata/Clin_process.txt",row.names = 1,check.names = F)
clin <- clin[clin$OS>0,]
clin <- clin[clin$Gender=="female",]
clin$status <- ifelse(clin$status=="Alive",1,0)
NK_expre <- NK_expre[,intersect(colnames(NK_expre),rownames(clin))]
clin<- clin[intersect(colnames(NK_expre),rownames(clin)),]
identical(colnames(NK_expre),rownames(clin))
###批量单因素
tmpdata<- data.frame(clin[,c(2,3)],t(NK_expre))
variable.names <- colnames(tmpdata)[3:ncol(tmpdata)]
OScoxSummary <- lapply(variable.names,function(x){
  print(x)
  OS <- na.omit(tmpdata[,c("status","OS",x)])
  res.cox = coxph(Surv(OS,status)~get(x),data=OS)
  coxSummary = summary(res.cox)
  coxSummary
  OScoxSummary <- data.frame(id = x,
                             N = nrow(OS),
                             HR = coxSummary$conf.int[,"exp(coef)"],
                             HR_95L = coxSummary$conf.int[,"lower .95"],
                             HR_95H = coxSummary$conf.int[,"upper .95"],
                             C_index = coxSummary$concordance["C"],
                             C_index_se = coxSummary$concordance["se(C)"],
                             pvalue = coxSummary$coefficients[,"Pr(>|z|)"])
  return(OScoxSummary)
})
UnivariateOScoxSummary <- do.call(rbind,OScoxSummary)
xx <- UnivariateOScoxSummary %>% filter(pvalue < 0.05)
dir.create("返修/Model")
fwrite(UnivariateOScoxSummary,file = "返修/Model/Gene_Univar_Cox_result.csv")
UnivariateOScoxSummary<- fread("返修/Model/Gene_Univar_Cox_result.csv",data.table = F)
xxx <- UnivariateOScoxSummary %>% filter(pvalue< 0.05)
xxx$change <- ifelse(xxx$HR > 1,"Risk factor","Protective factor")
xxx$Group <- ifelse(xxx$HR > 1,'darkred','lightblue')
p<-ggplot(xxx, aes(x =HR, y = reorder(id,HR),color=change)) +#映射
  geom_errorbarh(aes(xmin = HR_95L, xmax = HR_95H,),color=xxx$Group, height = 0,size=1.5) +#置信区间
  geom_point(alpha=1,size=6,shape=15) +#OR值，按照P值映射
  scale_color_manual(values = c('darkred','lightblue'))+#自定义颜色
  # xlim(min(or_ci), max(or_ci)) +#横坐标轴范围
  labs(x="Hazard ratio",y="") +#坐标轴标题
  geom_vline(xintercept = 1.0,color="darkblue",linetype=2,linewidth=1)+#添加无效线
  labs(color="")+
  xlim(0.7,1.5)+
  theme_bw()+#设置主题
  theme(axis.text = element_text(size = 12,color = "black"),
        panel.border = element_rect(linewidth = 1))+mytheme
p
save_plot(p,filename = "返修/Model/Gene_Univar_Sig_gene_ORplot",style = "g",width = 6)

devtools::install_local("D:/WinSoftwares/Download/Mime-main.zip")

library(Mime1)
library(tibble)
install.packages("rmeta")
data <- tmpdata %>% select(OS,status,xxx$id)
data<- data %>% rownames_to_column("ID")
colnames(data)[c(2,3)]<- c("OS.time","OS")
data$OS.time<- data$OS.time/365

##MEABRCA
tmp<- fread("brca_metabric/brca_metabric/data_clinical_patient.txt",data.table = F,skip=4)
tmp<- tmp[,c(1,14,15)]
colnames(tmp)<- c("ID","OS.time","OS")
tmp$OS<- ifelse(str_detect(tmp$OS,"LIVING"),0,1)
tmp$OS.time<- tmp$OS.time/12
tmp<- na.omit(tmp) %>% as.data.frame()
rownames(tmp)<- tmp$ID
tmp<- tmp[,-1]

tmp_exp<- fread("brca_metabric/brca_metabric/data_mrna_illumina_microarray.txt",data.table = F)
tmp_exp<- tmp_exp %>% filter(Hugo_Symbol %in% xxx$id)
tmp_exp<- tmp_exp %>% tibble::column_to_rownames("Hugo_Symbol")
tmp_exp<- tmp_exp[,-1] %>% t() %>% as.data.frame()
inn<- intersect(rownames(tmp),rownames(tmp_exp))

tmp<- tmp[inn,]
tmp_exp<- tmp_exp[inn,]
identical(rownames(tmp),rownames(tmp_exp))
meta<- cbind(tmp,tmp_exp)
meta<- meta %>% rownames_to_column("ID")

# geo <- read.delim("../../rawdata/GSE20685_exp.txt")
# geo_pha <- read.delim("../../rawdata/GSE20685_clin.txt")
# geo_pha$Status<- ifelse(geo_pha$Status==1,0,1)
# geo_pha <- geo_pha %>% filter(OS > 0)
# geo <- t(geo)
# geo_pha <- geo_pha[rownames(geo),]
# identical(rownames(geo_pha),row.names(geo))
# intersect(T20$Variable,colnames(geo))
# geo_gene <- geo[,intersect(xxx$id,colnames(geo))] %>% as.data.frame()
# exp3 <- exp2[,intersect(xxx$id,colnames(geo))]
# clin2 <- geo_pha[intersect(rownames(geo_gene),rownames(geo_pha)),]
# meta <- cbind(clin2[,c(4,3)],geo_gene)
# meta$Status <- as.numeric(meta$Status)
# meta<- meta %>% rownames_to_column("ID")
# colnames(meta)[c(2,3)]<- c("OS.time","OS")
# str(meta)

list_train_vali_Data<- list()
list_train_vali_Data$TCGA<- data
list_train_vali_Data$METABRIC<- meta
rm(expre,clin)
gc()
i="RSF"
ciall<- data.frame()
TCGA<- data[,1:3]
GEO<- meta[,1:3]
for(i in c("RSF", "Enet", "StepCox", "CoxBoost", "plsRcox", "superpc", "GBM", "survivalsvm", "Ridge", "Lasso")){
  res <- Mime1::ML.Dev.Prog.Sig(train_data = list_train_vali_Data$TCGA,
                                 list_train_vali_Data = list_train_vali_Data,
                                 unicox.filter.for.candi = F,
                                 # unicox_p_cutoff = 0.05,
                                 candidate_genes = xxx$id,
                                 mode = 'single', ## 'all', 'single', and 'double'
                                 nodesize =5,
                                 single_ml = i,
                                 seed = 123)
  cip<- res$Cindex.res
  
  ciall<- rbind(ciall,cip)
  TCGA<- cbind(TCGA,i=res[[3]][[1]][[1]][,4])
  colnames(TCGA)[ncol(TCGA)]=i
  GEO <- cbind(GEO,i=res[[3]][[1]][[2]][,4])
  colnames(GEO)[ncol(GEO)]=i
  }

GEO<- res$riskscore$RSF$METABRIC %>% as.data.frame()
library(pROC)
mm<- c("RSF", "Enet", "StepCox", "CoxBoost", "plsRcox", "superpc", "GBM", "survivalsvm", "Ridge", "Lasso")
roc1<- roc(TCGA$OS,as.numeric(TCGA$RSF),levels=c(1, 0), direction=">",smooth=T)
roc2<- roc(TCGA$OS,as.numeric(TCGA$Enet),levels=c(1, 0), direction="<",smooth=T)
roc3<- roc(TCGA$OS,as.numeric(TCGA$StepCox),levels=c(1, 0), direction="<",smooth=T)
roc4<- roc(TCGA$OS,as.numeric(TCGA$CoxBoost),levels=c(1, 0), direction="<",smooth=T)
roc5<- roc(TCGA$OS,as.numeric(TCGA$plsRcox),levels=c(1, 0), direction="<",smooth=T)
roc6<- roc(TCGA$OS,as.numeric(TCGA$superpc),levels=c(1, 0), direction="<",smooth=T)
roc7<- roc(TCGA$OS,as.numeric(TCGA$GBM),levels=c(1, 0), direction=">",smooth=T)
roc8<- roc(TCGA$OS,as.numeric(TCGA$survivalsvm),levels=c(1, 0), direction=">",smooth=T)
roc9<- roc(TCGA$OS,as.numeric(TCGA$Ridge),levels=c(1, 0), direction="<",smooth=T)
roc10<- roc(TCGA$OS,as.numeric(TCGA$Lasso),levels=c(1, 0), direction="<",smooth=T)
i=1
TCGA_AUC<- data.frame()
for(i in 1:10){
  auc <- round(auc(get(paste0("roc",i))),3)
  tmp<- data.frame(model=mm[i],auc=as.numeric(auc[1]))
  TCGA_AUC<- rbind(TCGA_AUC,tmp)
}
TCGA_AUC$data<- "TCGA"
roc1<- roc(GEO$OS,as.numeric(GEO$RS),levels=c(1, 0), direction=">",smooth=T)
roc2<- roc(GEO$OS,as.numeric(GEO$Enet),levels=c(1, 0), direction=">",smooth=T)
roc3<- roc(GEO$OS,as.numeric(GEO$StepCox),levels=c(1, 0), direction=">",smooth=T)
roc4<- roc(GEO$OS,as.numeric(GEO$CoxBoost),levels=c(1, 0), direction=">",smooth=T)
roc5<- roc(GEO$OS,as.numeric(GEO$plsRcox),levels=c(1, 0), direction=">",smooth=T)
roc6<- roc(GEO$OS,as.numeric(GEO$superpc),levels=c(1, 0), direction=">",smooth=T)
roc7<- roc(GEO$OS,as.numeric(GEO$GBM),levels=c(1, 0), direction=">",smooth=T)
roc8<- roc(GEO$OS,as.numeric(GEO$survivalsvm),levels=c(1, 0), direction="<",smooth=T)
roc9<- roc(GEO$OS,as.numeric(GEO$Ridge),levels=c(1, 0), direction=">",smooth=T)
roc10<- roc(GEO$OS,as.numeric(GEO$Lasso),levels=c(1, 0), direction=">",smooth=T)
i=1
GEO_AUC<- data.frame()
for(i in 1:10){
  auc <- round(auc(get(paste0("roc",i))),3)
  tmp<- data.frame(model=mm[i],auc=as.numeric(auc[1]))
  GEO_AUC<- rbind(GEO_AUC,tmp)
}
GEO_AUC$data<- "GSE20685"
AUC<- rbind(TCGA_AUC,GEO_AUC)

AUC<- fread("返修/Model/Model_10_AUC.csv")
AUC$data<- factor(AUC$data,levels = c("TCGA","GSE20685"))
  p<-ggplot(AUC,aes(x=reorder(model,-auc),y=auc,fill=data))+
  geom_bar(stat = "identity",position = position_dodge2())+
  labs(y="AUC",x="ML models",fill="")+
  scale_fill_manual(values = c("darkblue","lightblue"))+
  theme_classic()+mytheme+theme(axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1))
save_plot(p,filename = "返修/Model/Model_10_AUC",style = "g",width = 10)
fwrite(AUC,file = "返修/Model/Model_10_AUC.csv")

ciall$ID<- factor(ciall$ID,levels = c("TCGA","GSE20685"))
p<-ggplot(ciall,aes(x=reorder(Model,-Cindex),y=Cindex,fill=ID))+
  geom_bar(stat = "identity",position = position_dodge2(),alpha=0.7)+
  labs(y="Cindex",x="ML models",fill="")+
  scale_fill_manual(values = c("#D53E4F","#F46D43"))+
  theme_classic()+mytheme+theme(axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1))
p
save_plot(p,filename = "返修/Model/Model_10_Cindex",style = "g",width = 10)
fwrite(ciall,file = "返修/Model/Model_10_ciall.csv")






library(survival)
library(survminer)
res.cut <- surv_cutpoint(TCGA,time = "OS.time",
                         event = "OS",
                         variables = names(TCGA)[4:ncol(TCGA)],
                         minprop = 0.49)
res_cat<- surv_categorize(res.cut)
mysurv<- Surv(res_cat$OS.time,res_cat$OS)


Gene <- data.frame(Gene=vector(),Pval=vector())
i="RSF"
for(i in colnames(res_cat)[3:ncol(res_cat)]){
  group <- res_cat[,i]
  surv_dat <- data.frame(group=group)
  fit <- survfit(mysurv~group)
  
  group <-factor(group,levels = c("low","high"))
  data.survdiff <- survdiff(mysurv~group)
  p.val <- 1-pchisq(data.survdiff$chisq,length(data.survdiff$n)-1)
  HR <- (data.survdiff$obs[2]/data.survdiff$exp[2])/(data.survdiff$obs[1]/data.survdiff$exp[1])
  up95 <- exp(log(HR))+qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1])
  low95 <- exp(log(HR))-qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1])
  xx <- c(i,p.val)
  HR <- paste("Hazard Ratio = ",round(HR,2),sep = "")
  CI <- paste("95% CI: ",paste(round(low95,2),round(up95,2),sep = "-"),sep = "")
  svsort <- TCGA[order(TCGA[,i]),]
  cutoff <- paste("cutoff = ",round(svsort[fit$n[2],i],2),sep = "")
  
  
  p <- ggsurvplot(fit,data = surv_dat,
                  conf.int = F,
                  censor=F,palette = c("darkblue","lightblue"),
                   legend.title = i,size=2,
                  font.x=15,font.y=15,font.title=15,
                  # legend.labs = c(paste0("High ","(",fit$n[1],")"),
                  #                 paste0("Low ","(",fit$n[2],")")),
                  font.lenend=20,
                  legend ="top",
                  legend.abs = c("c2","c1"),
                  xlab = "OS_time(Years)",
                  ylab = "Survival probablity",
                  title = "Survival Curves by Risk Group",
                  break.x.by =5,
                  break.y.by = 0.2,ggtheme = theme_classic()+mytheme,
                  pval = paste(pval=ifelse(p.val < 0.001,"p < 0.001",paste("p = ",round(p.val,3),sep = "")),
                               HR,CI,sep = "\n")
  )

  Gene <- rbind(Gene,xx)
   if(p.val < 0.05){

     save_plot(p,filename = paste0("返修/Model/model_10_KMplot/",i,"_diff"),style = "x")
   }else{
     save_plot(p,filename = paste0("返修/Model/model_10_KMplot/",i,"_Nodiff"),style = "x")
   }
}
dir.create("返修/Model/model_10_KMplot")
colnames(Gene) <- c("Gene","Pval")
fwrite(Gene,file = "返修/Model/model_10_KMplot_resullt.csv")

res.cut <- surv_cutpoint(GEO,time = "OS.time",
                         event = "OS",
                         variables = names(GEO)[4:ncol(GEO)],
                         minprop = 0.2)
res_cat<- surv_categorize(res.cut)
mysurv<- Surv(res_cat$OS.time,res_cat$OS)
library(showtext)
showtext_auto()
i="RS"
Gene <- data.frame(Gene=vector(),Pval=vector())
for(i in colnames(res_cat)[3:ncol(res_cat)]){
  group <- res_cat[,i]
  surv_dat <- data.frame(group=group)
  fit <- survfit(mysurv~group)
  
  group <-factor(group,levels = c("low","high"))
  data.survdiff <- survdiff(mysurv~group)
  p.val <- 1-pchisq(data.survdiff$chisq,length(data.survdiff$n)-1)
  HR <- (data.survdiff$obs[2]/data.survdiff$exp[2])/(data.survdiff$obs[1]/data.survdiff$exp[1])
  up95 <- exp(log(HR))+qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1])
  low95 <- exp(log(HR))-qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1])
  xx <- c(i,p.val)
  HR <- paste("Hazard Ratio = ",round(HR,2),sep = "")
  CI <- paste("95% CI: ",paste(round(low95,2),round(up95,2),sep = "-"),sep = "")
  svsort <- GEO[order(GEO[,i]),]
  cutoff <- paste("cutoff = ",round(svsort[fit$n[2],i],2),sep = "")
  
  
  p <- ggsurvplot(fit,data = surv_dat,
                  conf.int = F,
                  censor=F,palette = c("#D53E4F","green"),
                  legend.title = "RSF",size=2,
                  font.x=15,font.y=15,font.title=15,
                  # legend.labs = c(paste0("High ","(",fit$n[1],")"),
                  #                 paste0("Low ","(",fit$n[2],")")),
                  font.lenend=20,
                  legend ="top",
                  legend.abs = c("c2","c1"),
                  xlab = "OS_time(Years)",
                  ylab = "Survival probablity",
                  title = "Survival Curves by Risk Group",
                  break.x.by =5,
                  break.y.by = 0.2,ggtheme = theme_classic()+mytheme,
                  pval = paste(pval=ifelse(p.val < 0.001,"p < 0.001",paste("p = ",round(p.val,3),sep = "")),
                               HR,CI,sep = "\n")
  )
  library(CyDataPro)
  Gene <- rbind(Gene,xx)
  if(p.val < 0.05){
    dir.create("返修2/METABRIC验证",recursive = T)
    CyDataPro::save_plot(p,filename = "返修2/METABRIC验证/KM_diff1",style = "x")
    save_plot(p,filename = paste0("返修/Model/GEO_10_KMplot/",i,"_diff"),style = "x")
  }else{
    save_plot(p,filename = paste0("返修/Model/GEO_10_KMplot/",i,"_Nodiff"),style = "x")
  }
}

dir.create("返修/Model/GEO_10_KMplot")
colnames(Gene) <- c("Gene","Pval")
fwrite(Gene,file = "返修/Model/GEO_10_KMplot_resullt.csv")
library(timeROC)
time_roc_res <- timeROC(T=TCGA$OS.time,delta = TCGA$OS,
                        marker = TCGA$RSF,
                        cause = 1,weighting = "marginal",
                        times = c(1 , 3, 5 ),
                        ROC = T,iid = T)
time_roc_res$AUC
dat = data.frame(fpr = as.numeric(time_roc_res$FP),
                 tpr = as.numeric(time_roc_res$TP),
                 time = rep(as.factor(c(1,3,5)),each = nrow(time_roc_res$TP)))
p <-ggplot() +
  geom_line(data = dat,aes(x = fpr, y = tpr,color = time),size = 2) +
  scale_color_manual(name = NULL,values = c("#92C5DE", "#F4A582", "#66C2A5"),
                     labels = paste0("AUC of ",c(1,3,5),"-y survival: ",
                                     format(round(time_roc_res$AUC,2),nsmall = 2)))+
  geom_line(aes(x=c(0,1),y=c(0,1)),color = "grey")+
   theme_bw()+
  theme(plot.title = element_text(hjust = 0.5,size = 22),
        panel.grid = element_blank(),
        legend.background = element_rect(linetype = 1, size = 0.2, colour = "black"),
        legend.position = c(0.765,0.125),
        axis.text.x=element_text(hjust = 0.5,colour="black",size=16),
        axis.text.y=element_text(hjust = 0.5,colour="black",size=16),
        axis.title.x=element_text(size = 16),
        axis.title.y=element_text(size = 16),
        panel.border = element_blank(),axis.line = element_line(colour = "black",size=1),
        legend.text=element_text(colour="black",size=14),
        legend.title=element_text(colour="black", size=15))+
  scale_x_continuous(expand = c(0.005,0.005))+
  scale_y_continuous(expand = c(0.005,0.005))+
  labs(x = "1 - Specificity",
       y = "Sensitivity",title = "TCGA")+
  coord_fixed()
p
save_plot(p,filename = "返修/Model/RSF_timeROC",style = "g")
tt<- GEO
time_roc_res <- timeROC(T=tt$OS.time,delta = tt$OS,
                        marker = tt$RS,
                        cause = 1,weighting = "aalen",
                        times = c(  21,25,28 ),
                        ROC = T)
time_roc_res$AUC
dat = data.frame(fpr = as.numeric(time_roc_res$FP),
                 tpr = as.numeric(time_roc_res$TP),
                 time = rep(as.factor(c(21,25,28)),each = nrow(time_roc_res$TP)))
p <-ggplot() +
  geom_line(data = dat,aes(x = fpr, y = tpr,color = time),size = 2) +
  scale_color_manual(name = NULL,values = c("#92C5DE", "#F4A582", "#66C2A5"),
                     labels = paste0("AUC of ",c(1,3,5),"-y survival: ",
                                     format(round(time_roc_res$AUC,2),nsmall = 2)))+
  geom_line(aes(x=c(0,1),y=c(0,1)),color = "grey")+
   theme_bw()+
  theme(plot.title = element_text(hjust = 0.5,size = 22),
        panel.grid = element_blank(),
        legend.background = element_rect(linetype = 1, size = 0.2, colour = "black"),
        legend.position = c(0.765,0.125),
        axis.text.x=element_text(hjust = 0.5,colour="black",size=16),
        axis.text.y=element_text(hjust = 0.5,colour="black",size=16),
        axis.title.x=element_text(size = 16),
        axis.title.y=element_text(size = 16),
        panel.border = element_blank(),axis.line = element_line(colour = "black",size=1),
        legend.text=element_text(colour="black",size=14),
        legend.title=element_text(colour="black", size=15))+
  scale_x_continuous(expand = c(0.005,0.005))+
  scale_y_continuous(expand = c(0.005,0.005))+
  labs(x = "1 - Specificity",
       y = "Sensitivity",title = "METABRIC")+
  coord_fixed()
p
xxx$id
library(CyDataPro)
save_plot(p,filename = "返修2/METABRIC验证/METABRIC_RSF_timeROC",style = "g")
save_plot(p,filename = "返修/Model/GEO_RSF_timeROC",style = "g")

fwrite(TCGA,file = "返修/Model/TCGA_model_riskscore.csv")
fwrite(GEO,file = "返修/Model/GEO_model_riskscore.csv")



#######病理参数比较#####
rm(list = ls())

library(gmodels)
library(ggpie)
library(aplot)
install.packages("ggpie")
TCGA<- fread("返修/Model/TCGA_model_riskscore.csv")
tmp<- TCGA %>% column_to_rownames("ID")

clin <- read.delim("../../rawdata/Clin_process.txt",row.names = 1,check.names = F)
clin <- clin[clin$OS>0,]
clin <- clin[clin$Gender=="female",]
clin$status <- ifelse(clin$status=="Alive",1,0)
clin<- clin[rownames(tmp),]
identical(rownames(clin),rownames(tmp))
clin$Group<- ifelse(tmp$RSF> as.numeric(median(tmp$RSF)),"High","Low")
clin$RiskScore<- tmp$RSF


data <- clin %>% mutate(Status = ifelse(status==1 ,"Alive","Dead"),)
colnames(data)
data <- data %>% dplyr::select(RiskGroup,Status,Stage,Age,ER,HER2,PR,PAM50)
data$Age<- ifelse(data$Age>45,">45","<=45")
i=8
ChisqTset <-function(data){
  info <- c("Status","Stage","Age","ER","HER2","PR","PAM50")
  result <- tibble(variable = character(), Pvalue = character(),method = character())
  for (i in 2:ncol(data)) {
    x = info[i-1]
    a = CrossTable(data[,1],data[,i],chisq = T,fisher = F)
    pvalue =a$chisq$p.value
    pcor = a$chisq.corr$p.value
    result = rbind(result,tibble(variable=x,Pvalue=pvalue,method="Chisq"))
  }
  return(result)
}

Chisq_result <- ChisqTset(data)

write.table(Chisq_result,"返修/Model/RS_Clin_chisqTest_result.txt",row.names = F,quote = F,sep = "\t")
High <- data %>% filter(RiskGroup=="High")
Low <- data %>% filter(RiskGroup=="Low")
p1 <-ggdonut(High,group_key = "Status",count_type = "full",
             r0 = 0.5,r1=1,
             border_color = "white",label_type = "none",
             label_pos = "out",label_size = 8)+
  scale_fill_manual(values = c("#999999","#4D4D4D"))+
  theme(text = element_text(size = 15,colour = "black",face = "bold"),
        legend.position = "none")
p2 <- ggdonut(High,group_key = "Stage",count_type = "full",
              r0 = 0.5,r1=1,
              border_color = "white",label_type = "none",
              label_pos = "out",label_size = 8)+
  scale_fill_manual(values = c("#9ECAE1","#C6DBEF","#2171B5","#4292C6","#6BAED6"))+
  theme(text = element_text(size = 15,colour = "black",face = "bold"),
        legend.position = "none")

p3 <- ggdonut(High,group_key = "Age",count_type = "full",
              r0 = 0.5,r1=1,
              border_color = "white",label_type = "none",
              label_pos = "out",label_size = 8)+
  scale_fill_manual(values = c("#33A02C","#B2DF8A"))+
  theme(text = element_text(size = 15,colour = "black",face = "bold"),
        legend.position = "none")
p4 <- ggdonut(High,group_key = "ER",count_type = "full",
              r0 = 0.5,r1=1,
              border_color = "white",label_type = "none",
              label_pos = "out",label_size = 8)+
  scale_fill_manual(values = c("#FF7F00","#FDBF6F","#FEE08B"))+
  theme(text = element_text(size = 15,colour = "black",face = "bold"),
        legend.position = "none")
p5 <-ggdonut(High,group_key = "HER2",count_type = "full",
             r0 = 0.5,r1=1,
             border_color = "white",label_type = "none",
             label_pos = "out",label_size = 8)+
  scale_fill_manual(values = c("#91003F","#E7298A","#C994C7"))+
  theme(text = element_text(size = 15,colour = "black",face = "bold"),
        legend.position = "none")
p6 <- ggdonut(High,group_key = "PR",count_type = "full",
              r0 = 0.5,r1=1,
              border_color = "white",label_type = "none",
              label_pos = "out",label_size = 8)+
  scale_fill_manual(values =  c("#3F007D","#807DBA","#BCBDDC","#9E9AC8","#FCFBFD"))+
  theme(text = element_text(size = 15,colour = "black",face = "bold"),
        legend.position = "none")
p7 <- ggdonut(High,group_key = "PAM50",count_type = "full",
              r0 = 0.5,r1=1,
              border_color = "white",label_type = "none",
              label_pos = "out",label_size = 8)+
  scale_fill_manual(values = c("#016C59","#1C9099","#67A9CF","#A6BDDB","#D0D1E6"))+
  theme(text = element_text(size = 15,colour = "black",face = "bold"),
        legend.position = "none")
p11 <-ggdonut(Low,group_key = "Status",count_type = "full",
             r0 = 0.5,r1=1,
             border_color = "white",label_type = "none",
             label_pos = "out",label_size = 8)+
  scale_fill_manual(values = c("#999999","#4D4D4D"))+labs(fill="Status")+
  theme(text = element_text(size = 15,colour = "black",face = "bold"),
        legend.position = "bottom")
p22 <- ggdonut(Low,group_key = "Stage",count_type = "full",
              r0 = 0.5,r1=1,
              border_color = "white",label_type = "none",
              label_pos = "out",label_size = 8)+
  scale_fill_manual(values = c("#9ECAE1","#C6DBEF","#2171B5","#4292C6","#6BAED6"))+labs(fill="Stage")+
  theme(text = element_text(size = 15,colour = "black",face = "bold"),
        legend.position = "bottom")

p33 <- ggdonut(Low,group_key = "Age",count_type = "full",
              r0 = 0.5,r1=1,
              border_color = "white",label_type = "none",
              label_pos = "out",label_size = 8)+
  scale_fill_manual(values = c("#33A02C","#B2DF8A"))+labs(fill="Age")+
  theme(text = element_text(size = 15,colour = "black",face = "bold"),
        legend.position = "bottom")
p44 <- ggdonut(Low,group_key = "ER",count_type = "full",
              r0 = 0.5,r1=1,
              border_color = "white",label_type = "none",
              label_pos = "out",label_size = 8)+
  scale_fill_manual(values = c("#FF7F00","#FDBF6F","#FEE08B"))+labs(fill="Gender")+labs(fill="ER")+
  theme(text = element_text(size = 15,colour = "black",face = "bold"),
        legend.position = "bottom")
p55 <-ggdonut(Low,group_key = "HER2",count_type = "full",
             r0 = 0.5,r1=1,
             border_color = "white",label_type = "none",
             label_pos = "out",label_size = 8)+
  scale_fill_manual(values = c("#91003F","#E7298A","#C994C7"))+labs(fill="HER2")+
  theme(text = element_text(size = 15,colour = "black",face = "bold"),
        legend.position = "bottom")
p66 <- ggdonut(Low,group_key = "PR",count_type = "full",
              r0 = 0.5,r1=1,
              border_color = "white",label_type = "none",
              label_pos = "out",label_size = 8)+
  scale_fill_manual(values =  c("#3F007D","#807DBA","#BCBDDC","#9E9AC8","#FCFBFD"))+labs(fill="PR")+
  theme(text = element_text(size = 15,colour = "black",face = "bold"),
        legend.position = "bottom")
p77 <- ggdonut(Low,group_key = "PAM50",count_type = "full",
              r0 = 0.5,r1=1,
              border_color = "white",label_type = "none",
              label_pos = "out",label_size = 8)+
  scale_fill_manual(values = c("#016C59","#1C9099","#67A9CF","#A6BDDB","#D0D1E6"))+labs(fill="PAM50")+
  theme(text = element_text(size = 15,colour = "black",face = "bold"),
        legend.position = "bottom")




all <-ggarrange(p1,p2,p3,p4,p5,p6,p7,p11,p22,p33,p44,p55,p66,p77,ncol = 7,nrow = 2)
save_plot(xx,filename = "ML/RS_Clin_chisq_PiePlot",style = "g",width = 15,height = 8)
all
xx <-annotate_figure(all,top = text_grob("High",color = "black",face = "bold",size = 14,x=0,y=0.6,hjust = 1),
                     left = text_grob("Low",color = "black",face = "bold",size = 14,y=0.5))

dir.create("返修/Model/Clin")
save_plot(xx,filename = "返修/Model/RS_Clin_chisq_PiePlot",style = "g",width = 30,height = 10)



data <- data %>% dplyr::select(RiskScore,Stage,T_stage,M_stage,N_stage)

data$Stage_n<- recode(data$Stage,"I"="I&II","II"="I&II",
                      "III"="III&IV","IV"="III&IV")
data$T_stage_n<- recode(data$T_stage,"T1a"="T1","T1b"="T1","T1c"="T1",
                      "T2a"="T2","T2b"="T2","T3a"="T3","T4b"="T4","T4d"="T4")
data$T_stage_nn<- recode(data$T_stage_n,"T1"="T1&T2","T2"="T1&T2","T3"="T3&T4",
                        "T4"="T3&T4")
data$M_stage_n<- recode(data$M_stage,"cM0 (i+)"="M0")
data$N_stage_n<- recode(data$N_stage,"N0 (i-)"="N0","N0 (i+)"="N0","N0 (mol+)"="N0",
                        "N1a"="N1","N1b"="N1","N1c"="N1","N1mi"="N1",
                        "N2a"="N2","N3a"="N3","N3b"="N3","N3c"="N3")

data$N_stage_nn<- recode(data$N_stage_n,"N0"="N0&N1","N1"="N0&N1",
                         "N2"="N2&N3","N3"="N2&N3")

data<- data %>% filter(Stage !="X",T_stage_n != "TX",M_stage !="MX",N_stage !="NX")

p<-ggplot(data,aes(x=N_stage_nn,y=RiskScore,fill=N_stage_nn))+
  geom_boxplot(size=1)+labs(x="N_stage",fill="")+
  stat_compare_means(size=6,label.x.npc = "center",label = "..p.format..")+
  theme_test(base_rect_size = 2)+mytheme
dir.create("返修2/风险模型Stage比较")
save_plot(p,filename = "返修2/风险模型Stage比较/N_stage_notdiff",style = "g")

p<-ggplot(data,aes(x=Stage_n,y=RiskScore,fill=Stage_n))+
  geom_boxplot(size=1)+labs(x="Stage",fill="")+
  stat_compare_means(size=6,label.x.npc = "center",label = "..p.format..")+
  theme_test(base_rect_size = 2)+mytheme
p

save_plot(p,filename = "返修2/风险模型Stage比较/Stage_notdiff",style = "g")

data$T_stage_nn<- ifelse(data$T_stage_nn=="T1&T2","T3&T4","T1&T2")
p<-ggplot(data,aes(x=T_stage_nn,y=RiskScore,fill=T_stage_nn))+
  geom_boxplot(size=1)+labs(x="T_stage",fill="")+
  stat_compare_means(size=6,label.x.npc = "center",label = "..p.format..")+
  theme_test(base_rect_size = 2)+mytheme
p

save_plot(p,filename = "返修2/风险模型Stage比较/T_stage_diff",style = "g")

data$M_stage_n<- ifelse(data$M_stage_n=="M0","M1","M0")
p<-ggplot(data,aes(x=M_stage_n,y=RiskScore,fill=M_stage_n))+
  geom_boxplot(size=1)+labs(x="T_stage",fill="")+
  stat_compare_means(size=6,label.x.npc = "center",label = "..p.format..")+
  theme_test(base_rect_size = 2)+mytheme
p

save_plot(p,filename = "返修2/风险模型Stage比较/M_stage_notdiff",style = "g")



####新辅助治疗
data<- fread("https://gdc-hub.s3.us-east-1.amazonaws.com/download/TCGA-BRCA.clinical.tsv.gz",data.table = F)
xx<- data %>% filter(treatment_or_therapy.treatments.diagnoses !="['not reported', 'not reported']")
group<- ifelse(str_detect(xx$treatment_or_therapy.treatments.diagnoses,"yes"),"Yes","No")
table(group)

tmp<- xx %>% dplyr::select(sample,treatment_type.treatments.diagnoses,treatment_or_therapy.treatments.diagnoses)
tmp$It1<-str_split_fixed(tmp$treatment_type.treatments.diagnoses,"[',]",4) [,2]
tmp$n1<-str_split_fixed(tmp$treatment_or_therapy.treatments.diagnoses,"[',]",4) [,2]

tmp$It2<-str_split_fixed(tmp$treatment_type.treatments.diagnoses,"[,]",4) [,3]
tmp$It2<-str_remove(tmp$It2,"'")
tmp$n2<-str_split_fixed(tmp$treatment_or_therapy.treatments.diagnoses,"[,']",6) [,5]
tmp<- tmp[,-c(2,3)]

tmp$It1<- paste0(tmp$It1,"-",tmp$n1)
tmp$It2<- paste0(tmp$It2,"-",tmp$n2)

rownames(tmp)<- tmp$sample
inn<- intersect(rownames(tmp),rownames(clin))
tmp<- tmp[inn,]
clin<- clin[inn,]
tmp$RiskScore<- clin$RiskScore

x1<- tmp %>% filter(It1 =="Radiation Therapy-yes")
x2<- tmp %>% filter(It2 =="Radiation Therapy-yes")

x1<- tmp %>% filter(It1 =="Pharmaceutical Therapy-yes")
x2<- tmp %>% filter(It2 =="Pharmaceutical Therapy-yes")
xx<- rbind(x1,x2)

tmp$Group<- ifelse(tmp$sample %in% xx$sample,"Yes","No")
tmp$Group<- ifelse(tmp$n1 =="yes" |tmp$n2=="yes","Yes","No")
table(tmp$Group)
#No Yes 
#108 870 
library(pROC)
roc<- roc(tmp$Group,tmp$RiskScore,levels=c("No", "Yes"), direction=">",smooth=F,ci=T)
roc$auc
roc$ci
aucs <- round(ci.auc(roc),3);aucs
lb = paste0("AUC:", aucs[2], "\n",
            "95%CI:", aucs[1], "-", aucs[3])
p<-ggroc(roc, color = "darkred", legacy.axes = T, linewidth = 2) +
  annotate("segment", x = 0, y = 0, xend = 1, yend = 1, linetype = "dashed") +
  coord_cartesian(xlim = c(-0.01, 1.01), ylim = c(-0.01, 1.01), expand = F) +
  annotate("text", x = 0.9, y = 0.25, label = lb, hjust = 1,size=6,color="red") +
  theme_test(base_rect_size = 2) +
  theme(panel.grid = element_blank())+mytheme
dir.create("返修2/新辅助医疗")
save_plot(p,filename = "返修2/新辅助医疗/ROC",style = "g")
####表达量比较###
rm(list = ls())
load("../../rawdata//DEG.Rdata")
exp <- read.delim("../../rawdata/exprmatx_without_normal.txt",check.names = F)
Top20 <- xxx$id

gene_exp <- exp[Top20,]
a <-sample_group %>% as.data.frame()
clin <- read.delim("ML/Clin_with_Risk.txt",row.names = 1)
group <- data.frame(Sample = rownames(clin),Group=clin$Group)
write.table(group,file = "../../rawdata/TCGA/TCGA-BRCA/Sample_Group.txt",row.names = F,sep = "\t",quote = F)
gene_exp <- gene_exp %>% rownames_to_column("Gene")
data<-reshape2::melt(gene_exp,id.vars = "Gene")
names(data) <- c("Gene","Sample","values")
data1 <- inner_join(data,group)
p <-ggplot(data1,aes(Gene,values,fill=Group))+
  geom_boxplot(color="black",alpha=0.6,outlier.shape = 21,outlier.size = 1.2)+
  scale_fill_manual(values = c("#d9352a","#4979b6"))+
  labs(x="Gene",y="Expression")+
  stat_compare_means(aes(group= Group),label = "p.signif",size=5)+theme_bw()+mytheme
  theme(axis.text.x = element_text(angle = 45))
p
dir.create("返修/TCGA")
save_plot(p,filename = "返修/TCGA/KeyGene_expre_RsGroup_boxplot",style = "g",width = 10,height = 6)

gene_exp<- exprmatix[gene$id,] %>% t() %>% as.data.frame()
gene_exp$Group <- ifelse(as.numeric(str_sub(rownames(gene_exp), 
                                            14, 15)) < 10, "Tumor", "Normal")

data<-reshape2::melt(gene_exp,id.vars = "Group")
names(data) <- c("Group","Gene","values")
data1 <- inner_join(data,group)
p <-ggplot(data,aes(Gene,values,fill=Group))+
  geom_boxplot(color="black",alpha=0.6,outlier.shape = 21,outlier.size = 1.2)+
  scale_fill_manual(values = c("#d9352a","#4979b6"))+
  labs(x="Gene",y="Expression")+
  stat_compare_means(aes(group= Group),label = "p.signif",size=5)+theme_bw()+mytheme
dir.create("返修/TCGA/expression")
save_plot(p,filename = "返修/TCGA/expression/KeyGene_expre_Normal_Tumor_boxplot",style = "g",width = 10,height = 6)
p

group <- group %>% column_to_rownames("Sample")
# library(ComplexHeatmap)
# anncol <- HeatmapAnnotation(Group = group$Group,
#                             col = list(RsGroup = c("Tumor"= "#4979b6", "normal"="#d9352a")))
# plot_data <- t(gene_exp)
# plot_data <- plot_data[,rownames(group)]
# identical(colnames(plot_data),rownames(group))
# data <- apply(plot_data,1,scale)
# rownames(data) <- colnames(plot_data)
# Heatmap(plot_data,na_col = "white",
#             show_column_names = F,
#             show_row_names = T,name = "GSVA Score",
#             column_order = c(colnames(plot_data)[c(grep("Tumor",group$Group),grep("normal",group$Group))]),
#             column_split = group$Group,
#             cluster_columns = F,column_title = NULL,
#             top_annotation = anncol)

##体细胞突变###
library(devtools)
library(data.table)
library(maftools)
library(tidyverse)
rm(list=ls())
maf <- fread("../../rawdata/TCGA.BRCA.mutect.995c0111-d90b-4140-bee7-3845436c3b42.DR-10.0.somatic (1).maf.gz")
maf$Tumor_Sample_Barcode <- substr(maf$Tumor_Sample_Barcode,1,16)
# clin <- read.delim("ML/Clin_with_Risk.txt",check.names = F,row.names = 1)
clin <- clin %>% rownames_to_column("Tumor_Sample_Barcode")
mafd <- read.maf(maf,clin)
vc_cols = c("#1F6CB6","red3","#70BC6B","#F4B9C5","#784194","#B81A7B","#A65628","#9E1F63","#005A32")  # RdYlBu
names(vc_cols) = c(
  'Missense_Mutation',
  'Multi_Hit',
  'Frame_Shift_Del',
  'Nonsense_Mutation',
  'Frame_Shift_Ins',
  'In_Frame_Ins',
  'Splice_Site',
  'In_Frame_Del',
  'Nonstop_Mutation'
)
gene <- fread("返修/Model/Gene_Univar_Cox_result.csv") %>% filter(pvalue < 0.05)
dir.create("返修/TCGA/Mutation")
pdf("返修/TCGA/Mutation/Key_gene_mutation_oncoplot.pdf",width = 8,height = 6)
oncoplot(mafd,clinicalFeatures = "Group",sortByAnnotation = T,genes = gene$id,
         showTumorSampleBarcodes = F,
         bgCol = "gray",colors = vc_cols,fontSize = .7)
dev.off()
pdf("返修/TCGA/Mutation/Top30_gene_mutation_oncoplot.pdf",width = 8,height = 6)
oncoplot(mafd,clinicalFeatures = "Group",sortByAnnotation = T,top = 30,
         showTumorSampleBarcodes = F,removeNonMutated = T,borderCol = NA,colors = vc_cols,fontSize = .7,sepwd_samples = 0)

dev.off()

###TMB###
library(dplyr)
tmb <- tmb(maf = mafd)
group <- clin %>% dplyr::select("Tumor_Sample_Barcode","Group")

tmbGroping<- inner_join(tmb,group) %>% filter(total != 0)
data_p <- tmbGroping %>% dplyr::select(total_perMB,Group) %>% melt()
p <-ggplot(data_p,aes(Group,value,fill=Group))+
  geom_violin(trim = T,aes(fill=Group),color="white",cex=1,alpha=0.4)+
  geom_boxplot(width=0.15,position = position_dodge(0.9),color="white")+
  # geom_signif(comparisons = list(c("High","Low")),step_increase = 0.1,size=1,test = "wilcox.test",textsize = 5)+
  #  geom_jitter()+
  stat_compare_means(aes(group=Group),size=5)+
  ylab("TMB")+xlab("Risk Group")+ylim(0.8,1.05)+
   theme_classic()+mytheme+
  scale_fill_manual(values = c("#CA0020","#7FBC41"))
p
dir.create("TCGA/TMB")
save_plot(p,filename = "TCGA/TMB/RsGroup_TMB_boxplot",style = "g")

##预后
library(survival)
library(survminer)
group <- clin %>% select("Tumor_Sample_Barcode","status","OS","Group","risk_score")
survtmb <-inner_join(tmb,group) %>% filter(total != 0)
survtmb$TMB_Group <-if_else(survtmb$total_perMB >= median(survtmb$total_perMB), "High", "Low")

fit <- survfit(Surv(OS,status)~TMB_Group,data = survtmb)

p <- ggsurvplot(fit,data = survtmb,
                censor.shape="|",censor.size=4,
                # conf.int = T,
                # con.int.style ="ribnbon",
                # con.in.alpha = 0.2,
                pval = TRUE,
                # pval.coord = c(800,0.1),
                # pval.method = T,
                # pval.method.coord = c(1,0.1),
                #palette = "lancet",
                surv.median.line = "hv",
                ggtheme = theme_prism(),
                legend =c(0.8,0.5),
                legend.abs = c("High","Low"),
                xlab = "OS_time(Days)",
                ylab = "Survival probablity",
                title = "Survival Curves by TMB Group",
                break.x.by =1000,
                break.y.by = 0.2,
                risk.table = F,
                # risk.table.col = "strata",
                # risk.table.height = 0.2,
                risk.table.y.text = FALSE)
p
save_plot(p ,filename = "TCGA/TMB/Surv_TMBGroup_kvPlot",style = "xx")

survtmb$ConGroup <- ifelse(survtmb$risk_score >= median(survtmb$risk_score)& survtmb$total_perMB >= median(survtmb$total_perMB),"riskHigh-TMBHigh",
                           ifelse(survtmb$risk_score < median(survtmb$risk_score)&survtmb$total_perMB >= median(survtmb$total_perMB),"riskLow-TMBHigh",
                                  ifelse(survtmb$risk_score >= median(survtmb$risk_score)&survtmb$total_perMB < median(survtmb$total_perMB),"riskHigh-TMBLow","riskLow-TMBLow")))
table(survtmb$ConGroup)
fit <- survfit(Surv(OS,status)~ConGroup,data = survtmb)

p <- ggsurvplot(fit,data = survtmb,
                censor.shape="|",censor.size=4,
                # conf.int = T,
                # con.int.style ="ribnbon",
                # con.in.alpha = 0.2,
                pval = TRUE,
                # pval.coord = c(800,0.1),
                # pval.method = T,
                # pval.method.coord = c(1,0.1),
                #palette = "lancet",
                surv.median.line = "hv",
                ggtheme = theme_prism(),
                legend =c(0.8,0.5),
                legend.abs = c("riskHigh-TMBHigh","riskHigh-TMBLow","riskLow-TMBHigh","riskLow-TMBLow"),
                xlab = "OS_time(Days)",
                ylab = "Survival probablity",
                title = "Survival Curves by riskScore-TMB Group",
                break.x.by =1000,
                break.y.by = 0.2,
                risk.table = F,
                # risk.table.col = "strata",
                # risk.table.height = 0.2,
                risk.table.y.text = FALSE)
p
save_plot(p,filename = "TCGA/TMB/Surv_risk-TMB_plot",style = "xx")
write.table(tmbGroping,file = "TCGA/TMB/TMB_result.txt",sep = "\t",quote = F,row.names = F)
##拷贝数变异
library(ChIPseeker)
library(GenomicRanges)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library(bigmemory)
library(tidyverse)
install.packages("bigmemory")
BiocManager::install("TxDb.Hsapiens.UCSC.hg38.knownGene")
setwd("D:/WORK/Project/P2")
rm(list = ls())
data <- read.delim(gzfile("../../rawdata/TCGA-BRCA.ggene_cnv.tsv.gz"),check.names = F,row.names = 1)

# Key_gene <- read.delim("Ml/rf_Importance.txt")[1:20,] 
# cnv_dat <- data[data$`Gene Symbol` %in% Key_gene$Variable,]
# 
# 
# p <- data
# p$Chrom <- paste0("chr",p$Chrom)
# peak <- GRanges(sample = p[,1],
#                 Segment_Mean = p[,5],
#                 seqnames = Rle(p[,2]),
#                 ranges = IRanges(p[,3],p[,4]),
#                 strand = rep(c("*"),nrow(p)))
# txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
# gc()
# memory.limit()
# rm(data)
# rm(p)
# object.size(peak)
# options(future.globals.maxSize = 1000 * 1024^3)
# peakAnno <- annotatePeak(peak,tssRegion = c(-3000,3000),
#                          TxDb =TxDb.Hsapiens.UCSC.hg38.knownGene ,annoDb = "org.Hs.eg.db")

id <- read.delim("../../rawdata/gencode.v22.annotation.gene.probeMap",row.names = 1)
geneinfo <- na.omit(id[,1:2])
int <- intersect(rownames(geneinfo),rownames(data))
geneinfo <- geneinfo[int,] %>% na.omit() #基因信息行数与表达矩阵匹配
geneinfo <- geneinfo[!duplicated(geneinfo$gene),] #去除重复
data1 <- data[rownames(geneinfo),]
which(is.na(geneinfo$gene))
tail(geneinfo)
class(geneinfo)
class(expr)
rownames(data1) <- geneinfo$gene 

gene_data <- gene$id
clin <- read.delim("ML/Clin_with_Risk.txt")
High <- clin %>% filter(Group == "High")
Low <- clin %>% filter(Group == "Low")
int <- intersect(gene_data,rownames(data1))
cnv_dat <- data1[int,]
cnv_dat[cnv_dat > 0.3] =0.4
cnv_dat[cnv_dat < -0.3] = -0.4
cnv_dat[cnv_dat < 0.3 & cnv_dat > -0.3 ]=0

cnv_dat[cnv_dat == 0.4] <- "Amplifications"
cnv_dat[cnv_dat == -0.4 ] <- "Deletions"
cnv_dat[cnv_dat == 0 ] <- "neutral"

cnv <- as.data.frame(t(cnv_dat))
cnv <- as.data.frame(t(apply(cnv,2,table)))
dat <- cnv[,1:2]
dat$Gene <- rownames(dat)
df <- reshape2::melt(dat)
library(ggalt)

p <-ggplot(dat,aes(y=Gene))+
  geom_point(data = df,aes(x=value,color=variable),size=3)+
  geom_dumbbell(aes(x=Deletions,xend =Amplifications),size=3,color = "#e3e2e1",
                colour_x = "green",colour_xend = "red",
                dot_guide = T,dot_guide_size = 0.25)+
  theme_bw()+coord_flip()+
  scale_color_manual(name = "",values = c("green","red"))+mytheme+
  theme(axis.text.x = element_text(angle = 45,hjust = 1))+
  xlab("CNV")
p
dir.create("返修/TCGA/CNV")
save_plot(p ,filename = "返修/TCGA/CNV/ALL_CNV_plot",style = "g")
write.table(dat,file = "返修/TCGA/CNV/All_CNV_result.txt",row.names = T,quote = F,sep = "\t")

cnv_High <- cnv_dat[,intersect(colnames(cnv_dat),High$Tumor_Sample_Barcode)]
cnv_Low <- cnv_dat[,intersect(colnames(cnv_dat),Low$Tumor_Sample_Barcode)]
dat_H <- as.data.frame(t(cnv_High))
cnv_H <- as.data.frame(t(apply(dat_H,2,table)))
cnv_H <- cnv_H[,1:2]
cnv_H$Gene <- rownames(cnv_H)
cnv_HD <- reshape2::melt(cnv_H)
cnv_HD$Group <- "High"

dat_L <- as.data.frame(t(cnv_Low))
cnv_L <- as.data.frame(t(apply(dat_L,2,table)))
cnv_L <- cnv_L[,1:2]
cnv_L$Gene <- rownames(cnv_L)
cnv_LD <- reshape2::melt(cnv_L)
cnv_LD$Group <- "High"

All <- rbind(cnv_H,cnv_L)
p <-ggplot(cnv_L,aes(y=Gene))+
  geom_point(data = cnv_LD,aes(x=value,color=variable),size=3)+
  geom_dumbbell(aes(x=Deletions,xend =Amplifications),size=3,color = "#e3e2e1",
                colour_x = "green",colour_xend = "red",
                dot_guide = T,dot_guide_size = 0.25)+
  theme_bw()+coord_flip()+
  scale_color_manual(name = "",values = c("green","red"))+mytheme+
  theme(axis.text.x = element_text(angle = 45,hjust = 1))+
  xlab("CNV")+
  ggtitle("Low Risk")
p
p <- save_plot(p ,filename = "返修/TCGA/CNV/LowRisk_CNV_plot",style = "g")
write.table(cnv_L,file = "返修/TCGA/CNV/LowRisk_CNV_result.txt",row.names = T,quote = F,sep = "\t")

p <-ggplot(cnv_H,aes(y=Gene))+
  geom_point(data = cnv_HD,aes(x=value,color=variable),size=3)+
  geom_dumbbell(aes(x=Deletions,xend =Amplifications),size=3,color = "#e3e2e1",
                colour_x = "green",colour_xend = "red",
                dot_guide = T,dot_guide_size = 0.25)+
  theme_bw()+coord_flip()+
  scale_color_manual(name = "",values = c("green","red"))+mytheme+
  theme(axis.text.x = element_text(angle = 45,hjust = 1))+
  xlab("CNV")+
  ggtitle("High Risk")
p
save_plot(p ,filename = "返修/TCGA/CNV/HighRisk_CNV_plot",style = "g")
write.table(cnv_L,file = "返修/TCGA/CNV/HighRisk_CNV_result.txt",row.names = T,quote = F,sep = "\t")

####甲基化##########

rm(list = ls())
setwd("D:/WORK/Project/乳腺癌NC细胞/")
library(data.table)
library(tidyverse)
library(ggsignif)
library(ggprism)
library(ggpubr)
gc()
rm()
TCGA<- fread("返修/Model/TCGA_model_riskscore.csv")
tmp<- TCGA %>% column_to_rownames("ID")

clin <- read.delim("../../rawdata/Clin_process.txt",row.names = 1,check.names = F)
clin <- clin[clin$OS>0,]
clin <- clin[clin$Gender=="female",]
clin$status <- ifelse(clin$status=="Alive",1,0)
clin<- clin[rownames(tmp),]
identical(rownames(clin),rownames(tmp))
clin$Group<- ifelse(tmp$RSF> as.numeric(median(tmp$RSF)),"High","Low")

get_miss_rate <- function(x){sum(is.na(x))/length(x)*100}
data <- fread("../../rawdata/TCGA-BRCA.methylation450.tsv.gz",data.table = F)
idmap <- fread("../../rawdata/illuminaMethyl450_hg38_GDC",data.table = F)
key_gene <- fread("返修/Model/Gene_Univar_Cox_result.csv") %>% filter(pvalue<0.05)

rownames(data)<- data[,1]
data <- data[,-1]
data <- data %>% rownames_to_column("id")
miss_prob <- data %>% apply(1, get_miss_rate) %>% as.data.frame() %>% rename(x=".")
prob <- miss_prob %>% filter(x<10) %>% rownames()
data<- data %>% na.omit()
data_clean <- data %>% rownames_to_column("id")

Pro_ID <- select(idmap,c("#id","gene")) 
Pro_ID <-Pro_ID[Pro_ID$gene!='.',] #去掉空值

Pro_ID$gene<- sapply(Pro_ID$gene, function(x){
  unlist(str_split(x,","))[1]
})
colnames(Pro_ID) <- c("id","Geneid")
exp_symbol <- merge(Pro_ID,data,by = "id") %>%
  select(-id) %>% group_by(Geneid) %>% 
  summarise_all(max) %>% filter(!Geneid =='')
Methy_data <- exp_symbol %>% column_to_rownames("Geneid")
s <- unique(exp_symbol$Geneid[duplicated((exp_symbol$Geneid))])
inner<- intersect(key_gene$id,rownames(Methy_data))
Methy_key <- Methy_data[inner,]
miss_prob <- Methy_key %>% apply(1, get_miss_rate) %>% as.data.frame() %>% rename(x=".")
prob <- miss_prob %>% filter(x<10) %>% rownames()
Methy_key <- Methy_key[prob,] 

cluster <- clin
inner <- intersect(rownames(cluster),colnames(Methy_key)) 
data_f <- Methy_key[,inner] %>% rownames_to_column("Geneid")
cluster1 <- cluster[inner,] %>% rownames_to_column("variable")
cluster1 <- cluster1 %>% select(variable,Group)
data1 <- melt(data_f)
data_p <- inner_join(data1,cluster1) %>% na.omit()
which(is.na(data_p$value))

p <-ggplot(data_p,aes(Geneid,value,fill=Group))+
  geom_boxplot(color="black",alpha=0.6,outlier.shape = 21,outlier.size = 1.2)+
  scale_fill_manual(values = c("#E41A1C","#4DAF4A"))+
  labs(x="Gene",y="Methylation Leve")+
  stat_compare_means(aes(group= Group),label = "p.signif",size=5)+theme_bw()+mytheme

p
library(CyDataPro)
dir.create("返修/TCGA/Methy")
save_plot(p,filename = "返修/TCGA/Methy//Gene_Methy_boxplot",style = "g")

p <-ggplot(data_p,aes(Group,value,fill=Group))+
  geom_violin(trim = F,aes(file=Group),color="white",cex=1,alpha=0.5)+
  geom_boxplot(width=0.1,position = position_dodge(0.9),color="white")+
  geom_signif(comparisons = list(c("High","Low")),step_increase = 0.1,size=1,test = "wilcox.test",textsize = 5)+
  #  geom_jitter()+
  #  stat_compare_means(aes(group=Cluster),label.x = 0.6,size=5,)+
  ylab("Methylation Level")+xlab("Risk Subtype")+ylim(0.8,1.05)+
  theme_bw()+mytheme+
  scale_fill_manual(values = c("#B3CDE3","#DECBE4"))
p
save_plot(p,filename = "返修/TCGA/Methy/Group_keyGene_Methy_boxplot",style = "g")  
write.table(data_p,file = "返修/TCGA/Methy/KeyGene_MethyLevel.txt",row.names = F,quote = F,sep = "\t")



###独立预后###
rm(list = ls())
library(survival)
library(survminer)
library(plyr)
library(ezcox)
library(tidyverse)
install.packages("ezcox")
setwd("D:/WORK/Project/P2/")

edit(UCSC_clinData_Normalize)
phe<- fread("../../rawdata/TCGA-BRCA.GDC_phenotype.tsv.gz",data.table = F)
phe$age_at_initial_pathologic_diagnosis
phe<- phe %>% column_to_rownames("submitter_id.samples")
colnames_num <- grep('receptor_status',colnames(phe))
phenotype_colnames <- colnames(phe)[colnames_num]
eph <- phe[,c(5,colnames_num[1:3])]
table(eph$metastatic_breast_carcinoma_progesterone_receptor_status)
library(TCGAbiolinks)
subtypes <- PanCancerAtlas_subtypes()
DT::datatable(subtypes,
              filter = 'top',
              options = list(scrollX = TRUE, keys = TRUE, pageLength = 10),
              rownames = FALSE)
BRCA_subtypes <- as.data.frame(subtypes[subtypes$cancer.type == "BRCA",]) 
clinical3 <- data.frame(SampleNames = BRCA_subtypes$pan.samplesID,
                        PAM50 = BRCA_subtypes$Subtype_mRNA,
                        SampleType = sapply(strsplit(BRCA_subtypes$pan.samplesID,'-'),"[",4))
clinical3$PatientNames <- substr(clinical3$SampleNames,1,16)
clinical3$Group <- ifelse(clinical3$SampleType > 10,"Normal","Tumor")
clinical3 <- clinical3[clinical3$Group == "Tumor",] 
clinical3<- clinical3 %>% column_to_rownames("PatientNames")
rownames(clinical3)<- clinical3$PatientNames

inn<- intersect(rownames(eph),rownames(clinical3))
clinical3<- clinical3[inn,]
eph<- eph[inn,]
colnames(eph)<- c("Age1","ER","HER2","PR")
eph$PAM50<- clinical3$PAM50

table(clinical3$PAM50)
TCGA<- fread("返修/Model/TCGA_model_riskscore.csv")
tmp<- TCGA %>% column_to_rownames("ID")

clin <- read.delim("../../rawdata/Clin_process.txt",row.names = 1,check.names = F)
clin <- clin[clin$OS>0,]
clin <- clin[clin$Gender=="female",]
clin$status <- ifelse(clin$status=="Alive",1,0)
clin<- clin[rownames(tmp),]
identical(rownames(clin),rownames(tmp))
clin$RiskGroup<- ifelse(tmp$RSF> as.numeric(median(tmp$RSF)),"High","Low")
clin$RiskScore<- tmp$RSF

inn<- intersect(rownames(clin),rownames(eph))
clin<- clin[inn,]
eph<- eph[inn,]
clin<- cbind(clin,eph)
clin$ER<- ifelse(clin$ER=="","No",clin$ER)
clin$PR<- ifelse(clin$PR=="","No",clin$PR)
clin$HER2<- ifelse(clin$HER2=="","No",clin$HER2)
clin$Age<- clin$Age1
fwrite(clin,file = "BRCA_clin.csv",row.names = T)


tmp<- clin
tmp$Stage <- toupper(tmp$Stage)
tmp <- tmp %>% mutate(Stage = case_when(Stage %in% c("I", 
                                                     "II") ~ "I&II", Stage %in% 
                                          c("STAGE II", "STAGE IIA", "STAGE IIB", "STAGE IIC") ~ 
                                          "I&II", Stage %in% c("III", "IV", "STAGE IIIB", 
                                                             "STAGE IIIC") ~ "III&IV", Stage %in% c("STAGE IV", 
                                                                                                 "STAGE IVA", "STAGE IVB", "STAGE IVC") ~ "III&IV", Stage %in% 
                                          c("STAGE X", "NOT REPORTED") ~ "X", .default = Stage))
tmp <- tmp %>% filter(Stage != "X")
table(tmp$Stage)




  tmp$T_stage <- toupper(tmp$T_stage)
  tmp <- tmp %>% mutate(T_stage = case_when(T_stage %in% 
                                              c("T1", "T1A", "T1B", "T1C") ~ "T1&T2", T_stage %in% 
                                              c("T2", "T2A", "T2B", "T2C") ~ "T1&T2", T_stage %in% 
                                              c("T3", "T3A", "T3B", "T3C") ~ "T3&T4", T_stage %in% 
                                              c("TX", "NOT REPORTED", "") ~ "TX", T_stage %in% 
                                              c("T4","T4B","T4D")~ "T3&T4",
                                            .default = T_stage))
  tmp <- tmp %>% filter(T_stage != "TX")
table(tmp$T_stage)

tmp$M_stage <- toupper(tmp$M_stage)
tmp <- tmp %>% mutate(M_stage = case_when(M_stage %in% 
                                            c("CM0 (I+)", "M0") ~ "M0", M_stage %in% 
                                            c("M1", "MX") ~ "M1", 
                                          .default = M_stage))
tmp$N_stage <- toupper(tmp$N_stage)
tmp <- tmp %>% mutate(N_stage = case_when(N_stage %in% 
                                            c("N0", "N0 (I-)", "N0 (I+)", "N0 (MOL+)") ~ "N0&N1", N_stage %in% 
                                            c("N1", "N1A", "N1B", "N1C","N1MI") ~ "N0&N1", N_stage %in% 
                                            c("N2", "N2A") ~ "N2&N3", N_stage %in% 
                                            c("N3", "N3A", "N3B","N3C") ~ "N2&N3",
                                          .default = N_stage))

tmp <- tmp %>% filter(N_stage != "NX")
colnames(tmp)
table(tmp$N_stage)

variable.names<- colnames(tmp)[c(4,6:9,12)]

result  <- ezcox(tmp,time = "OS",status = "status",covariates = variable.names,return_models = T)
dir.create("返修/TCGA/IndependentFactor")
mod <- get_models(result)
a <- result$res %>% as.data.frame()
write.table(a,file = "返修/TCGA/IndependentFactor/Single_cox_result.txt",row.names = F,sep = "\t",quote = F)
pdf("返修/TCGA/IndependentFactor/Singlle.pdf",width = 8,height = 6)
show_models(mod)
dev.off()
png("返修/TCGA/IndependentFactor/Singlle.png",width = 8,height = 6,units = "in",res = 300)
show_models(mod)
dev.off()

td <- cbind(tmp[,2:3],tmp[,variable.names])
tdmultiCox=coxph(Surv(OS, status) ~ ., data = td)
x <-summary(tdmultiCox) 
pvalue=signif(as.matrix(x$coefficients)[,5],2)
HR=signif(as.matrix(x$coefficients)[,2],2)
low=signif(x$conf.int[,3],2)
high=signif(x$conf.int[,4],2)
multi_res=data.frame(p.value=pvalue,
                     HR=paste(HR," (",low,"-",high,")",sep=""),
                     stringsAsFactors = F
)
write.table(multi_res,file = "返修/TCGA/IndependentFactor/multi_cox_result.txt",row.names = F,quote = F,sep = "\t")
p <-ggforest(tdmultiCox,
             main = "Hazard ratio",
             cpositions = c(0.02,0.15, 0.35), #前三列的位置，第二列是样品数，设了个负值，相当于隐藏了
             fontsize = 1.3, #字体大小
             refLabel = "reference", 
             noDigits = 2)

p
save_plot(p,filename = "返修/TCGA/IndependentFactor/Multi",style = "g",width = 12)

variable.names <- c("N_stage","RiskScore")
td <- cbind(tmp[,2:3],tmp[,variable.names])
tdmultiCox=coxph(Surv(OS, status) ~ ., data = td)




library(ggprism)
library(ggDCA)
library(regplot)



pdf(file = "TCGA/IndependentFactor/Nomogram_plot",width = 8,height = 8,onefile = FALSE)
regplot(tdmultiCox,observation = T,
        failtime = c(365,365*2),
        prfail = T,showP = T,
        droplines=F,
        rank = "sd",
        interval = "confidence",points = T,title = "Nomogram")
  dev.off()
save_plot(p,filename = "TCGA/Lasso_cox/Nomogram_plot",style = "x",width = 8,height = 6)


library(nomogramFormula)

library(rms)
library(survivalROC)
td$N_stage <- ifelse(td$N_stage=="N0&N1",1,0)

td$RiskGroup <- ifelse(td$RiskScore> median(td$RiskScore),1,0)

td$OS_t <- td$OS/365
dd <- datadist(td)
options(datadist="dd")
str(td)
f <- cph(formula(Surv(OS_t, status) ~N_stage+RiskScore), data = td, x = TRUE, y = TRUE,
         surv = TRUE)

surv <- Surv(f)
nomo <- nomogram(f)
plot(nomo)
result <- formula_rd(nomo)
result$formula
td$points <- points_cal(formula = result$formula,lp=f$linear.predictors)

Score <- coxph(Surv(OS, status) ~ RiskScore, data = td)
Stage <- coxph(Surv(OS, status) ~ N_stage, data = td)
Nomo<- coxph(Surv(OS, status) ~ points, data = td)
dca_cox <- dca(Nomo,Score,Stage,model.names = c("Nomogram","RiskScore","N_stage"))
p <-ggplot(dca_cox,
           linetype = F,lwd=1.2)+
  theme_classic()+  
  # theme_prism(base_size =17)+
  theme_bw()+mytheme
scale_x_continuous(
  limits = c(0.25, 1),
  guide = "prism_minor") +
  scale_y_continuous(
    limits = c(-0.01, 0.4),
    guide = "prism_minor")
p
save_plot(p,filename = "返修/TCGA/IndependentFactor/DCA_plot",style = "g",width = 10,height = 6)


library(survivalROC)
sRocFuction=function(td=null,gene=null){
  par(mar= c(5,5,1,1),cex.lab=1.2,cex.axis= 1.2)
  sROC=survivalROC(Stime=td$OS_t, status=td$status, marker = gene, predict.time =1, method="KM")
  plot(sROC$FP, sROC$TP, type="l", xlim=c(0,1), ylim=c(0,1),col="red", 
       xlab="False positive rate", ylab="True positive rate",
       lwd = 2, cex.main=1.3, cex.lab=1.5, cex.axis=1.2, font=1.2)
  abline(0,1)
  aucText=paste0("1 years"," (AUC=",sprintf("%.3f",sROC$AUC),")")
  sROC3=survivalROC(Stime=td$OS_t, status=td$status, marker = gene, predict.time =3, method="KM")
  lines(sROC3$FP, sROC3$TP, type="l", xlim=c(0,1), ylim=c(0,1),col="green",lwd = 2)
  aucText3=paste0("3 years"," (AUC=",sprintf("%.3f",sROC3$AUC),")") 
  sROC5=survivalROC(Stime=td$OS_t, status=td$status, marker = gene, predict.time =5, method="KM")
  lines(sROC5$FP, sROC5$TP, type="l", xlim=c(0,1), ylim=c(0,1),col="blue",lwd = 2)
  aucText5=paste0("5 years"," (AUC=",sprintf("%.3f",sROC5$AUC),")") 
  legend("bottomright", c(aucText,aucText3,aucText5),
         lwd=2,bty="n",col=c("red","green","blue"))
}
pdf("返修/TCGA/IndependentFactor/NomoPoint_TimeROC.pdf",width = 8,height = 6)
sRocFuction(td=td,gene=td$points)
dev.off()
png("返修/TCGA/IndependentFactor/NomoPoint_ROC.png",width = 8,height = 6,units = "in",res = 300)
sRocFuction(td=td,gene=td$points)
dev.off()
write.table(td,file = "返修TCGA/IndependentFactor/Nomopoint_result.txt",quote = F,sep = "\t",row.names = T)

##GEO验证

geo_pha <- read.delim("../../rawdata/GSE20685_clin.txt",row.names = 1)
geo_pha$Status<- ifelse(geo_pha$Status==1,0,1)

geo_pha$M_stage<- str_split_i(geo_pha$M_stage,":",2)
geo_pha$M_stage %>% table()
geo_pha <- geo_pha %>% mutate(N_stage = case_when(M_stage %in% 
                                            c(" 0"," 1") ~ "N0&N1", M_stage %in% 
                                            c(" 2"," 3") ~ "N2&N3",
                                          .default = M_stage))

geo_score<- read.csv("返修/Model/GEO_model_riskscore.csv",row.names = 1)
geo_score$N_stage <- geo_pha$N_stage
identical(rownames(geo_pha),rownames(geo_score))
colnames(geo_score)[c(1,2,9)]<- c("OS","status","RiskScore")
str(geo_score)
geo_score$N_stage <- ifelse(geo_score$N_stage=="N0&N1",1,0)

geo_score$RiskGroup <- ifelse(geo_score$RiskScore> median(geo_score$RiskScore),1,0)

geo_score$OS_t <- geo_score$OS
dd <- datadist(geo_score)
options(datadist="dd")
str(td)
f <- cph(formula(Surv(OS_t, status) ~N_stage+RiskScore), data = geo_score, x = TRUE, y = TRUE,
         surv = TRUE, time.inc = 3)

surv <- Surv(f)
nomo <- nomogram(f)
plot(nomo)
result <- formula_rd(nomo)
result$formula
td$points <- points_cal(formula = result$formula,rd=td)

geo_score$points <- points_cal(formula = result$formula,rd=geo_score)

Score <- coxph(Surv(OS, status) ~ RiskScore, data = xx)
Stage <- coxph(Surv(OS, status) ~ N_stage, data = xx)
Nomo<- coxph(Surv(OS, status) ~ points, data = xx)
dca_cox <- dca(Nomo,Score,Stage,model.names = c("Nomogram","RiskScore","N_stage"))
p <-ggplot(dca_cox,
           linetype = F,lwd=1.2)+
  theme_classic()+  
  # theme_prism(base_size =17)+
  theme_bw()+mytheme
scale_x_continuous(
  limits = c(0.25, 1),
  guide = "prism_minor") +
  scale_y_continuous(
    limits = c(-0.01, 0.4),
    guide = "prism_minor")
p
save_plot(p,filename = "返修/TCGA/IndependentFactor/GEO_DCA_plot",style = "g",width = 10,height = 6)

library(timeROC)
xx<- geo_score
xx$status<- ifelse(xx$status==1,0,1)


library(survivalROC)
sRocFuction=function(td=null,gene=null){
  par(mar= c(5,5,1,1),cex.lab=1.2,cex.axis= 1.2)
  sROC=survivalROC(Stime=td$OS, status=td$status, marker = gene, predict.time =1, method="KM")
  plot(sROC$FP, sROC$TP, type="l", xlim=c(0,1), ylim=c(0,1),col="red", 
       xlab="False positive rate", ylab="True positive rate",
       lwd = 2, cex.main=1.3, cex.lab=1.5, cex.axis=1.2, font=1.2)
  abline(0,1)
  aucText=paste0("1 years"," (AUC=",sprintf("%.3f",sROC$AUC),")")
  sROC3=survivalROC(Stime=td$OS, status=td$status, marker = gene, predict.time =3, method="KM")
  lines(sROC3$FP, sROC3$TP, type="l", xlim=c(0,1), ylim=c(0,1),col="green",lwd = 2)
  aucText3=paste0("3 years"," (AUC=",sprintf("%.3f",sROC3$AUC),")") 
  sROC5=survivalROC(Stime=td$OS, status=td$status, marker = gene, predict.time =5, method="KM")
  lines(sROC5$FP, sROC5$TP, type="l", xlim=c(0,1), ylim=c(0,1),col="blue",lwd = 2)
  aucText5=paste0("5 years"," (AUC=",sprintf("%.3f",sROC5$AUC),")") 
  legend("bottomright", c(aucText,aucText3,aucText5),
         lwd=2,bty="n",col=c("red","green","blue"))
}
pdf("返修/TCGA/IndependentFactor/GEO_NomoPoint_TimeROC.pdf",width = 8,height = 6)
sRocFuction(td=xx,gene=xx$points)
dev.off()
png("返修/TCGA/IndependentFactor/GEO_NomoPoint_ROC.png",width = 8,height = 6,units = "in",res = 300)
sRocFuction(td=xx,gene=xx$points)
dev.off()
write.table(td,file = "返修TCGA/IndependentFactor/Nomopoint_result.txt",quote = F,sep = "\t",row.names = T)

chis
library(rms)
ca1 <- cph(Surv(OS_t,status)~points,data = xx,x=T,y=T,time.inc = 1,surv = T)
call1 <- calibrate(ca1,
                   cmethod = "KM",
                   method = "boot",
                   u=1,
                   m=90,
                   B=1000)
ca2 <- cph(Surv(OS_t,status)~points,data = xx,x=T,y=T,time.inc = 3,surv = T)
call2 <- calibrate(ca2,
                   cmethod = "KM",
                   method = "boot",
                   u=3,
                   m=90,
                   B=1000)

ca3 <- cph(Surv(OS_t,status)~points,data = xx,x=T,y=T,time.inc = 5,surv = T)
call3 <- calibrate(ca3,
                   cmethod = "KM",
                   method = "boot",
                   u=5,
                   m=90,
                   B=1000)
pdf("返修/TCGA/IndependentFactor/GEO_Nomo校准曲线.pdf",width = 8,height = 8)
plot(call1,lwd = 2,lty = 1,errbar.col = c("#EA5455"),
     bty = "l", #只画左边和下边框
     xlim = c(0,1),ylim= c(0,1),
     xlab = "Nomogram-prediced OS (%)",ylab = "Observed OS (%)",
     col = c("#EA5455"),subtitles = F,
     cex.lab=1.2,cex.axis=1, cex.main=1.2, cex.sub=0.6)
lines(call1[,c('mean.predicted',"KM")],
      type = 'b', lwd = 1, col = c("#EA5455"), pch = 16)

plot(call2,lwd = 2,lty = 1,errbar.col = c("#63D07F"),
     xlim = c(0,1),ylim= c(0,1),col = c("#63D07F"),add = T)
lines(call2[,c('mean.predicted',"KM")],
      type = 'b', lwd = 1, col = c("#63D07F"), pch = 16)

plot(call3,lwd = 2,lty = 1,errbar.col = c("#2166AC"),
     xlim = c(0,1),ylim= c(0,1),col = c("#2166AC"),add = T)
lines(call3[,c('mean.predicted',"KM")],
      type = 'b', lwd = 1, col = c("#2166AC"), pch = 16)

abline(0,1, lwd = 2, lty = 3, col = c("#224444"))

legend("bottomright", #图例的位置
       legend = c("1-year","3-year","5-year"), #图例文字
       col =c("#EA5455","#63D07F","#2166AC"), #图例线的颜色，与文字对应
       lwd = 2,#图例中线的粗细
       cex = 1.2,#图例字体大小
       bty = "n")#不显示图例边框
dev.off()


ca1 <- cph(Surv(OS_t,status)~points,data = xx,x=T,y=T,time.inc = 1,surv = T)
call1 <- calibrate(ca1,
                   cmethod = "KM",
                   method = "boot",
                   u=1,
                   m=90,
                   B=1000)
ca2 <- cph(Surv(OS_t,status)~points,data = xx,x=T,y=T,time.inc = 3,surv = T)
call2 <- calibrate(ca2,
                   cmethod = "KM",
                   method = "boot",
                   u=3,
                   m=90,
                   B=1000)

ca3 <- cph(Surv(OS_t,status)~points,data = xx,x=T,y=T,time.inc = 5,surv = T)
call3 <- calibrate(ca3,
                   cmethod = "KM",
                   method = "boot",
                   u=5,
                   m=90,
                   B=1000)
pdf("返修/TCGA/IndependentFactor/GEO_Nomo校准曲线.pdf",width = 8,height = 8)
plot(call1,lwd = 2,lty = 1,errbar.col = c("#EA5455"),
     bty = "l", #只画左边和下边框
     xlim = c(0,1),ylim= c(0,1),
     xlab = "Nomogram-prediced OS (%)",ylab = "Observed OS (%)",
     col = c("#EA5455"),subtitles = F,
     cex.lab=1.2,cex.axis=1, cex.main=1.2, cex.sub=0.6)
lines(call1[,c('mean.predicted',"KM")],
      type = 'b', lwd = 1, col = c("#EA5455"), pch = 16)

plot(call2,lwd = 2,lty = 1,errbar.col = c("#63D07F"),
     xlim = c(0,1),ylim= c(0,1),col = c("#63D07F"),add = T)
lines(call2[,c('mean.predicted',"KM")],
      type = 'b', lwd = 1, col = c("#63D07F"), pch = 16)

plot(call3,lwd = 2,lty = 1,errbar.col = c("#2166AC"),
     xlim = c(0,1),ylim= c(0,1),col = c("#2166AC"),add = T)
lines(call3[,c('mean.predicted',"KM")],
      type = 'b', lwd = 1, col = c("#2166AC"), pch = 16)

abline(0,1, lwd = 2, lty = 3, col = c("#224444"))

legend("bottomright", #图例的位置
       legend = c("1-year","3-year","5-year"), #图例文字
       col =c("#EA5455","#63D07F","#2166AC"), #图例线的颜色，与文字对应
       lwd = 2,#图例中线的粗细
       cex = 1.2,#图例字体大小
       bty = "n")#不显示图例边框
dev.off()


1027/3
ca1 <- cph(Surv(OS_t,status)~points,data = td,x=T,y=T,time.inc = 1,surv = T)
call1 <- calibrate(ca1,
                   cmethod = "KM",
                   method = "boot",
                   u=1,
                   m=290,
                   B=1000)
ca2 <- cph(Surv(OS_t,status)~points,data = td,x=T,y=T,time.inc = 3,surv = T)
call2 <- calibrate(ca2,
                   cmethod = "KM",
                   method = "boot",
                   u=3,
                   m=290,
                   B=1000)

ca3 <- cph(Surv(OS_t,status)~points,data = td,x=T,y=T,time.inc = 5,surv = T)
call3 <- calibrate(ca3,
                   cmethod = "KM",
                   method = "boot",
                   u=5,
                   m=290,
                   B=1000)
pdf("返修/TCGA/IndependentFactor/TCGA_Nomo校准曲线.pdf",width = 8,height = 8)
plot(call1,lwd = 2,lty = 1,errbar.col = c("#EA5455"),
     bty = "l", #只画左边和下边框
     xlim = c(0,1),ylim= c(0,1),
     xlab = "Nomogram-prediced OS (%)",ylab = "Observed OS (%)",
     col = c("#EA5455"),subtitles = F,
     cex.lab=1.2,cex.axis=1, cex.main=1.2, cex.sub=0.6)
lines(call1[,c('mean.predicted',"KM")],
      type = 'b', lwd = 1, col = c("#EA5455"), pch = 16)

plot(call2,lwd = 2,lty = 1,errbar.col = c("#63D07F"),
     xlim = c(0,1),ylim= c(0,1),col = c("#63D07F"),add = T)
lines(call2[,c('mean.predicted',"KM")],
      type = 'b', lwd = 1, col = c("#63D07F"), pch = 16)

plot(call3,lwd = 2,lty = 1,errbar.col = c("#2166AC"),
     xlim = c(0,1),ylim= c(0,1),col = c("#2166AC"),add = T)
lines(call3[,c('mean.predicted',"KM")],
      type = 'b', lwd = 1, col = c("#2166AC"), pch = 16)

abline(0,1, lwd = 2, lty = 3, col = c("#224444"))

legend("bottomright", #图例的位置
       legend = c("1-year","3-year","5-year"), #图例文字
       col =c("#EA5455","#63D07F","#2166AC"), #图例线的颜色，与文字对应
       lwd = 2,#图例中线的粗细
       cex = 1.2,#图例字体大小
       bty = "n")#不显示图例边框
dev.off()


ca1 <- cph(Surv(OS_t,status)~points,data = xx,x=T,y=T,time.inc = 1,surv = T)
call1 <- calibrate(ca1,
                   cmethod = "KM",
                   method = "boot",
                   u=1,
                   m=90,
                   B=1000)
ca2 <- cph(Surv(OS_t,status)~points,data = xx,x=T,y=T,time.inc = 3,surv = T)
call2 <- calibrate(ca2,
                   cmethod = "KM",
                   method = "boot",
                   u=3,
                   m=90,
                   B=1000)

ca3 <- cph(Surv(OS_t,status)~points,data = xx,x=T,y=T,time.inc = 5,surv = T)
call3 <- calibrate(ca3,
                   cmethod = "KM",
                   method = "boot",
                   u=5,
                   m=90,
                   B=1000)
pdf("返修/TCGA/IndependentFactor/GEO_Nomo校准曲线.pdf",width = 8,height = 8)
plot(call1,lwd = 2,lty = 1,errbar.col = c("#EA5455"),
     bty = "l", #只画左边和下边框
     xlim = c(0,1),ylim= c(0,1),
     xlab = "Nomogram-prediced OS (%)",ylab = "Observed OS (%)",
     col = c("#EA5455"),subtitles = F,
     cex.lab=1.2,cex.axis=1, cex.main=1.2, cex.sub=0.6)
lines(call1[,c('mean.predicted',"KM")],
      type = 'b', lwd = 1, col = c("#EA5455"), pch = 16)

plot(call2,lwd = 2,lty = 1,errbar.col = c("#63D07F"),
     xlim = c(0,1),ylim= c(0,1),col = c("#63D07F"),add = T)
lines(call2[,c('mean.predicted',"KM")],
      type = 'b', lwd = 1, col = c("#63D07F"), pch = 16)

plot(call3,lwd = 2,lty = 1,errbar.col = c("#2166AC"),
     xlim = c(0,1),ylim= c(0,1),col = c("#2166AC"),add = T)
lines(call3[,c('mean.predicted',"KM")],
      type = 'b', lwd = 1, col = c("#2166AC"), pch = 16)

abline(0,1, lwd = 2, lty = 3, col = c("#224444"))

legend("bottomright", #图例的位置
       legend = c("1-year","3-year","5-year"), #图例文字
       col =c("#EA5455","#63D07F","#2166AC"), #图例线的颜色，与文字对应
       lwd = 2,#图例中线的粗细
       cex = 1.2,#图例字体大小
       bty = "n")#不显示图例边框
dev.off()
###GSEA功能差异#######
rm(list = ls())
library(GSVA)
library(GSEABase)
library(rstatix)
library(reshape2)
exp <- read.delim("../../rawdata/exprmatx_without_normal.txt",check.names = F)
Group <- clin
exp <- exp[,rownames(Group)]
gmt <- getGmt("../../rawdata/c5.go.v2023.2.Hs.symbols.gmt")
gsva_result <- gsva(as.matrix(exp),gmt,method="gsva",
                    kcdf = "Poisson",parallel.sz=15)

gsva <- t(gsva_result) %>% as.data.frame()
dir.create("返修/TCGA/GSVA")
write.table(gsva,"返修/TCGA/GSVA/gsva_BP_result.txt",sep = "\t",quote = F,row.names = T)
inner <- intersect(rownames(Group),rownames(gsva))
g_result <- gsva[inner,] 
clin_data <- Group[inner,]
identical(rownames(g_result),rownames(clin_data))
g_result$Group <- clin_data$RiskGroup

result <- melt(g_result) %>% group_by(variable) %>% 
  wilcox_test(value~Group) %>% 
  adjust_pvalue(method = "fdr") %>% 
  add_significance("p.adj")
write.table(result,"返修/TCGA/GSVA/gsva_wilcox_test.txt",sep = "\t",quote = F)
sig_gsva <- result %>% filter(p.adj.signif != "ns")
top30<- arrange(sig_gsva,p.adj)[1:30,]
write.table(top30,"返修/TCGA/GSVA/top30_gsva_result.txt",quote = F,sep = "\t",row.names = F)
top30<- read.delim("返修/TCGA/GSVA/top30_gsva_result.txt")
top_gsva <- g_result[,top30$variable]

library(ComplexHeatmap)
anncol <- HeatmapAnnotation(Group = Group$RiskGroup,
                            col = list(Group = c("Low"= "#1A9641", "High"="#D7191C")),
                            annotation_legend_param = list(Group = list(nrow=1)))
plot_data <- t(top_gsva)
rownames(plot_data) <- str_remove(rownames(plot_data),pattern = "GOBP_")
p <-Heatmap(plot_data,na_col = "white",
            show_column_names = F,
            show_row_names = T,name = "GSVA Score",row_names_gp = grid::gpar(fontsize=11),
            column_order = c(colnames(plot_data)[c(grep("High",Group$RiskGroup),grep("Low",Group$RiskGroup))]),
            column_split = Group$Group,
            cluster_columns = F,column_title = NULL,
            top_annotation  = anncol,
            heatmap_legend_param = list(direction="horizontal"),
            column_names_rot = 90,row_names_max_width = unit(20,"cm"),
)
pdf("返修/TCGA/GSVA/gsva_top30_heatmap.pdf",width = 15,height = 10)
draw(p,heatmap_legend_side="bottom",
     annotation_legend_side="bottom",merge_legend=T)
dev.off()
png("返修/TCGA/GSVA/gsva_top30_heatmap.png",width = 15,height = 10,units = "in",res = 300)
draw(p,heatmap_legend_side="bottom",
     annotation_legend_side="bottom",merge_legend=T)
dev.off()

###免疫检查点#########
rm(list = ls())
IGG <- readRDS("../../rawdata/IGG.Rds")
exp <- read.delim("../../rawdata/exprmatx_without_normal.txt",row.names = 1,check.names = F)
Rs_g <- clin
inner <- intersect(rownames(Rs_g),colnames(exp))
exp1 <- exp[,inner]
Rs_g1 <- Rs_g[inner,]
plot_data <- exp1[intersect(IGG,rownames(exp1)),]

data <- t(plot_data) %>% as.data.frame()
identical(rownames(data),rownames(Rs_g1))
data$Group<- Rs_g1$RiskGroup
Cluster_data <- melt(data)
result <- melt(data) %>% group_by(variable) %>% 
  wilcox_test(value~Group) %>% 
  adjust_pvalue(method = "fdr") %>% 
  add_significance("p.adj")
dir.create("返修/TCGA//IGGs")
write.table(result,file = "返修/TCGA//IGGs/Iggs_wilcox_test.txt",row.names = F,sep = "\t",quote = F)
sig_ICC <- result %>% filter(p.adj.signif != "ns")
plot_data <- Cluster_data %>% filter(variable %in% sig_ICC$variable)
p <- ggplot(plot_data,aes(variable,value,fill=Group))+
  geom_boxplot(color="black",alpha=0.6,outlier.shape = 21,outlier.size = 1.2)+
  scale_fill_manual(values = c("#D7191C","#1A9641"))+
  labs(x="IGGs",y="Expression level")+
  stat_compare_means(aes(group= Group),label = "p.signif",size=5)+
  theme_bw()+mytheme+theme(axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1))
p

dir.create("TCGA/Lasso_cox/IGG")

save_plot(p,filename = "返修/TCGA/IGGs/Risk_group_IGG_plot",width =20 ,height = 6,style = "g")

###TME ssGSEA###########
library(ggpubr)
library(lemon)
rm(list = ls())
tme <- read.delim("../乳腺癌前哨淋巴结mRNA+lncRNA/P1/mRNA/TCGA/ConsensusCluster/TME_Score/Estamite_result.txt",row.names = 1)
RS_group <-clin
tme1 <- tme[rownames(RS_group),]
tme1 <- tme1[,-5]

identical(rownames(tme1),rownames(RS_group))
tme1$RS_group <- RS_group$RiskGroup
plot_data <- melt(tme1)
dir.create("返修/TCGA/TME")
p <-ggplot(plot_data,aes(RS_group,value,fill=RS_group))+
  geom_violin(trim = F,aes(fill=RS_group),color="white",cex=1,alpha=0.5)+
  geom_boxplot(width=0.1,position = position_dodge(0.9),color="gray")+
  # geom_signif(comparisons = list(c("High","Low")),step_increase = 0.1,size=1,test = "wilcox.test",textsize = 5)+
  facet_rep_wrap(.~variable,scales = 'free',repeat.tick.labels = 'left',ncol = 4)+
  stat_compare_means(comparisons = combn(unique(plot_data$RS_group), 2, simplify =FALSE),
                     method = 't.test',size=6)+
  theme_bw()+
  scale_fill_manual(values = c("#d9352a","#4979b6"))+ylab("TME Score")+xlab("")+
  theme(panel.grid=element_blank())+mytheme

p
write.table(tme1,"返修/TCGA/TME/RS_TME_result.txt",quote = F,sep = "\t",row.names = T)
save_plot(p,filename = "返修/TCGA/TME/RS_group_estimate_boxplot",style = "g",width = 12,height = 8)
##ssGSEA###
g_result <- read.delim("../乳腺癌前哨淋巴结mRNA+lncRNA/P1/mRNA/TCGA/ConsensusCluster/GSVA/ssGSEA_result,txt",check.names = F,row.names = 1)
g_result <- select(g_result,-RSGroup,-Cluster)
clin_data <- clin
inner <- intersect(rownames(clin_data),rownames(g_result))
g_result <- g_result[inner,] 
clin_data <- clin_data[inner,]
identical(rownames(g_result),rownames(clin_data))
g_result$RS_group <- clin_data$RiskGroup

result <- melt(g_result) %>% group_by(variable) %>% 
  wilcox_test(value~RS_group) %>% 
  adjust_pvalue(method = "fdr") %>% 
  add_significance("p.adj")
write.table(result,"返修/TCGA/TME/ssgsea_wilcoxTest.txt",row.names = F,sep = "\t",quote = F)
write.table(g_result,"返修/TCGA/TME/ssGSEA_result.txt",row.names = T,quote = F,sep = "\t")
plot_data <- melt(g_result)
data_arrange <- plot_data %>% group_by(variable) %>% 
  summarise(de = median(value)) %>% 
  arrange(desc(de)) %>% pull("variable")

plot_data$variable <- factor(plot_data$variable,levels = unique(data_arrange))
p <- ggplot(plot_data,aes(reorder(variable,-value),value,fill=RS_group))+
  geom_boxplot(color="black",alpha=0.6,outlier.shape = 21,outlier.size = 1.2)+
  scale_fill_manual(values = c("#E7298A","#1B9E77"))+
  labs(x="Cell type",y="ssGSEA Score")+
  stat_compare_means(aes(group= RS_group),label = "p.signif",size=5)+
 theme_bw()+mytheme+
  theme(axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1))
p
save_plot(p ,filename = "返修/TCGA/TME/ssgsea_RSGroup_boxplot",style = "g",width = 10,height = 6)


##相关性
library(IOBR)
exp <- read.delim("../../rawdata/exprmatx_without_normal.txt",row.names = 1,check.names = F)
cibersort <- deconvo_tme(eset = exp, 
                         method = "cibersort", 
                         arrays = FALSE, 
                         perm = 500)
g_result<- cibersort %>% filter(`P-value_CIBERSORT`<0.05)
g_result<- g_result[,-c(24,25,26)] %>% as.data.frame()
colnames(g_result)<- str_remove(colnames(g_result),"_CIBERSORT")
g_result<- g_result %>% column_to_rownames("ID")

g_result <- read.delim("../乳腺癌前哨淋巴结mRNA+lncRNA/P1/mRNA/TCGA/ConsensusCluster/GSVA/ssGSEA_result,txt",check.names = F,row.names = 1)
g_result <- dplyr::select(g_result,-RSGroup,-Cluster)

risk<- fread("返修/Model/TCGA_model_riskscore.csv",data.table = F)
risk<- risk[,c(1,4)]
risk<- risk %>% column_to_rownames("ID")
inn<- intersect(rownames(risk),rownames(g_result))
risk<- risk[inn,]
g_result<- g_result[inn,]
g_result$RiskScore<- risk
library(psych)
library(corrplot)
res<- corr.test(g_result,method = "pearson")

R<- res$r
P<- res$p
dir.create("返修2/风险模型免疫相关性")
fwrite(as.data.frame(R),"返修2/风险模型免疫相关性/corr_r.csv",row.names = T)
fwrite(as.data.frame(P),"返修2/风险模型免疫相关性/corr_p.csv",row.names = T)


corr <- as.matrix(R)
pvalue <- as.matrix(P)

addcol <- colorRampPalette(c("blue", "white", "red"))
pdf("返修2/风险模型免疫相关性/corr_heatmap.pdf", width = 10, height = 10) 
png("返修2/风险模型免疫相关性/corr_heatmap.png", width = 10, height = 10,units = "in",res = 300) 
corrplot(corr, # 相关性矩阵
         method = "color", # 表示用颜色表示相关性大小
         col = addcol(100), 
         tl.col = "black", # 文本标签的颜色
         tl.cex = 1, # 文本标签的字符大小
         tl.srt = 90, #  文本标签的旋转角度
          tl.pos = "td", # 文本标签位置，td表示顶部和对角线 
         p.mat = pvalue, #  P 值矩阵
         diag = F, # 是否显示对角线上的相关性值
         type = 'upper', # 只绘制上三角部分
         sig.level = c(0.05), # 设置显著性水平阈值，可设置多个
         pch.cex = 1.5,  # 显著性标记字符大小
         pch.col = 'grey20',  # 显著性标记字符颜色
         insig = 'label_sig',
         order = 'AOE', #设置一种排序方式
)
dev.off()

###药物敏感性####
rm(list = ls())
Drug <- read.delim("../../rawdata/calcPhenotype_Output/calcPhenotype_Output/BRCA_DrugPredictions.csv",sep = ",",row.names = 1)
Cluster <- clin
inner <- intersect(rownames(Cluster),rownames(Drug))
Drug1 <-Drug[inner,] 
identical(rownames(Drug1),rownames(Cluster))

data <- Drug1
identical(rownames(data),rownames(Cluster))
data$Group <- Cluster$RiskGroup
group_data <- melt(data)
result <- melt(data) %>% group_by(variable) %>% 
  wilcox_test(value~Group) %>% 
  adjust_pvalue(method = "fdr") %>% 
  add_significance("p.adj")
dir.create("返修/TCGA/Drug")
write.table(result,"返修/TCGA/Drug/Drug_wilcox_test.txt",row.names = F,sep = "\t",quote = F)
sig_Drug <- result %>% filter(p.adj.signif != "ns")
plot_data <- group_data %>% filter(variable %in% sig_Drug$variable)
p <- ggplot(plot_data,aes(variable,value,fill=Group))+
  geom_boxplot(width=1,position = position_dodge(0.9),color="black",outlier.shape = 21)+
  scale_fill_manual(values = c("#E7298A","#CAB2D6"))+
  labs(x="",y="Drug response level")+
  stat_compare_means(aes(group= Group),label = "p.signif",size=5)+
  theme_bw()+mytheme+theme(axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1))+
  geom_hline(yintercept = 0,lty=4,col="black",lwd=0.8)

p
save_plot(p ,filename = "返修/TCGA//Drug/Drug_Sign_boxplot",style = "g",width = 25,height = 6)

plot_data <- data[,sig_Drug$variable] %>% t()
exp <- apply(plot_data, 1, scale)
rownames(exp) <- colnames(plot_data)
exp <- t(exp)
anncol <- HeatmapAnnotation(RsGroup = Cluster$RiskGroup,
                            col = list(RsGroup = c("Low"= "#1A9641", "High"="#D7191C")))
p <-Heatmap(exp,na_col = "white",
            show_column_names = F,
            # width = ncol(mat)*unit(5, "mm"), 
            height = nrow(exp)*unit(4, "mm"),
            show_row_names = T,name = "Drug response level",row_names_gp = grid::gpar(fontsize=11),
            column_order = c(colnames(exp)[c(grep("High",Cluster$RiskGroup),grep("Low",Cluster$RiskGroup))]),
            column_split = Cluster$Group,
            cluster_columns = F,column_title = NULL,
            top_annotation = anncol,)
pdf("返修/TCGA/Drug/SigDrug_heatmap.pdf",width = 10,height = 20)
draw(p,heatmap_legend_side="right",
     annotation_legend_side="right",merge_legend=T)
dev.off()
png("返修/TCGA/Drug/SigDrug_heatmap.png",width = 10,height = 20,units = "in",res = 300)
draw(p,heatmap_legend_side="right",
     annotation_legend_side="right",merge_legend=T)
dev.off()

data_arrange <- sig_Drug %>% group_by(variable) %>% 
  arrange(desc(p.adj)) %>% pull(variable)
sig_Drug <- sig_Drug %>% arrange(p.adj)
Top10 <- sig_Drug[1:10,]
plot_data <- group_data %>% filter(variable %in% Top10$variable)
p <- ggplot(plot_data,aes(variable,value,fill=Group))+
  geom_violin(trim = F,aes(fill=Group),color="white",alpha=0.5)+
  geom_boxplot(width=0.2,position = position_dodge(0.9),color="gray",outlier.shape = 21)+
  scale_fill_manual(values = c("#E7298A","#1A9641"))+
  labs(x="",y="Drug response level")+
  stat_compare_means(aes(group= Group),label = "p.signif",size=5)+
  theme_bw()+mytheme+theme(axis.text.x = element_text(angle = 45,hjust  = 0.9,vjust = 0.9))+
  geom_hline(yintercept = 0,lty=4,col="black",lwd=0.8)
p
save_plot(p ,filename = "返修/TCGA//Drug/Drug_Sign10_boxplot",style = "g",width = 8,height = 6)

###HPA下载#########
library(BiocStyle)
install.packages("HPAanalyze")
BiocManager::install("HPAanalyze")
library(HPAanalyze)
library(dplyr)
library(tibble)
library(readr)
library(tidyr)
dir.create("返修/TCGA/HPA")
genes <- fread("返修/Model/Gene_Univar_Cox_result.csv") %>% filter(pvalue<0.05)
genes <- genes$id
tissue="Breast"
filgene <- genes
gene <-"AREG"
for (gene in filgene) {
  #获得HPA网站中该基因的xml文件
  hpa_target_gene<-try(hpaXmlGet(gene))
  if(class(hpa_target_gene)[1] == "try-error"){
    print(paste0(gene,"不在数据库！"))
    next
  }else{
    #将xml中组织染色的信息提取出来
    
    hpa_target_gene_fig_url<-hpaXmlTissueExpr(hpa_target_gene)
    hpa_target_gene_fig_url_1<-try(as.data.frame(hpa_target_gene_fig_url[[1]]))
    if(class(hpa_target_gene_fig_url_1) == "try-error"){
      print(paste0(gene," 下载不了"))
      next
    }else{
      tmp <- try(hpa_target_gene_fig_url_1[1:6,1:18])
      if(class(tmp) == "try-error"){
        print(paste0(gene," 下载不了"))
        next
      }else{
        #选择自己感兴趣的组织
        hpa_target_gene_fig_url_tissue<-hpa_target_gene_fig_url_1[hpa_target_gene_fig_url_1$tissueDescription2==tissue,]
        # hpa_target_gene_fig_url_tissue<-hpa_target_gene_fig_url_2[hpa_target_gene_fig_url_2$tissueDescription2==tissue,]
        
        # hpa_target_gene_fig_url_2<-as.data.frame(hpa_target_gene_fig_url[[2]])
        # hpa_target_gene_fig_url_2[1:6,1:18]#为该组织该基因单独建个文件夹储存
        picDir <- paste('返修/TCGA/HPA/',gene, tissue,"IHC-2/", sep = "_")
        if (!dir.exists(picDir)) {
          dir.create(picDir)
        }
        
        
        for (i in 1:nrow(hpa_target_gene_fig_url_tissue)) {
          file_url<-hpa_target_gene_fig_url_tissue$imageUrl[i]
          file_dir<-paste(picDir,gene,tissue,hpa_target_gene_fig_url_tissue$patientId[i],hpa_target_gene_fig_url_tissue$tissueDescription1[i],hpa_target_gene_fig_url_tissue$tissueDescription2[i],".tiff",sep = "_")
          print(paste0("正在下载：",gene))
          download.file(url = file_url,destfile = file_dir,mode = "wb")
        }
      }
    }
  }
  
}

library(KEGGREST)
pathways <- keggList("pathway","hsa")
one.pathway <- keggGet("hsa05200")[[1]]
one_pathway$GENE
for(pathway in names(pathways)){
  one.pathway <- keggGet(pathway)[[1]]
  if(!is.null(one.pathway$GENE)){
    pathway.genes <- one.pathway$GENE
    indices <- seq(2,length(pathway.genes),2)
    indices_de <- NULL
    for(i in indices){
      if (grepl(";", pathway.genes[i])) {
        next
      }else{
        indices_de <- c(indices_de,i)
      }
    }
    if(!is.null(indices_de)){
      indices_de <- c(indices_de,indices_de-1)
      pathway.genes <- pathway.genes[-c(indices_de)]
    }
    genes.info <- unlist(lapply(pathway.genes,function(x) strsplit(x,";")))
    genes <- genes.info[1:length(genes.info)%%3 ==2]
    genes.reshape <- paste(genes,collapse = ",")
    result <- paste(c(one.pathway$ENTRY,one.pathway$PATHWAY_MAP,genes.reshape),collapse = "\t")
    write.table(result, file = 'hsa_kegg_pathway.txt', sep = "\t",row.names = F,col.names=F, append=T,quote=FALSE)
  }
}
