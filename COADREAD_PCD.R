dir.create("D:/WORK/Project/结直肠癌_程序性死亡_模型25")
setwd("D:/WORK/Project/结直肠癌_程序性死亡_模型25/")
dir.create("data")
library(data.table)
library(tidyverse)
library(CyDataPro)
PCD <- read.csv("data/PCD_genelist.csv")
load("../结直肠癌_RNAmSNP_MR_Singlecell16/data/TCGA_CRC.Rdata")
load("data/TCGA_Tumor_surv.Rdata")
PCD <- fread("data/PCD_genelist.csv",check.names = T) %>% as.data.frame()

Sample <- TCGA_Tumor_Surv %>% filter(OS_time > 0)
exp <- dectect_log(count_matix)
inn <- intersect(colnames(exp),rownames(Sample))
exp <- exp[,inn]
Sample <- Sample[inn,]
meta<- Sample %>% select(OS_time,Event)
meta$OS_time <- as.numeric(meta$OS_time)
pFilter <- 0.05 
PCDFliter <- list()
SigGene = vector()
i=colnames(PCD)[1]
for(i in colnames(PCD)){
  SigGene = vector()
  tmpE<- exp[intersect(rownames(exp),PCD[,i]),] %>% t() %>% as.data.frame()
  svdata <- cbind(meta,tmpE)
  for(j in colnames(svdata[3:ncol(svdata)])){
    cox <- coxph(Surv(OS_time,Event) ~ svdata[,j], data = svdata)
    coxSummary <- summary(cox)
    uniCox <- rbind(uniCox, data.frame(gene = j,
                                       HR = coxSummary$conf.int[,"exp(coef)"],
                                       HR.95L = coxSummary$conf.int[,"lower .95"],
                                       HR.95H = coxSummary$conf.int[,"upper .95"],
                                       P_value = coxSummary$coefficients[,"Pr(>|z|)"]))
    coxP <- coxSummary$coefficients[,"Pr(>|z|)"]
    if(coxP<pFilter){
      SigGene <- c(SigGene,j)
    }
  }
  xx <- list(SigGene)
  names(xx)<-i
  PCDFliter <- append(PCDFliter,xx)
}
xx <- data.frame(variables = names(PCDFliter),
                 SigGene_N = lengths(PCDFliter))

fwrite(xx,file = "PCDFilter_result.csv")
save(PCDFliter,file = "data/PCD_COX_filter.Rdata")
# add prognostic label

##ssGSEA
library(GSVA)
exp_log <- dectect_log(count_matix)
Tumor <- exp_log[,as.numeric(str_sub(colnames(exp_log),14,15)) <10]
PCD <- PCDFliter[lengths(PCDFliter) >=5]
xx <- as.data.frame(t(sapply(PCD, "[", i = 1:max(sapply(PCD, length)))))
xx <- t(xx)
fwrite(xx,"PCDFilter_genelist.csv")
genelist <- as.list(PCD)
dir.create("ssGSEA")
gsea<- GSVA::gsva(as.matrix(Tumor),genelist,method="ssgsea")
ES <- as.data.frame(t(gsea))
ES_zscore <- scale(ES) %>% data.frame() 
pheatmap::pheatmap(t(ES_zscore),scale = "none",show_colnames = F,
                   fontsize = 15,main = "ssGSEA Score",
                   color = colorRampPalette(c("#6A3D9A","white","#E31A1C"))(1000),
                   filename = "ssGSEA/ssGSEA_result_heatmap.pdf",
                   width = 10,height = 6)
fwrite(ES,file = "ssGSEA/ssGSEA_result.csv",row.names = T)
fwrite(ES_zscore,file ="ssGSEA/ssGSEA_result_zscore.csv")
library(rstatix)
library(Hmisc)
library(corrplot)
cor <- rcorr(as.matrix(ES),type = "spearman")
R <- cor$r
P <- cor$P
P[is.na(P)] =0
pdf("ssGSEA/ssGSEA_corr_heatmap.pdf",width = 8,height = 6)
corrplot(as.matrix(R),
            p.mat = as.matrix(P),
            insig = "label_sig",
             sig.level = .05,
            pch.cex = 1,
            method = "square",
            order = "alphabet",
            type = "lower",
            tl.col = "black",
            tl.pos = "ld",
            tl.srt=45,
            is.corr=T,
            col = rev(COL2("PiYG",10)),
            cl.pos = "b")

dev.off()
dev.new()
fwrite(as.data.frame(R),file = "ssGSEA/ssGSEA_corr_spearman_R.csv",row.names = T)
fwrite(as.data.frame(R),file = "ssGSEA/ssGSEA_corr_spearman_P.csv",row.names = T)

###ssGAES Score与肿瘤不同阶段
dir.create("ssGSEA/clinCorr")
library(ComplexHeatmap)
clin <- as.data.frame(clin)
table(tmp1$Event)
rownames(clin) <- clin$Sample_id
table(clin$OS_time)
tmp <- clin %>% filter(OS_time !="[Discrepancy]" & OS_time !="")
tmp$Stage <- toupper(tmp$Stage)
tmp <- tmp %>% mutate(Stage = case_when(Stage %in% c("STAGE I", 
                                                     "STAGE IA", "STAGE IB", "STAGE IC") ~ "I", Stage %in% 
                                          c("STAGE II", "STAGE IIA", "STAGE IIB", "STAGE IIC") ~ 
                                          "II", Stage %in% c("STAGE III", "STAGE IIIA", "STAGE IIIB", 
                                                             "STAGE IIIC") ~ "III", Stage %in% c("STAGE IV", 
                                                                                                 "STAGE IVA", "STAGE IVB", "STAGE IVC") ~ "IV", Stage %in% 
                                          c("STAGE X", "NOT REPORTED","[DISCREPANCY]",'') ~ "X", .default = Stage))
tmp <- tmp %>% filter(Stage !="X")
tmp <- tmp %>% filter(T_stage !="Tis")
tmp <- tmp %>% filter(M_stage !="" & M_stage !="MX")
table(tmp$M_stage)
TCGA_Tumor_Surv <- tmp
save(TCGA_Tumor_Surv,file = "data/TCGA_Tumor_surv.Rdata")
inn<- intersect(rownames(ES),rownames(tmp))
ES1 <- ES_zscore[inn,]
tmp1 <- tmp[inn,] 
identical(rownames(ES1),rownames(tmp1))
i=colnames(ES1)[1]
for(i in colnames(ES1)){
  meta <- tmp1
  meta[,i] <- ES1[,i]
  meta$OS_year <- as.numeric(meta$OS_time)/365 %>% round()
  data <- meta[,c(6,7,13:ncol(meta))]
  data<-data[order(data[,i]),]
  data%>%glimpse()
  
  survdata<-data.frame(row.names = rownames(data),OS=data[,ncol(data)])
  data <- data[,-ncol(data)]
  data$Event <- factor(data$Event,levels = c("1","0"),labels = c("Alive","Death"))
  colnames(data)[8] <- "ssGSEA Score"
  
  ha=HeatmapAnnotation(survival_time = anno_points(survdata,size =unit(0.1,'cm'),#点注释
                                                   gp=gpar(col="grey")),df=data,#点注释属性及其他临床特征注释
                       col = list(Age=c("<65"="#66C2A5",">=65"="#FC8D62"),#设置部分临床特征注释颜色
                                  Event=c("Alive"="#A6D854","Death"="#E78AC3"),
                                  Stage = c("I"="#E64B35B2","II"="#4DBBD5B2","III"="#3C5488B2","IV"="#00A087B2"),
                                  Gender = c("FEMALE" = "#BC3C29FF","MALE" = "#0072B5FF"),
                                  `ssGSEA Score`=circlize::colorRamp2(c(-2,0,3), c("blue", "white", "red"))))
  heat<-Heatmap(matrix(nrow=0, ncol=nrow(data)),top_annotation = ha,
                column_title = i,
                column_title_gp = gpar(fontsize = 20,
                                       fontface = 'bold', 
                                       # fill = 'red',
                                       col = 'black'),
                rect_gp = gpar(col= "#FFFFFF00"))
  
  pdf(paste0("ssGSEA/clinCorr/",i,"_clin_heatmap.pdf"),width = 8,height = 6)
  draw(heat,annotation_legend_side="bottom")
  dev.off()
}
i <- colnames()
library(ggpubr)
ES1 <- ES[inn,]
clin <- clins[1]
clins <- colnames(tmp1)[c(6,7,14,15,16,17)]
chisq_result <- data.frame()
for(i in colnames(ES1)){
  tmp1$Group <- ifelse(as.numeric(ES1[,i]) > median(as.numeric(ES1[,i])),"High","Low")
  for(clin in clins){
    tab <- xtabs(~tmp1$Group+tmp1[,clin])
    tt<-chisq.test(tab)
    tmp <- c(i,clin,tt$p.value)
    chisq_result <- rbind(chisq_result,tmp)
    P <- ifelse(tt$p.value < 0.001,"< 0.001",round(tt$p.value,3))
    p<-ggplot(tmp1, aes(x = Group, fill = tmp1[,clin])) +
      geom_bar(width = 0.9, position = "fill") + # 百分比柱状图
      scale_fill_brewer(palette = "PuBu")  +
       scale_y_continuous(labels = scales::percent) +
      labs(title = paste0(i," and ",clin),
           y = "Percentage",fill="")+
      theme(panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),
            panel.background = element_blank(),
            axis.text.x = element_text(face = "bold",colour = "black",size=15),
            axis.text.y = element_text(face = "bold",colour = "black",size=15),
            axis.ticks.x = element_blank(),
            axis.title.y = element_blank(),
            axis.title.x = element_blank(),
            plot.title = element_text(hjust = 0.5, face="bold"))+
      coord_cartesian(clip = 'off', ylim = c(0,1)) + 
      theme(plot.margin = margin(0.5,0.5,1.4,0.5,'cm'))+ #自定义图片上左下右的边框宽度
      annotate( "text",
                cex=6,
                x=1.5, y=-0.15, # 根据自己的数据调节p value的位置
                label=paste0("Pval= ",P), # 添加P值
                color="black"
      )+
      annotate("rect", xmin = 0.55, xmax = 1.45, ymin = -0.1, ymax = -0.05, 
               fill = "#40a1cf")+
      annotate("rect", xmin = 1.55, xmax = 2.45, ymin = -0.1, ymax = -0.05, 
               fill = "#dd816d")
    p
    dir.create(paste0("ssGSEA//clinCorr/batplot/",i),recursive = T)
    if(P < 0.05){
      ggsave(filename = paste0("ssGSEA//clinCorr/batplot/",i,"/",clin,"_Diff_barplot.pdf"),plot = p,width = 8,height = 6)
    }else{
      ggsave(filename = paste0("ssGSEA//clinCorr/batplot/",i,"/",clin,"_NoDiff_barplot.pdf"),plot = p,width = 8,height = 6)
    }
    
  }
  
}
-log10(0.05)
colnames(chisq_result)<- c("PCD","Clin","Pvalue")
min(chisq_result$Pvalue)
xx <- chisq_result %>% filter(Pvalue < 0.05)
chisq_result$logP <- -log10(as.numeric(chisq_result$Pvalue))
plot <- chisq_result %>% group_by(Clin) %>% arrange(Pvalue)
plot$PCD <- factor(plot$PCD,levels = unique(plot$PCD))
plot$Clin <- factor(plot$Clin,levels = unique(plot$Clin))
p<-ggplot(plot,aes(PCD,Clin))+
  geom_point(aes(size=logP,color=logP))+ #以-Log10(FDR)代表点的大小，Log2FC代表点的颜色
  scale_color_gradient2(low="#3F8CD1",mid = "#F5F5F5",high = "#F08080",midpoint = 1.31)+#设置颜色值
  scale_size_continuous(range = c(2,8))+ #调整点的大小为4-10
  theme_bw()+ #删除背景灰色
  theme(legend.position = "bottom")+ #修改图例的位置
  labs(color = "-log10(Pvalue)", size = "-log10(Pvalue)")+ # 修改图例的名字
  xlab(NULL)+ #删除横坐标
  ylab(NULL)+ #删除纵坐标
  mytheme+theme(axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1))
p
save_plot(p,filename = "ssGSEA/clinCorr/chiqs_test_pointPlot",style = "g",width = 10)

fwrite(chisq_result,file = "ssGSEA/clinCorr/chiqs_test_result.csv")

##LAsso
library(glmnet)
library(survival)
library(broom)
library(RColorBrewer)
identical(rownames(tmp1),rownames(ES1))
svdata <- cbind(tmp1,ES1)
svdata$OS_time <- as.numeric(svdata$OS_time)
svdata <- svdata %>% filter(OS_time > 0)
clint <- colnames(ES1)
x <- svdata[,clint] %>% as.matrix()
class(svdata$OS_time)
which(meta2$OS_time==0)
y<- data.matrix(Surv(time = svdata$OS_time,event = svdata$Event))
fit_lasso <- glmnet(x,y,alpha = 1,family = "cox")
pdf("4/Model/lambda_fit.pdf",width = 8,height = 6)
plot(fit_lasso,xvar="norm",label=T)
dev.off()
#交叉验证
cv_lasso <- cv.glmnet(x,y,family="cox",type.measure = "deviance")
pdf("4/Model/lambda_cv.pdf",width = 8,height = 6)
plot(cv_lasso)
abline(x=min)
dev.off()
min = cv_lasso$lambda.min#选择lambda min  0.01544244
se <- cv_lasso$lambda.1se
tidy_df <- broom::tidy(fit_lasso)
tidy_cvdf <- broom::tidy(cv_lasso)
mypalette <- c(brewer.pal(11,"BrBG"),brewer.pal(11,"Spectral"),brewer.pal(5,"Accent"))
p<-ggplot(tidy_df, aes(lambda, estimate, group = term, color = term)) +
  geom_line(size=1.2)+
  geom_hline(yintercept = 0)+
  geom_vline(xintercept = c(min,se),color=c("red","blue"),lty=4,lwd=1)+
  scale_x_log10(name = "Log Lambda")+
  ylab("Coefficients")+
  scale_color_manual(name="variable",values = mypalette)+
  # annotate(geom = "text", x = 0.01, y = 0.3, label = paste0("min Lambda",round(min,4)), size = 5)+
  theme_bw()+mytheme
p
dir.create("ssGSEA/Model/Lasso",recursive = T)
save_plot(p,filename = "ssGSEA/Model/Lasso/lambda_fit",style = "g")
p<-ggplot()+
  geom_point(data=tidy_cvdf, aes(lambda,estimate),color="red")+
  geom_errorbar(data = tidy_cvdf, aes(x=lambda,ymin=conf.low,ymax=conf.high))+
  geom_vline(xintercept = c(min,se),color=c("red","blue"),lty=4,lwd=1)+
  scale_x_log10(name = "Log Lambda")+
  ylab("Coefficients")+
  theme_bw()+mytheme
p
save_plot(p,filename = "ssGSEA/Model/Lasso/lambda_cv",style = "g")
model_lambda_min <- glmnet(x,y,alpha = 1,family = "cox",lambda = min)
best <- broom::tidy(model_lambda_min)
best_gene <- model_lambda_min$beta

gene_min<-rownames(best_gene)[as.numeric(best_gene)!=0]#lasso 基因
lambda_gene <- gene_min %>% as.data.frame()
colnames(lambda_gene)<- "Geneid"
fwrite(best,file = "ssGSEA/Model/Lasso/Lasso_optVar.csv")
lam <- data.frame(min=min,se=se)
fwrite(lam,file = "ssGSEA/Model/Lasso/lambda_value.csv")
intersect(best$term,lambda_gene$Geneid)

##单变量COX
library(survival)
dir.create("ssGSEA/Model/COX")
svdata <- cbind(tmp1,ES1)
svdata$OS_time <- as.numeric(svdata$OS_time)
svdata <- svdata %>% filter(OS_time > 0)
clint <- colnames(ES1)
pFilter <- 0.05 
uniCox <- data.frame()
uniCoxSig <- data.frame()
for(j in colnames(svdata[,clint])){
  cox <- coxph(Surv(OS_time,Event) ~ svdata[,j], data = svdata)
  coxSummary <- summary(cox)
  uniCox <- rbind(uniCox, data.frame(gene = j,
                                     HR = coxSummary$conf.int[,"exp(coef)"],
                                     HR.95L = coxSummary$conf.int[,"lower .95"],
                                     HR.95H = coxSummary$conf.int[,"upper .95"],
                                     P_value = coxSummary$coefficients[,"Pr(>|z|)"]))
  coxP <- coxSummary$coefficients[,"Pr(>|z|)"]
  if(coxP<pFilter){
    uniCoxSig <- rbind(uniCoxSig, data.frame(variable = j, 
                                             HR = coxSummary$conf.int[,"exp(coef)"],
                                             HR.95L = coxSummary$conf.int[,"lower .95"],
                                             HR.95H = coxSummary$conf.int[,"upper .95"],
                                             P_value = coxSummary$coefficients[,"Pr(>|z|)"]))
  }
}

# add prognostic label
uniCox$Risk_group <- ifelse(uniCox$HR >= 1 & uniCox$P_value < 0.05, "Risk",
                            ifelse(uniCox$HR < 1 & uniCox$P_value < 0.05, "Protect","Not_sig"))
uniCoxSig$Risk_group <- ifelse(uniCoxSig$HR >= 1, "Risk", "Protect")
uniCox <- arrange(uniCox,HR)
uniCox$gene<- factor(uniCox$gene,levels = uniCox$gene)
xx <- scale(uniCox$HR)
table(svdata$Event)
p<-ggplot(uniCox, aes(HR, gene)) + 
  geom_errorbarh(aes(xmax = HR.95L, xmin = HR.95H,col= Risk_group),height=0.3,size=0.7) +
  geom_point(shape=18, aes(col=Risk_group, size = -log10(P_value)))  +
  geom_vline(aes(xintercept = 1),color="black",linetype = "dashed",size=1) + 
  xlim(0,15000)+
  labs(title = "uniCox",y="")+
  scale_color_manual(values = c("gray","red","red"))+
  theme_bw()+mytheme+theme(plot.title = element_text(hjust = 0.5))
p
CyDataPro::save_plot(p,filename = "ssGSEA/Model/COX/uniCox_foresplot",style = "g")
fwrite(uniCox,file = "ssGSEA/Model/COX/uniCOX_result.csv")

Iterm <- data.frame(Variables = intersect(uniCoxSig$variable,best$term))
fwrite(Iterm,file = "ssGSEA/Model/Lasso_uniCOX_Sig_result.csv")
save(svdata,file = "data/TCGA_surv.Rdata")
library(ggvenn)
library(Cycolors)
mycol <- Disc5(2)
data2 <-list(Lasso = best$term,uniCOX = uniCoxSig$variable)

p1<-ggvenn(data2,
           columns = NULL,
           show_elements = F,
           label_sep = "\n",
           show_percentage = T,
           digits = 1,
           fill_color = mycol,
           fill_alpha = 0.5,
           stroke_color = "white",
           stroke_alpha = 0.5,
           stroke_size = 0.5,
           stroke_linetype = "solid",
           set_name_color = "black",
           set_name_size = 8,
           text_color = "white",
           text_size = 8)

p1
save_plot(p1,filename = "ssGSEA/Model/Lasso_uniCOX_Sig_Venn",style = "g")

##预后
library(survminer)
Sdata <- svdata
Sdata$OS_time <-Sdata$OS_time/365
res.cut <- surv_cutpoint(Sdata,time = "OS_time",
                         event = "Event",
                         variables = Iterm$Variables,
                         minprop = 0.3)

res_cat<- surv_categorize(res.cut)

mysurv<- Surv(res_cat$OS_time,res_cat$Event)

i=Iterm$Variables[3]
Gene <- data.frame(Gene=vector(),Pval=vector())
for(i in colnames(res_cat)[3:ncol(svdata)]){
  group <- res_cat[,i]
  surv_dat <- data.frame(group=group)
  fit <- survfit(mysurv~group)
  
  group <-factor(group,levels = c("low","high"))
  data.survdiff <- survdiff(mysurv~group)
  p.val <- 1-pchisq(data.survdiff$chisq,length(data.survdiff$n)-1)
  HR <- (data.survdiff$obs[2]/data.survdiff$exp[2])/(data.survdiff$obs[1]/data.survdiff$exp[1])
  up95 <- exp(log(HR))+qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1])
  low95 <- exp(log(HR))-qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1])
  xx <- c(i,p.val,HR,up95,low95)
  HR <- paste("Hazard Ratio = ",round(HR,2),sep = "")
  CI <- paste("95% CI: ",paste(round(low95,2),round(up95,2),sep = "-"),sep = "")
  svsort <- svdata[order(svdata[,i]),]
  cutoff <- paste("cutoff = ",round(svsort[fit$n[2],i],2),sep = "")
  
  p <- ggsurvplot(fit,data = surv_dat,
                  conf.int = F,
                  censor=F,palette = c("#f87669","#2fa1dd"),
                  legend.title = i,
                  font.x=15,font.y=15,font.title=15,
                  legend.labs = c(paste0("High ","(",fit$n[1],")"),
                                  paste0("Low ","(",fit$n[2],")")),
                  font.lenend=20,
                   xlim=c(0,12),
                  xlab = "OS_time(Years)",
                  ylab = "Survival probablity",
                  break.x.by =2,
                  break.y.by = 0.2,ggtheme = theme_classic()+mytheme,
                  pval = paste(pval=ifelse(p.val < 0.001,"p < 0.001",paste("p = ",round(p.val,3),sep = "")),
                               HR,CI,sep = "\n")
  )
  if(p.val < 0.05){

    CyDataPro::save_plot(p,filename = paste0("ssGSEA/Model/SigVar_KM/",i,"_diff"),style = "x")
  }else{
    CyDataPro::save_plot(p,filename = paste0("ssGSEA/Model/SigVar_KM//",i,"_Nodiff"),style = "x")
  }
  Gene <- rbind(Gene,xx)
  message(i)
}
colnames(Gene) <- c("Variables","Pvalue","HR","Up95","Low95")
fwrite(Gene,file = "ssGSEA/Model/SigVar_KM/KM_result.csv")

##ROC
library(pROC)
Gene$Variables
plotd <- svdata
plotd$Event <- factor(plotd$Event,levels = c(1,0),labels = c("Alive","Death"))
res <- roc(Event~Ferroptosis+Autophagy+Cuproptosis+NETosis+Paraptosis,data=plotd,smooth=T)
p<-ggroc(res,legacy.axes = T,size=1)+
      geom_segment(aes(x=0,xend=1,y=0,yend=1),color="black",linetype=4,linewidth=1)+
  theme_bw()+ggtitle("ssGSEA Score ROC")+ggsci::scale_color_lancet()+labs(color="Variables")+
  annotate("text",x=0.75,y=0.14,label=paste0("Ferroptosis AUC: ",round(res$Ferroptosis$auc,3)),size=6)+
  annotate("text",x=0.75,y=0.26,label=paste0("Autophagy AUC: ",round(res$Autophagy$auc,3)),size=6)+
  annotate("text",x=0.75,y=0.32,label=paste0("Cuproptosis AUC: ",round(res$Cuproptosis$auc,3)),size=6)+
  annotate("text",x=0.75,y=0.08,label=paste0("NETosis AUC: ",round(res$NETosis$auc,3)),size=6)+
  annotate("text",x=0.75,y=0.2,label=paste0("Paraptosis AUC: ",round(res$Paraptosis$auc,3)),size=6)+
  mytheme+theme(plot.title = element_text(hjust = 0.5))
save_plot(p,filename = "ssGSEA/Model/Sig_ROC_plot",style = "g")
library(timeROC)
col <- Disc5(5)
a3 <- plotd
a3$OS_year<- a3$OS_time/365
result <- with(a3,timeROC(T=OS_year,
                          delta = Event,
                          marker = Cuproptosis,
                          cause = 1,
                          times = c(0,0.5,1,1.5,2,2.5,3,3.5,4,4.5,5)),idd=T)
result1 <- with(a3,timeROC(T=OS_year,
                          delta = Event,
                          marker = Ferroptosis,
                          cause = 1,
                          times = c(0,0.5,1,1.5,2,2.5,3,3.5,4,4.5,5)),idd=T)
result2 <- with(a3,timeROC(T=OS_year,
                           delta = Event,
                           marker = Autophagy,
                           cause = 1,
                           times = c(0,0.5,1,1.5,2,2.5,3,3.5,4,4.5,5)),idd=T)
result3 <- with(a3,timeROC(T=OS_year,
                           delta = Event,
                           marker = NETosis,
                           cause = 1,
                           times = c(0,0.5,1,1.5,2,2.5,3,3.5,4,4.5,5)),idd=T)
result4 <- with(a3,timeROC(T=OS_year,
                           delta = Event,
                           marker = Paraptosis,
                           cause = 1,
                           times = c(0,0.5,1,1.5,2,2.5,3,3.5,4,4.5,5)),idd=T)
pdf("ssGSEA/Model/Sig_Time_AUCPlot.pdf",width = 8,height = 6)
plotAUCcurve(result,col = col[1])
plotAUCcurve(result1,col = col[2],add = T)
plotAUCcurve(result2,col = col[3],add = T)
plotAUCcurve(result3,col = col[4],add = T)
plotAUCcurve(result4,col = col[5],add = T)
legend("topright",c("Cuproptosis","Ferroptosis","Autophagy","NETosis","Paraptosis"),
       col=col,
       bty="n",
       lty=2,
       lwd=3,
       cex=1)
dev.off()
###GEO验证
library(GEOquery)
library(AnnoProbe)
library(limma)
GEODownload("GSE39582",outdir = "data/",GPL = "../自身免疫病_骨关节炎_MR20/data/GEO/GPL570-55999.txt",GPLskip = 16,Symbol_col = "Gene Symbol")
getOption("timeout")
options(timeout = 1000)
load("data/GSE39582.Rdata")
table(Sampleinfo$characteristics_ch1.1)
exp_log<- dectect_log(GEO_Smb_expr)
TuS <- Sampleinfo %>% filter(!(str_detect(Sampleinfo$characteristics_ch1.1,"Non")))
TuS <- TuS %>% filter(`os.delay:ch1` !="N/A")
TuS <- TuS %>% filter(`os.event:ch1` !="N/A")
inn <- intersect(rownames(TuS),colnames(GEO_Smb_expr))
TuS <- TuS[inn,]
Texp <- GEO_Smb_expr[,inn]
TuS$Event <- ifelse(TuS$`os.event:ch1` ==1,0,1)
TuS$OS_time <- as.numeric(TuS$`os.delay:ch1`)
TuS$Gender <- TuS$`Sex:ch1`
TuS$Stage <- TuS$`tnm.stage:ch1`
TuS$T_stage <- TuS$`tnm.t:ch1`
TuS$M_stage <- TuS$`tnm.m:ch1`
TuS$N_stage <- TuS$`tnm.n:ch1`
TuS$Age <- TuS$`age.at.diagnosis (year):ch1` %>% as.numeric()
TuS <- TuS %>% filter(T_stage !="N/A" & T_stage !="Tis")
TuS <- TuS %>% filter(M_stage !="N/A" & M_stage !="MX")
TuS <- TuS %>% filter(N_stage !="N/A" & N_stage !="N+")
TuS <- TuS %>% filter(Stage !=0)
meta <- TuS %>% select(Event,OS_time,Gender,Stage,T_stage,M_stage,N_stage,Age)

Item <- Gene$Variables
PCD <- PCDFliter[lengths(PCDFliter) >=5]
genelist <- as.list(PCD)
pcd <- PCD[Item]
genelist <- as.list(pcd)
gsea<- GSVA::gsva(as.matrix(Texp),genelist,method="ssgsea")
ES <- as.data.frame(t(gsea))
ES_zscore <- scale(ES) %>% data.frame() 
dir.create("ssGSEA/GSE39582")
pheatmap::pheatmap(t(ES_zscore),scale = "none",show_colnames = F,
                   fontsize = 15,main = "ssGSEA Score",
                   color = colorRampPalette(c("#6A3D9A","white","#E31A1C"))(1000),
                   filename = "ssGSEA/GSE39582/ssGSEA_result_heatmap.pdf",
                   width = 10,height = 6)
fwrite(ES,file = "ssGSEA/GSE39582/ssGSEA_result.csv",row.names = T)
fwrite(ES_zscore,file ="ssGSEA/GSE39582/ssGSEA_result_zscore.csv")
library(rstatix)
library(Hmisc)
library(corrplot)
cor <- rcorr(as.matrix(ES),type = "spearman")
R <- cor$r
P <- cor$P
P[is.na(P)] =0
COL
pdf("ssGSEA/GSE39582/ssGSEA_corr_heatmap.pdf",width = 8,height = 6)
corrplot(as.matrix(R),
         p.mat = as.matrix(P),
         insig = "label_sig",
         sig.level = .05,
         pch.cex = 1,
         method = "square",
         order = "alphabet",
         type = "lower",
         tl.col = "black",
         tl.pos = "ld",
         tl.srt=45,
         is.corr=T,
         col = rev(COL2("PRGn",10)),
         cl.pos = "b")

dev.off()
dev.new()
fwrite(as.data.frame(R),file = "ssGSEA/GSE39582/ssGSEA_corr_spearman_R.csv",row.names = T)
fwrite(as.data.frame(R),file = "ssGSEA/GSE39582/ssGSEA_corr_spearman_P.csv",row.names = T)

###ssGAES Score与肿瘤不同阶段
dir.create("ssGSEA/GSE39582/clinCorr")
library(ComplexHeatmap)
tmp<- meta
inn<- intersect(rownames(ES),rownames(tmp))
ES1 <- ES_zscore[inn,]
tmp1 <- tmp[inn,] 
identical(rownames(ES1),rownames(tmp1))
i=colnames(ES1)[1]
for(i in colnames(ES1)){
  meta1 <- tmp1
  meta1[,i] <- ES1[,i]
  meta1$OS_year <- as.numeric(meta1$OS_time)/12
  data <- meta1[,c(1,3:ncol(meta1))]
  data<-data[order(data[,i]),]
  data%>%glimpse()
  
  survdata<-data.frame(row.names = rownames(data),OS=data[,ncol(data)])
  data <- data[,-ncol(data)]
  data$Event <- factor(data$Event,levels = c("1","0"),labels = c("Alive","Death"))
  data$Stage <- factor(data$Stage,levels = c(1,2,3,4),labels = c("I","II","III","IV"))
  data$Age <- ifelse(data$Age  < 65,"<65",">=65")
  colnames(data)[8] <- "ssGSEA Score"
  table(data$Stage)
  ha=HeatmapAnnotation(survival_time = anno_points(survdata,size =unit(0.1,'cm'),#点注释
                                                   gp=gpar(col="grey")),df=data,#点注释属性及其他临床特征注释
                       col = list(Age=c("<65"="#66C2A5",">=65"="#FC8D62"),#设置部分临床特征注释颜色
                                  Event=c("Alive"="#A6D854","Death"="#E78AC3"),
                                  Stage = c("I"="#E64B35B2","II"="#4DBBD5B2","III"="#3C5488B2","IV"="#00A087B2"),
                                  Gender = c("Female" = "#BC3C29FF","Male" = "#0072B5FF"),
                                  `ssGSEA Score`=circlize::colorRamp2(c(-2,0,3), c("blue", "white", "red"))))
  heat<-Heatmap(matrix(nrow=0, ncol=nrow(data)),top_annotation = ha,
                column_title = i,
                column_title_gp = gpar(fontsize = 20,
                                       fontface = 'bold', 
                                       # fill = 'red',
                                       col = 'black'),
                rect_gp = gpar(col= "#FFFFFF00"))
  
  pdf(paste0("ssGSEA/GSE39582/clinCorr/",i,"_clin_heatmap.pdf"),width = 8,height = 6)
  draw(heat,annotation_legend_side="bottom")
  dev.off()
}
dev.new()
i <- colnames()
library(ggpubr)
clin <- clins[1]
tmp1$Age <- ifelse(tmp1$Age  < 65,"<65",">=65")
clins <- colnames(tmp1)[c(-1,-2)]
chisq_result <- data.frame()
for(i in colnames(ES1)){
  tmp1$Group <- ifelse(as.numeric(ES1[,i]) > median(as.numeric(ES1[,i])),"High","Low")
  for(clin in clins){
    tab <- xtabs(~tmp1$Group+tmp1[,clin])
    tt<-chisq.test(tab)
    tmp <- c(i,clin,tt$p.value)
    chisq_result <- rbind(chisq_result,tmp)
    P <- ifelse(tt$p.value < 0.001,"< 0.001",round(tt$p.value,3))
    p<-ggplot(tmp1, aes(x = Group, fill = tmp1[,clin])) +
      geom_bar(width = 0.9, position = "fill") + # 百分比柱状图
      scale_fill_brewer(palette = "PuBu")  +
      scale_y_continuous(labels = scales::percent) +
      labs(title = paste0(i," and ",clin),
           y = "Percentage",fill="")+
      theme(panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),
            panel.background = element_blank(),
            axis.text.x = element_text(face = "bold",colour = "black",size=15),
            axis.text.y = element_text(face = "bold",colour = "black",size=15),
            axis.ticks.x = element_blank(),
            axis.title.y = element_blank(),
            axis.title.x = element_blank(),
            plot.title = element_text(hjust = 0.5, face="bold"))+
      coord_cartesian(clip = 'off', ylim = c(0,1)) + 
      theme(plot.margin = margin(0.5,0.5,1.4,0.5,'cm'))+ #自定义图片上左下右的边框宽度
      annotate( "text",
                cex=6,
                x=1.5, y=-0.15, # 根据自己的数据调节p value的位置
                label=paste0("Pval= ",P), # 添加P值
                color="black"
      )+
      annotate("rect", xmin = 0.55, xmax = 1.45, ymin = -0.1, ymax = -0.05, 
               fill = "#40a1cf")+
      annotate("rect", xmin = 1.55, xmax = 2.45, ymin = -0.1, ymax = -0.05, 
               fill = "#dd816d")
    p
    dir.create(paste0("ssGSEA/GSE39582//clinCorr/batplot/",i),recursive = T)
    if(P < 0.05){
      ggsave(filename = paste0("ssGSEA/GSE39582//clinCorr/batplot/",i,"/",clin,"_Diff_barplot.pdf"),plot = p,width = 8,height = 6)
    }else{
      ggsave(filename = paste0("ssGSEA/GSE39582//clinCorr/batplot/",i,"/",clin,"_NoDiff_barplot.pdf"),plot = p,width = 8,height = 6)
    }
    
  }
  
}
-log10(0.05)
colnames(chisq_result)<- c("PCD","Clin","Pvalue")
chisq_result$logP <- -log10(as.numeric(chisq_result$Pvalue))
plot <- chisq_result %>% group_by(Clin) %>% arrange(Pvalue)
plot$PCD <- factor(plot$PCD,levels = unique(plot$PCD))
plot$Clin <- factor(plot$Clin,levels = unique(plot$Clin))
p<-ggplot(plot,aes(PCD,Clin))+
  geom_point(aes(size=logP,color=logP))+ #以-Log10(FDR)代表点的大小，Log2FC代表点的颜色
  scale_color_gradient2(low="#3F8CD1",mid = "#F5F5F5",high = "#F08080",midpoint = 1.31)+#设置颜色值
  scale_size_continuous(range = c(2,8))+ #调整点的大小为4-10
  theme_bw()+ #删除背景灰色
  theme(legend.position = "bottom")+ #修改图例的位置
  labs(color = "-log10(Pvalue)", size = "-log10(Pvalue)")+ # 修改图例的名字
  xlab(NULL)+ #删除横坐标
  ylab(NULL)+ #删除纵坐标
  mytheme+theme(axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1))
p
save_plot(p,filename = "ssGSEA/GSE39582/clinCorr/chiqs_test_pointPlot",style = "g",width = 10)

fwrite(chisq_result,file = "ssGSEA/GSE39582/clinCorr/chiqs_test_result.csv")

##单变量COX
library(survival)
dir.create("ssGSEA/GSE39582/Model/COX",recursive = T)
pFilter <- 0.05 
uniCox <- data.frame()
uniCoxSig <- data.frame()
identical(rownames(meta),rownames(ES1))
svdata<- cbind(meta,ES1)
table(svdata$Event)
svdata$Event <- ifelse(svdata$Event==1,0,1)
Iterm
save(svdata,file = "data/GSE39582_Tumor_Surv.Rdata")
for(j in colnames(svdata[,Iterm$Variables])){
  cox <- coxph(Surv(OS_time,Event) ~ svdata[,j], data = svdata)
  coxSummary <- summary(cox)
  uniCox <- rbind(uniCox, data.frame(gene = j,
                                     HR = coxSummary$conf.int[,"exp(coef)"],
                                     HR.95L = coxSummary$conf.int[,"lower .95"],
                                     HR.95H = coxSummary$conf.int[,"upper .95"],
                                     P_value = coxSummary$coefficients[,"Pr(>|z|)"]))
  coxP <- coxSummary$coefficients[,"Pr(>|z|)"]
  if(coxP<pFilter){
    uniCoxSig <- rbind(uniCoxSig, data.frame(variable = j, 
                                             HR = coxSummary$conf.int[,"exp(coef)"],
                                             HR.95L = coxSummary$conf.int[,"lower .95"],
                                             HR.95H = coxSummary$conf.int[,"upper .95"],
                                             P_value = coxSummary$coefficients[,"Pr(>|z|)"]))
  }
}

# add prognostic label
uniCox$Risk_group <- ifelse(uniCox$HR >= 1 & uniCox$P_value < 0.05, "Risk",
                            ifelse(uniCox$HR < 1 & uniCox$P_value < 0.05, "Protect","Not_sig"))
uniCoxSig$Risk_group <- ifelse(uniCoxSig$HR >= 1, "Risk", "Protect")
uniCox <- arrange(uniCox,HR)
uniCox$gene<- factor(uniCox$gene,levels = uniCox$gene)
p<-ggplot(uniCox, aes(HR, gene)) + 
  geom_errorbarh(aes(xmax = HR.95L, xmin = HR.95H,col= Risk_group),height=0.3,size=0.7) +
  geom_point(shape=18, aes(col=Risk_group, size = -log10(P_value)))  +
  geom_vline(aes(xintercept = 1),color="black",linetype = "dashed",size=1) + 
  xlim(0,1.5)+labs(title = "uniCox",y="")+
  scale_color_manual(values = c("gray","blue","red"))+
  theme_bw()+mytheme+theme(plot.title = element_text(hjust = 0.5))
p
CyDataPro::save_plot(p,filename = "ssGSEA/GSE39582/Model/COX/uniCox_foresplot",style = "g")
fwrite(uniCox,file = "ssGSEA/GSE39582/Model/COX/uniCOX_result.csv")


##预后
library(survminer)
Sdata <- svdata
Sdata$Event <- ifelse(Sdata$Event ==1,0,1)
Sdata$OS_time <-Sdata$OS_time/12
res.cut <- surv_cutpoint(Sdata,time = "OS_time",
                         event = "Event",
                         variables = Iterm$Variables,
                         minprop = 0.3)

res_cat<- surv_categorize(res.cut)

mysurv<- Surv(res_cat$OS_time,res_cat$Event)

i=Iterm$Variables[3]
Gene <- data.frame(Gene=vector(),Pval=vector())
for(i in colnames(res_cat)[3:ncol(svdata)]){
  group <- res_cat[,i]
  surv_dat <- data.frame(group=group)
  fit <- survfit(mysurv~group)
  
  group <-factor(group,levels = c("high","low"))
  data.survdiff <- survdiff(mysurv~group)
  p.val <- 1-pchisq(data.survdiff$chisq,length(data.survdiff$n)-1)
  HR <- (data.survdiff$obs[2]/data.survdiff$exp[2])/(data.survdiff$obs[1]/data.survdiff$exp[1])
  up95 <- exp(log(HR))+qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1])
  low95 <- exp(log(HR))-qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1])
  xx <- c(i,p.val,HR,up95,low95)
  HR <- paste("Hazard Ratio = ",round(HR,2),sep = "")
  CI <- paste("95% CI: ",paste(round(low95,2),round(up95,2),sep = "-"),sep = "")
  svsort <- svdata[order(svdata[,i]),]
  cutoff <- paste("cutoff = ",round(svsort[fit$n[2],i],2),sep = "")
  
  p <- ggsurvplot(fit,data = surv_dat,
                  conf.int = F,
                  censor=F,palette = c("#f87669","#2fa1dd"),
                  legend.title = i,
                  font.x=15,font.y=15,font.title=15,
                  legend.labs = c(paste0("High ","(",fit$n[1],")"),
                                  paste0("Low ","(",fit$n[2],")")),
                  font.lenend=20,
                   xlim=c(0,12),
                  xlab = "OS_time(Years)",
                  ylab = "Survival probablity",
                  break.x.by =1,
                  break.y.by = 0.2,ggtheme = theme_classic()+mytheme,
                  pval = paste(pval=ifelse(p.val < 0.001,"p < 0.001",paste("p = ",round(p.val,3),sep = "")),
                               HR,CI,sep = "\n")
  )
  p
  if(p.val < 0.05){
    
    CyDataPro::save_plot(p,filename = paste0("ssGSEA/GSE39582/Model/SigVar_KM/",i,"_diff"),style = "x")
  }else{
    CyDataPro::save_plot(p,filename = paste0("ssGSEA/GSE39582/Model/SigVar_KM//",i,"_Nodiff"),style = "x")
  }
  Gene <- rbind(Gene,xx)
  message(i)
}
dir.create("ssGSEA/GSE39582/Model/SigVar_KM")
colnames(Gene) <- c("Variables","Pvalue","HR","Up95","Low95")
fwrite(Gene,file = "ssGSEA/GSE39582/Model/SigVar_KM/KM_result.csv")

devtools::install_github("gpli/DNB")
library(DNB)
rm(list = ls())
load("../结直肠癌_RNAmSNP_MR_Singlecell16/data/TCGA_CRC.Rdata")
load("data/TCGA_Tumor_surv.Rdata")
load("data/PCD_COX_filter.Rdata")

Iterm <- fread("ssGSEA/Model/Lasso_uniCOX_Sig_result.csv")
PCD <- PCDFliter[Iterm$Variables]

exp <- dectect_log(count_matix)
inn <- intersect(colnames(exp),rownames(TCGA_Tumor_Surv))
exp <- exp[,inn]
Sample <- TCGA_Tumor_Surv[inn,]
Sample<- Sample %>% mutate(time = recode(Stage,"I"=1,"II"=2,"III"=3,"IV"=4))
identical(rownames(Sample),colnames(exp))
i=names(PCD)[2]
class(Iterm)
for(i in names(PCD)){
  dir.create(paste0("DNB/",i))
  gene <- PCD[i][[1]]
  exp1 <- exp[intersect(gene,rownames(exp)),]
  dnb <- new_DNB(data = exp1,
                 time = as.factor(Sample$time),)
  dnb <- cal_cor(dnb)
  dnb <- cal_cv(dnb)
  dnb <- search_candidates(dnb,with_ctrl = F,min_size = 3)
  dnb <- cal_final(dnb)
  dnb_genes <- get_DNB_genes(dnb)
  dnb_genes <- data.frame(Gene=dnb_genes)
  fwrite(dnb_genes,file = paste0("DNB/",i,"/DNB_gene.csv"))
  candidates <- get_candidates(dnb)
  plot_DNB(candidates)
  final <- get_final(dnb)
  final$CI <- abs(final$score)
  final <- final %>% dplyr::mutate(Stage = recode(time,"1" ="I","2" ="II","3"="III","4"="IV"))
  a<-top_n(final,1,CI)
 
  CI<-ggplot(final,aes(Stage,CI,group=1))+
    geom_point(color="red",size=4)+
    geom_line(color = "red",linewidth=1.5)+
    labs(title = "CI",y="")+
    geom_segment(x=a[1,8],xend=a[1,8],y=0,yend=a[1,7],linetype="dashed",linewidth=1.5,color="blue")+
    theme_bw()+mytheme+theme(plot.title = element_text(hjust = 0.5))
  save_plot(CI,filename = paste0("DNB/",i,"/CI_plot"),style = "g")
  
  CV<-ggplot(final,aes(Stage,cv_in,group=1))+
    geom_point(color="red",size=4)+
    geom_line(color = "red",linewidth=1.5)+
    labs(title = "CV",y="")+
    theme_bw()+mytheme+theme(plot.title = element_text(hjust = 0.5))  
  save_plot(CV,filename = paste0("DNB/",i,"/CV_plot"),style = "g")
  PCCin<-ggplot(final,aes(Stage,abs(cor_in),group=1))+
    geom_point(color="red",size=4)+
    geom_line(color = "red",linewidth=1.5)+
    labs(title = "PCCin",y="")+
    theme_bw()+mytheme+theme(plot.title = element_text(hjust = 0.5))  
  save_plot(PCCin,filename = paste0("DNB/",i,"/PCCin_plot"),style = "g")
  Pccout<-ggplot(final,aes(Stage,abs(cor_out),group=1))+
    geom_point(color="red",size=4)+
    geom_line(color = "red",linewidth=1.5)+
    labs(title = "PCCout",y="")+
    theme_bw()+mytheme+theme(plot.title = element_text(hjust = 0.5)) 
  save_plot(Pccout,filename = paste0("DNB/",i,"/Pccout_plot"),style = "g")
  fwrite(final,file = paste0("DNB/",i,"/final_result.csv"))
  dir.create(paste0("DNB/",i,"/corr_result"))
  StageI <- dnb[["correlation"]][["1"]] %>% as.data.frame()
  StageII <- dnb[["correlation"]][["2"]] %>% as.data.frame()
  StageIII <- dnb[["correlation"]][["3"]] %>% as.data.frame()
  StageIV <- dnb[["correlation"]][["4"]] %>% as.data.frame()
  fwrite(StageI,file = paste0("DNB/",i,"/corr_result/StageI.csv"),row.names = T)
  fwrite(StageII,file = paste0("DNB/",i,"/corr_result/StageI.csv"),row.names = T)
  fwrite(StageIII,file = paste0("DNB/",i,"/corr_result/StageI.csv"),row.names = T)
  fwrite(StageIV,file = paste0("DNB/",i,"/corr_result/StageI.csv"),row.names = T)
  print(paste0(i," 完成"))
  
}

###多变量cOX
library(ggprism)
remotes::install_github("yikeshu0611/ggrisk")
install.packages("ggrisk")
library(ggrisk)
library(ggDCA)
library(Cycolors)
library(ezcox)
ggrisk
rm(list = ls())
load("data/TCGA_surv.Rdata")
load("../结直肠癌_RNAmSNP_MR_Singlecell16/data/TCGA_CRC.Rdata")
exp <- dectect_log(count_matix)
meta <- svdata[,c(12,13,6,7,14:17)]
inn<- intersect(rownames(meta),colnames(exp))
meta<- meta[inn,]
exp1<- exp[,inn] %>% t() %>% as.data.frame()
identical(rownames(meta),rownames(exp1))
Items <- fread("ssGSEA/Model/Lasso_uniCOX_Sig_result.csv")
Items <- Items$Variables
i=Items[2]
class(meta$OS_time)
edit(ggrisk)
mycol2 <- Disc5(3)

mycol <- Disc9(3)
for(i in Items){
  gene <- fread(paste0("DNB/",i,"/DNB_gene.csv"))
  dir.create(paste0("DNB/",i,"/multiCox/"))
  inn <- intersect(colnames(exp1),gene$Gene)
  exp_tmp <- exp1[,inn]
  meta1 <- cbind(meta,exp_tmp)
  meta2 <- meta1 %>% select(all_of(c("OS_time","Event",inn)))
  cox <- coxph(Surv(OS_time,Event)~., data = meta2)
  p<-two_scatter(cox,
         cutoff.value = "median",
         cutoff.y = 1,size.cutoff = 5,
         size.dashline = 1.2)
  p
  save_plot(p,filename = paste0("DNB/",i,"/multiCox/Risk_Plot"),style = "g")
  RS <- predict(cox,type = "risk",newdata = meta2) %>% as.data.frame()
  RS<- RS %>% rownames_to_column()
  colnames(RS) <- c("sample","risk_score")
  fwrite(RS,paste0("DNB/",i,"/multiCox/RiskScore.csv"))
  clin1 <- meta1 %>% rownames_to_column("sample")
  RS_meta <-RS %>% inner_join(clin1)
  RS_meta$Group <- ifelse(RS_meta$risk_score > median(RS_meta$risk_score),"High","Low")
  fit <- survfit(Surv(OS_time,Event)~Group,data = RS_meta)
  data.survdiff <- survdiff(Surv(OS_time,Event)~Group,data = RS_meta)
  p.val <- 1-pchisq(data.survdiff$chisq,length(data.survdiff$n)-1)
  HR <- (data.survdiff$obs[2]/data.survdiff$exp[2])/(data.survdiff$obs[1]/data.survdiff$exp[1])
  up95 <- exp(log(HR))+qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1])
  low95 <- exp(log(HR))-qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1])
  
  HR <- paste("Hazard Ratio = ",round(HR,2),sep = "")
  CI <- paste("95% CI: ",paste(round(low95,2),round(up95,2),sep = "-"),sep = "")
  p<-ggsurvplot(fit,data = RS_meta,
                censor.shape="|",censor.size=6,
                conf.int = F,
                con.int.style ="ribnbon",
                con.in.alpha = 0.2,
                pval = paste(pval=ifelse(p.val < 0.001,"p < 0.001",paste("p = ",round(p.val,3),sep = "")),
                             HR,CI,sep = "\n"),
                palette = "npg",
                surv.median.line = "hv",
                ggtheme = theme_bw()+mytheme+theme(plot.title = element_text(hjust = 0.5)),
                legend ="top",
                legend.abs = c("High","Low"),
                xlab = "OS_time(Days)",
                ylab = "Survival probablity",
                title = paste0(i," Gene Risk Model"),
                break.x.by =1000,
                break.y.by = 0.2,
                risk.table = F)
  p
  save_plot(p,filename = paste0("DNB/",i,"/multiCox/KM_plot"),style = "xx")
  fdca <- dca(cox,model.names=c("RiskScores"))
  p<-ggplot(fdca,
            linetype = F,lwd=1.5)+
    labs(title = i)+
    theme_classic()+  
    theme(legend.position="bottom",plot.title = element_text(hjust = 0.5))+
    # scale_x_continuous(
    #   limits = c(0.25, 1)) +
    scale_y_continuous(
      limits = c(-0.05, 0.4))+
    scale_color_manual(values = mycol)+
    mytheme
  p
  save_plot(p,filename = paste0("DNB/",i,"/multiCox/RiskScore_DCAPlot"),style = "g")
  coxph(Surv(OS_time,Event)~risk_score, data = RS_meta)
  
  variable.names <- c("Gender","AgeGroup","Stage","risk_score")
  result  <- ezcox(RS_meta,time = "OS_time",status = "Event",covariates = variable.names,return_models = T)
  dir.create("TCGA/IndependentFactor")
  mod <- get_models(result)
  a <- result$res %>% as.data.frame()
  write.table(a,file = paste0("DNB/",i,"/multiCox/Clin_Unicox_result.txt"),row.names = F,sep = "\t",quote = F)
  p<-show_models(mod,model_names = paste0("Model ", 1:4))
  save_plot(p,filename = paste0("DNB/",i,"/multiCox/Clin_Unicox_result"),style = "xx")
  
  td <- cbind(RS_meta[,3:4],RS_meta[,variable.names])
  tdmultiCox=coxph(Surv(OS_time,Event) ~ ., data = td)
  x <-summary(tdmultiCox) 
  pvalue=signif(as.matrix(x$coefficients)[,5],2)
  HR=signif(as.matrix(x$coefficients)[,2],2)
  low=signif(x$conf.int[,3],2)
  high=signif(x$conf.int[,4],2)
  multi_res=data.frame(p.value=pvalue,
                       HR=paste(HR," (",low,"-",high,")",sep=""),
                       stringsAsFactors = F
  )
  write.table(multi_res,file =paste0("DNB/",i,"/multiCox/Clin_multicox_result.txt"),row.names = F,quote = F,sep = "\t")
  p <-ggforest(tdmultiCox,
               main = "Hazard ratio",
               cpositions = c(0.02,0.15, 0.35), #前三列的位置，第二列是样品数，设了个负值，相当于隐藏了
               fontsize = 1.3, #字体大小
               refLabel = "reference", 
               noDigits = 2)
  p
  
  save_plot(p,filename = paste0("DNB/",i,"/multiCox/Clin_multicox_foresplot"),style = "g",width = 12)
  dca_cox <- dca(tdmultiCox,model.names = c("AgeGroup+Gender+riskScore+Stage"))

  p <-ggplot(dca_cox,
             linetype = F,lwd=1.2)+
    labs(title = i)+
    theme_classic()+  
    theme(legend.position="bottom",plot.title = element_text(hjust = 0.5))+
    # scale_x_continuous(
    #   limits = c(0.25, 1)) +
    # scale_y_continuous(
    #   limits = c(-0.05, 0.4))+
    scale_color_manual(values = mycol2)+
    mytheme
  save_plot(p,filename = paste0("DNB/",i,"/multiCox/Clin_multicox_DCA"),style = "g",width = 10,height = 6)
  
}
data <- fread("DNB/Ferroptosis/multiCox/RiskScore.csv")
colnames(data)[2]<- "Ferroptosis"
for(i in Items[-1]){
  rsk <- fread(paste0("DNB/",i,"/multiCox/RiskScore.csv"))
  colnames(rsk)[2]<- i
  data <- inner_join(data,rsk)
}


data <- data %>% column_to_rownames("sample")
identical(rownames(data),rownames(meta))
data$Event <- meta$Event
data$OS_time <- meta$OS_time
plotd <- data
plotd$Event <- factor(plotd$Event,levels = c(1,0),labels = c("Alive","Death"))
coll <- Disc9(5)
res <- roc(Event~Ferroptosis+Autophagy+Cuproptosis+NETosis+Paraptosis,data=plotd,smooth=T)
p<-ggroc(res,legacy.axes = T,size=1)+
  geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1), color="black", linetype=8,size=2)+
  theme_bw()+ggtitle('Risk Score ROC')+
  scale_color_manual(name = NULL,values = coll)+
  annotate("text",x=0.75,y=0.2,label=paste0("Ferroptosis AUC: ",round(res$Ferroptosis$auc,3)),size=5)+
  annotate("text",x=0.75,y=0.26,label=paste0("Autophagy AUC: ",round(res$Autophagy$auc,3)),size=5)+
  annotate("text",x=0.75,y=0.14,label=paste0("Cuproptosis AUC: ",round(res$Cuproptosis$auc,3)),size=5)+
  annotate("text",x=0.75,y=0.08,label=paste0("NETosis AUC: ",round(res$NETosis$auc,3)),size=5)+
  annotate("text",x=0.75,y=0.32,label=paste0("Paraptosis AUC: ",round(res$Paraptosis$auc,3)),size=5)+
  mytheme+theme(plot.title = element_text(hjust = 0.5))
p

save_plot(p,filename = "DNB/RiskScore_ROC",style = "g")
library(timeROC)
col <- Disc9(5)
a3 <- plotd
a3$OS_year<- a3$OS_time/365
result <- with(a3,timeROC(T=OS_year,
                          delta = Event,
                          marker = Cuproptosis,
                          cause = 1,
                          times = c(0,0.5,1,1.5,2,2.5,3,3.5,4,4.5,5)),idd=T)
result1 <- with(a3,timeROC(T=OS_year,
                           delta = Event,
                           marker = Ferroptosis,
                           cause = 1,
                           times = c(0,0.5,1,1.5,2,2.5,3,3.5,4,4.5,5)),idd=T)
result2 <- with(a3,timeROC(T=OS_year,
                           delta = Event,
                           marker = Autophagy,
                           cause = 1,
                           times = c(0,0.5,1,1.5,2,2.5,3,3.5,4,4.5,5)),idd=T)
result3 <- with(a3,timeROC(T=OS_year,
                           delta = Event,
                           marker = NETosis,
                           cause = 1,
                           times = c(0,0.5,1,1.5,2,2.5,3,3.5,4,4.5,5)),idd=T)
result4 <- with(a3,timeROC(T=OS_year,
                           delta = Event,
                           marker = Paraptosis,
                           cause = 1,
                           times = c(0,0.5,1,1.5,2,2.5,3,3.5,4,4.5,5)),idd=T)
pdf("DNB/RiskScore_Time_AUCPlot.pdf",width = 8,height = 6)
plotAUCcurve(result,col = col[1])
plotAUCcurve(result1,col = col[2],add = T)
plotAUCcurve(result2,col = col[3],add = T)
plotAUCcurve(result3,col = col[4],add = T)
plotAUCcurve(result4,col = col[5],add = T)
legend("topright",c("Cuproptosis","Ferroptosis","Autophagy","NETosis","Paraptosis"),
       col=col,
       bty="n",
       lty=2,
       lwd=3,
       cex=1)
dev.off()
library(VennDetail)
xx <- fread(("DNB/DNB_gene.csv"))
x <- as.list(xx)
ven <- VennDetail::venndetail(x)
plot(ven,type = "upset")
a<- result(ven)
install.packages("tinyarray")
library(TCGAplot)
tpm<-get_all_tpm()
Disc7(7)
"#32A852" "#63D07F" "#6DB1EE" "#1877C9" "#C9244C"
"#6BB592" "#36ACA2" "#C4A751" "#EC748B" "#B67FB3"
"#D55E86"  "#C9244C" "#32A852" "#87447A" "#63D07F" "#1877C9" "#6DB1EE"
"#f87669","#2fa1dd"
gene_ano <- data.frame(Gene = c(xx$Autophagy,xx$Cuproptosis[1:3],xx$Ferroptosis[1:4],xx$NETosis[1:4],xx$Paraptosis[1:4]),
                       Class = c(rep("Autophagy",39),rep("Cuproptosis",3),rep("Ferroptosis",4),rep("NETosis",4),rep("Paraptosis",4)))
gene_ano<- gene_ano %>% column_to_rownames("Gene") %>% as.data.frame()
exp <- tpm %>% filter(Group =="Tumor")
aa<-table(tpm$Cancer,tpm$Group) %>% as.data.frame() %>% filter(Var2 =="Normal")
aa <- aa %>% filter(Freq> 0)
exp <- tpm %>% filter(Cancer %in% aa$Var1[1])
Sampano <- data.frame(Sample = rownames(exp),
                      Group = exp$Group)
Sampano<- Sampano %>% column_to_rownames("Sample") %>% as.matrix()
exp1 <- exp %>% select(all_of(c(rownames(gene_ano))))
ROW_A <- HeatmapAnnotation(
  Class=gene_ano,which = "row",show_annotation_name = F,
  col = list(Class=c("Autophagy"="#32A852","Cuproptosis"="#6DB1EE","Ferroptosis"="#C9244C","NETosis"="#C4A751","Paraptosis"="#87447A"))
)
COL_A<- HeatmapAnnotation(Group = Sampano,which = "column",show_annotation_name = F,
                          col = list(Group = c("Tumor" = "#f87669","Normal" = "#2fa1dd")))
dev.off()
dev.new()
Heatmap(t(exp1),cluster_rows = T,cluster_columns = T,use_raster = T,
        row_names_gp = gpar(fontsize=8),column_names_gp = gpar(fontsize=8),
        row_split = gene_ano,column_split = Sampano,
        row_title_gp = gpar(fontsize=10),column_title_gp = gpar(fontsize=10),
        right_annotation = ROW_A,show_column_names = F,
        top_annotation = COL_A)

table(exp$Cancer)
ggplot(exp,aes(x=Cancer,y=PIK3C3,fill=Group))+
  geom_boxplot()
pheatmap::pheatmap(exp,show_rownames = F)
pan_tumor_boxplot(rownames(gene_ano)[1])
p<-tcga_kmplot("BLCA","PIK3C3")
edit(tcga_kmplot)
a<-p$data.survplot
edit(gene_deg_heatmap)
gene = rownames(gene_ano)
##泛癌DEG

cancer = aa$Var1[1]
for(cancer in aa$Var1){
  exp = subset(tpm,  Cancer == cancer)
  Group = exp$Group
  Group = factor(Group, levels = c("Normal", "Tumor"))
  exp = exp[, -c(1:2)]
  exp = as.matrix(t(exp))
  dat = normalizeBetweenArrays(exp)
  design = model.matrix(~Group)
  fit = lmFit(dat, design)
  fit = eBayes(fit)
  options(digits = 4)
  DEG = topTable(fit, coef = 2, adjust = "BH", n = Inf)
  DEG = na.omit(DEG)
  DEG = subset(DEG, P.Value < 0.05)
  DEG = DEG %>% rownames_to_column("Gene")
  K1 <- (DEG$P.Value<0.05)&(DEG$logFC < -0.7)
  K2 <- (DEG$P.Value<0.05)&(DEG$logFC > 0.7)
  DEG$change <- ifelse(K1,"Down",ifelse(K2,"Up","not"))
  G_DEG =DEG %>% filter(Gene %in% gene)
  G_DEG = G_DEG %>% select(Gene,logFC,P.Value,change)
  
  colnames(G_DEG)[c(2,3,4)] = c(paste0(cancer,"_logFC"),paste0(cancer,"_Pvalue"),cancer)
  fwrite(G_DEG,file = paste0("TCGA_PAN/",cancer,"_DEG.csv"))
  print(paste0(cancer,"wanc"))
}
cancer = aa$Var1[1]
data <- fread(paste0("TCGA_PAN/DEG/",cancer,"_DEG.csv"))
data <- data %>% select(Gene,all_of(cancer))
for(cancer in aa$Var1[-1]){
  tmp <- fread(paste0("TCGA_PAN/DEG/",cancer,"_DEG.csv"))
  tmp <- tmp%>% select(Gene,all_of(cancer))
  data <- full_join(data,tmp)
}

Gene <- dft %>% filter(value != "not expression" & value !="not change")
table(Gene$Cancer,Gene$Gene)
m  <- aggregate(Gene$Gene, by = list(Gene$name), FUN = length)
colnames(m)<-c("Cancer","GeneNumber")
m<-m %>% arrange(GeneNumber)
m$Cancer<- factor(m$Cancer,levels = m$Cancer)
mycol = heat.colors(3)
brewer.pal(8,"Greens")
"#F7FCF5" "#E5F5E0" "#C7E9C0" "#A1D99B" "#74C476" "#41AB5D" "#238B45" "#005A32"
p<-ggplot(m,aes(y=Cancer,x=GeneNumber,fill=GeneNumber))+
  geom_bar(stat = "identity")+
  geom_text(aes(label = GeneNumber), vjust = 0.5,hjust=0.1,size=5,colour = "black")+
  labs(title = "DEG Number")+
  scale_fill_gradient(high = "#238B45",low = "#E5F5E0")+
  labs(x="Gene Number",fill="",y="Cancer")+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5))+
  mytheme
p
save_plot(p,filename = "TCGA_PAN/Pan_Cancer_DEG_Number_barplot",style = "g")
fwrite(m,file = "TCGA_PAN/DEG/DEG_number.csv")

data[is.na(data)]="not expression"
data[data=="not"] = "not change"
inn <- intersect(gene_ano$Gene,data$Gene)
dft <- data %>% filter(Gene %in% inn)
GeneA <- gene_ano %>% filter(Gene %in% inn)
identical(GeneA$Gene,dft$Gene)
dft <- data %>% pivot_longer(-Gene)


dft$Gene <- factor(dft$Gene,levels = GeneA$Gene)
p<-ggplot(data = dft,aes(x = name,y = Gene))+
  geom_tile(aes(fill = value),color = "grey")+
  theme_minimal()+
  theme(panel.grid = element_blank(),
        legend.position = "right",
        legend.text = element_text(size = 18),
        legend.title = element_text(size = 18),
        axis.ticks.y = element_blank(),
        axis.title = element_blank(),
        axis.text.x = element_text(angle = 45,hjust = 1,size = 15,color="black",vjust = 1),
        axis.text.y = element_text( size = 15,color="black"))+
  scale_fill_manual(values = c('#2fa1dd','gray',"white","#f87669"))+
  labs(fill = "Diff change")
p

Group <- data.frame(Gene = GeneA$Gene,
                    Xaxis = "Annotation",
                    Class = GeneA$Class)

Group$Gene <- factor(Group$Gene,levels = GeneA$Gene)
Block_Stage <- ggplot(data = Group,aes(Xaxis,Gene,fill = Class))+
  geom_tile()+
   scale_fill_manual(values = c("#32A852","#6DB1EE","#C9244C","#C4A751","#87447A"))+
  theme_void()+mytheme2
library(aplot)
p_block1 <- p %>% insert_right(Block_Stage,width = 0.07)
save_plot(p_block1,filename = "TCGA_PAN/PanCancer_DEG_heatmap",width = 10,height = 12,style = "g")
mytheme2 <- theme(panel.grid = element_blank(),
                  legend.position = "right",
                  legend.text = element_text(size = 18),
                  legend.title = element_text(size = 18),
                  axis.ticks.y = element_blank(),
                  axis.title = element_blank(),
                  axis.text.x =  element_blank(),
                  axis.text.y =  element_blank())

###泛癌预后
meta <- get_all_meta()
cancer = aa$Var1[1]
dir.create("TCGA_PAN/Surv/KMPlot",recursive = T)
Gene <- data.frame(Gene=vector(),Pval=vector())
for(cancer in aa$Var1){
  dir.create(paste0("TCGA_PAN/Surv/KMPlot/",cancer))
  exprSet = subset(tpm, Group == "Tumor" & Cancer == cancer) %>% 
    dplyr::select(all_of(gene)) %>% tibble::add_column(ID = stringr::str_sub(rownames(.), 
                                                                             1, 12)) %>% dplyr::filter(!duplicated(ID)) %>% tibble::remove_rownames(.) %>% 
    tibble::column_to_rownames("ID") %>% dplyr::filter(rownames(.) %in% 
                                                         rownames(subset(meta, Cancer == cancer))) %>% t() %>% 
    as.matrix()
  inn <-intersect(gene_ano$Gene,rownames(exprSet))
  exp <- exprSet[inn,]
  cl = meta[colnames(exp), ]
  exp <- t(exp)
  identical(rownames(exp),rownames(cl))
  Sdata <- cbind(cl[,c(2,3)],exp)
  res.cut <- surv_cutpoint(Sdata,time = "time",
                           event = "event",
                           variables = inn,
                           minprop = 0.3)
  res_cat<- surv_categorize(res.cut)
  mysurv<- Surv(res_cat$time,res_cat$event)
  for(i in colnames(res_cat)[3:ncol(Sdata)]){
    group <- res_cat[,i]
    surv_dat <- data.frame(group=group)
    fit <- survfit(mysurv~group)
    group <-factor(group,levels = c("low","high"))
    data.survdiff <- survdiff(mysurv~group)
    p.val <- 1-pchisq(data.survdiff$chisq,length(data.survdiff$n)-1)
    HR <- (data.survdiff$obs[2]/data.survdiff$exp[2])/(data.survdiff$obs[1]/data.survdiff$exp[1])
    up95 <- exp(log(HR))+qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1])
    low95 <- exp(log(HR))-qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1])
    xx <- c(cancer,i,p.val,HR,up95,low95)
    HR <- paste("Hazard Ratio = ",round(HR,2),sep = "")
    CI <- paste("95% CI: ",paste(round(low95,2),round(up95,2),sep = "-"),sep = "")
    svsort <- Sdata[order(Sdata[,i]),]
    cutoff <- paste("cutoff = ",round(svsort[fit$n[2],i],2),sep = "")
    
    p <- ggsurvplot(fit,data = surv_dat,
                    conf.int = F,
                    censor=T,palette = c("#f87669","black"),
                    legend.title = i,
                    font.x=15,font.y=15,font.title=15,
                    legend.labs = c(paste0("High ","(",fit$n[1],")"),
                                    paste0("Low ","(",fit$n[2],")")),
                    font.lenend=20,
                    # xlim=c(0,150),
                    xlab = "OS_time(Month)",
                    ylab = "Survival probablity",
                    # break.x.by =30,
                    break.y.by = 0.2,ggtheme = theme_classic()+mytheme,
                    pval = paste(pval=ifelse(p.val < 0.001,"p < 0.001",paste("p = ",round(p.val,3),sep = "")),
                                 HR,CI,sep = "\n")
    )
    if(p.val < 0.05){
      
      CyDataPro::save_plot(p,filename = paste0("TCGA_PAN/Surv/KMPlot/",cancer,"/",i,"_diff"),style = "x")
    }else{
      CyDataPro::save_plot(p,filename = paste0("TCGA_PAN/Surv/KMPlot/",cancer,"/",i,"_Ndiff"),style = "x")
    }
    Gene <- rbind(Gene,xx)
    message(i)
  }
  message(cancer)
}
colnames(Gene)<- c("Cancer","Gene","Pvalue","HR","Up95","Low95")
K1 <- (Gene$Pvalue<0.05)&(Gene$HR < 1)
K2 <- (Gene$Pvalue<0.05)&(Gene$HR > 1)
Gene$Change <- ifelse(K1,"Protective",ifelse(K2,"Risk","NotSig"))
table(Gene$Change)
Gene$Gene <- factor(Gene$Gene,levels = gene_ano$Gene)
p<-ggplot(data = Gene,aes(x = Cancer,y = Gene))+
  geom_tile(aes(fill = Change),color = "grey")+
  theme_minimal()+
  theme(panel.grid = element_blank(),
        legend.position = "right",
        legend.text = element_text(size = 18),
        legend.title = element_text(size = 18),
        axis.ticks.y = element_blank(),
        axis.title = element_blank(),
        axis.text.x = element_text(angle = 45,hjust = 1,size = 15,color="black",vjust = 1),
        axis.text.y = element_text( size = 15,color="black"))+
  scale_fill_manual(values = c("white",'#32A852',"#f87669"))+
  labs(fill = "Diff change")
p

Group <- data.frame(Gene = gene_ano$Gene,
                    Xaxis = "Annotation",
                    Class = gene_ano$Class)

Group$Gene <- factor(Group$Gene,levels = gene_ano$Gene)
Block_Stage <- ggplot(data = Group,aes(Xaxis,Gene,fill = Class))+
  geom_tile()+
  scale_fill_manual(values = c("#32A852","#6DB1EE","#C9244C","#C4A751","#87447A"))+
  theme_void()+mytheme2
library(aplot)
p_block1 <- p %>% insert_right(Block_Stage,width = 0.06)
save_plot(p_block1,filename = "TCGA_PAN/Surv/Surv_KM_heatmap",width = 10,height = 12,style = "g")
fwrite(Gene,file = "TCGA_PAN/Surv/Surv_KM_resylt.csv")
Gene <- fread("TCGA_PAN/Surv/Surv_KM_resylt.csv")
Gene <- Gene %>% filter(Change != "NotSig")
table(Gene$Cancer,Gene$Gene)
m  <- aggregate(Gene$Gene, by = list(Gene$Cancer), FUN = length)
colnames(m)<-c("Cancer","GeneNumber")
m<-m %>% arrange(GeneNumber)
m$Cancer<- factor(m$Cancer,levels = m$Cancer)
mycol = heat.colors(24)
p<-ggplot(m,aes(y=Cancer,x=GeneNumber,fill=GeneNumber))+
  geom_bar(stat = "identity")+
  geom_text(aes(label = GeneNumber), vjust = 0.5,hjust=0.1,size=5,colour = "black")+
  labs(title = "KM Sig Gene Number")+
  scale_fill_gradient2()+
  labs(x="Gene Number",fill="",y="Cancer")+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5))+
  mytheme
p
save_plot(p,filename = "TCGA_PAN/Surv/Pan_Cancer_KMSig_GeneNumber_barplot",style = "g")  
fwrite(m,file = "TCGA_PAN/Surv/Sure_Gene_number.csv")

n <- fread("TCGA_PAN/DEG/DEG_number.csv")
n <- top_n(n,15,GeneNumber)
m <- top_n(m,15,GeneNumber)
intersect(m$Cancer,n$Cancer)

###
rm(list = ls())
setwd("D:/WORK/Project/结直肠癌_程序性死亡_模型25/")
library(Seurat)
library(data.table)
library(tidyverse)
library(clustree)
library(ComplexHeatmap)

##空间转录组
library(SpatialExperiment)
library(Matrix)
library(SingleCellExperiment)
library(BayesSpace)
library(tidyverse)
library(data.table)
library(CyDataPro)
rm(list = ls())
gc()
sce <-readVisium("data/ST/CRC//")
sce <- scater::logNormCounts(sce)
sce <- spatialPreprocess(sce, platform="Visium", 
                         n.PCs=7, n.HVGs=2000, log.normalize=FALSE)
# cluster
sce <- qTune(sce, qs=seq(2, 10), platform="Visium", d=7)
qPlot(sce)
sce <- spatialCluster(sce, q=10, platform="Visium", d=7,
                      init.method="mclust", model="t", gamma=2,
                      nrep=10000, burn.in=100,    #  nrep建议大于10000 ，奈何我们为了演示只能牺牲点了
                      save.chain=TRUE)
gc()
sce.enhanced <- spatialEnhance(sce, q=10, platform="Visium", d=7,
                               model="t", gamma=2,
                               jitter_prior=0.3, jitter_scale=3.5,
                               nrep=10000, burn.in=100,  #   为了向内存妥协
                               save.chain=T)
clusterPlot(sce.enhanced,label = "spatial.cluster" )
xx <- fread("DNB/DNB_gene.csv") %>% as.data.frame()
gene_ano <- data.frame(Gene = c(xx$Autophagy,xx$Cuproptosis[1:3],xx$Ferroptosis[1:4],xx$NETosis[1:4],xx$Paraptosis[1:4]),
                       Class = c(rep("Autophagy",39),rep("Cuproptosis",3),rep("Ferroptosis",4),rep("NETosis",4),rep("Paraptosis",4)))
gene_ano<- gene_ano %>% column_to_rownames("Gene") %>% as.data.frame()
inn  <- intersect(rownames(gene_ano),rownames(sce.enhanced))
gene <- gene_ano[inn,] %>% as.data.frame()

markers <- inn
sce.enhanced <- enhanceFeatures(sce.enhanced, sce,
                                feature_names=markers,
                                nrounds=0)
it <- colnames(xx)[5]
for(it in colnames(xx)){
  dir.create(paste0("ST/CRC/",it),recursive = T)
  Gene<- intersect(rownames(sce.enhanced),xx[,it])
  for(i in Gene){
    p <-featurePlot(sce.enhanced, i)
    save_plot(p,filename = paste0("ST/CRC/",it,"/",i),style = "g")
  }
}


p<-clusterPlot(sce.enhanced,label = "spatial.cluster" )
p
save_plot(p,filename = "ST/CRC///Cluster_plot",style = "g")
save(sce.enhanced,file = "data/ST/CRC//sce_enhance.Rdata")
save(sce,file = "data/ST/CRC//sce.Rdata")
dir.create(paste0("ST/PRAD",i))
SpatialFeaturePlot(brac, features = c("PIK3C3", "ATP7A"))

###单细胞
library(Seurat)
count <- Read10X_h5("data/SingCell/UCEC//UCEC_GSE139555_expression.h5")
sce <- CreateSeuratObject(counts = count, min.cells = 3, min.features = 100, project = "CRC")

aa <- fread("data/SingCell/UCEC//UCEC_GSE139555_CellMetainfo_table.tsv",data.table = F)
aa <- aa %>% column_to_rownames("Cell")
aa$`Celltype (minor-lineage)` %>% table()
sce <- AddMetaData(sce,aa)
sce<- ScaleData(sce)
sce<-FindVariableFeatures(sce, selection.method = "vst", nfeatures = 3000) 
gc()
sce <- RunPCA(sce, features =VariableFeatures(sce))
ElbowPlot(sce,ndims = 40)
sce=RunUMAP(sce,  dims = 1:15, 
            reduction = "pca")
sce=RunTSNE(sce,  dims = 1:15, 
            reduction = "pca")
sce = readRDS("data/SingCell/LIHC.RDS")
p_umap=DimPlot(sce, reduction = "umap", group.by = "Celltype..minor.lineage.",label = T,label.size = 5)+mytheme+
  theme(legend.text = element_text(size = 15))+labs(title = "LIHC")
p_umap
dir.create("SingleCell/UCEC",recursive = T)
CyDataPro::save_plot(p_umap,filename = 'SingleCell/LIHC///umap_by_celltype',style = "g")


p_tene=DimPlot(sce, reduction = "tsne", group.by = "Celltype..minor.lineage.",label = T,label.size = 5)+mytheme+
  theme(legend.text = element_text(size = 15))+labs(title = "LIHC")
p_tene
CyDataPro::save_plot(p_tene,filename = 'SingleCell/LIHC////tsne_by_celltype',style = "g")
xx <- fread("DNB/DNB_gene.csv") %>% as.data.frame()

Idents(sce)<- sce$Celltype..minor.lineage.
i=Gene[1]
it <- colnames(xx)[5]
for(it in colnames(xx)){
  dir.create(paste0("SingleCell/UCEC///gene_exp/",it),recursive = T)
  Gene<- intersect(rownames(sce),xx[,it])
  for(i in Gene){
    p<-FeaturePlot(sce, features = i, min.cutoff = "q9",reduction = "tsne",label = T,label.size = 5,label.color = "black")+
      scale_colour_gradient(low = "gray",high = "red")+mytheme
    CyDataPro::save_plot(p,filename = paste0("SingleCell/UCEC//gene_exp/",it,"/",i,"_FeaturePlot_Tsne"),style = "g")
    p<-FeaturePlot(sce, features = i, min.cutoff = "q9",reduction = "umap",label = T,label.size = 5,label.color = "black")+
      scale_colour_gradient(low = "gray",high = "red")+mytheme
    CyDataPro::save_plot(p,filename = paste0("SingleCell/UCEC//gene_exp/",it,"/",i,"_FeaturePlot_Umap"),style = "g")
  }
}





saveRDS(sce,file = "data/SingCell/UCEC.RDS")

data <- read.delim("data/SingCell/UCEC///UCEC_GSE139555_expression_Celltype_minorlineage.txt",row.names = 1)
gene_ano <- data.frame(Gene = c(xx$Autophagy,xx$Cuproptosis[1:3],xx$Ferroptosis[1:4],xx$NETosis[1:4],xx$Paraptosis[1:4]),
                       Class = c(rep("Autophagy",39),rep("Cuproptosis",3),rep("Ferroptosis",4),rep("NETosis",4),rep("Paraptosis",4)))
gene_ano<- gene_ano %>% column_to_rownames("Gene") %>% as.data.frame()
inn  <- intersect(rownames(gene_ano),rownames(data))
gene <- gene_ano[inn,] %>% as.data.frame()
rownames(gene)<-inn
colnames(gene)<-"Class"
exp <- data[inn,]
COL_A<- HeatmapAnnotation(Group = gene$Class,which = "row",show_annotation_name = F,
                          col = list(Class=c("Autophagy"="#32A852","Cuproptosis"="#6DB1EE","Ferroptosis"="#C9244C","NETosis"="#C4A751","Paraptosis"="#87447A")))
dev.off()
dev.new()
Heatmap(exp,cluster_rows = T,cluster_columns = T,use_raster = T,
        row_names_gp = gpar(fontsize=8),column_names_gp = gpar(fontsize=8),
        row_split = gene,
        row_title_gp = gpar(fontsize=10),
        right_annotation = COL_A,show_column_names = F)

dev.new()
library(pheatmap)
coul <- colorRampPalette(c("navy", "white", "firebrick3"))(50)
p<-pheatmap(exp,annotation_row = gene,scale = "row",
         cluster_rows = T,color = coul,border=F,
         fontsize = 10,annotation_names_row = F,
         angle_col = 45,main = "CRC")

CyDataPro::save_plot(p,filename = "SingleCell/UCEC////gene_exp/gene_avrExpre",style = "xx",width = 8,height = 9)
library(monocle3)

data <-GetAssayData(sce, assay = 'RNA', slot = 'counts')
cell_metadata <- sce@meta.data
colnames(cell_metadata)[8] <- "CellType"
gene_annotation <- data.frame(gene_short_name = rownames(data),row.names = row.names(data))
cds <- new_cell_data_set(data,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)

cell_metadata$Celltype..minor.lineage.
cds <- preprocess_cds(cds, num_dim = 15)
plot_pc_variance_explained(cds)

cds <- reduce_dimension(cds,reduction_method = "UMAP")
plot_cells(cds)
dev.new()
p1 <- plot_cells(cds, reduction_method="UMAP", color_cells_by="CellType") + ggtitle('cds.umap')
p1
cds.embed <- cds@int_colData$reducedDims$UMAP
int.embed <- Embeddings(sce, reduction = "umap")
int.embed <- int.embed[rownames(cds.embed),]
cds@int_colData$reducedDims$UMAP <- int.embed
p2 <- plot_cells(cds, reduction_method="UMAP", color_cells_by="CellType") + ggtitle('int.umap')
p2


cds <- cluster_cells(cds,cluster_method = "louvain")
plot_cells(cds, color_cells_by = "CellType")
cds <- learn_graph(cds)
p = plot_cells(cds, color_cells_by = "CellType", label_groups_by_cluster=FALSE,
               label_leaves=FALSE, label_branch_points=FALSE)
p

p=plot_cells(cds,
             color_cells_by = "CellType",
             label_cell_groups=F,
             label_leaves=TRUE,
             label_branch_points=TRUE,
             graph_label_size=3)+mytheme
p
dev.off()
dev.new()
dir.create("SingleCell/UCEC////Monocle")
CyDataPro::save_plot(p,filename = "SingleCell/UCEC////Monocle/time_bin",style = "g")
x <- colnames(colData(cds))
myselect<-function(cds,select.classify,my_select){
  cell_ids<-which(colData(cds)[,select.classify]==my_select)
  closest_vertex<-
    cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex<-as.matrix(closest_vertex[colnames(cds),])
  root_pr_nodes<-
    igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names
                                                              (which.max(table(closest_vertex[cell_ids,]))))]
  root_pr_nodes}
cds<-order_cells(cds,root_pr_nodes=myselect(cds,select.classify='CellType',my_select="Treg"))


p<-plot_cells(cds, color_cells_by = "pseudotime",label_cell_groups=T,label_leaves=T,
           label_branch_points=T,graph_label_size=3)+mytheme
p
CyDataPro::save_plot(p,filename = "SingleCell/UCEC////Monocle/pseudotime",style = "g")
##差异
ciliated_cds_pr_test_res <- graph_test(cds, neighbor_graph="principal_graph", cores=8)
pr_deg_ids <- ciliated_cds_pr_test_res %>% filter(q_value <0.01 & morans_I > 0.05) %>% arrange(-morans_I)

xx <- fread("DNB/DNB_gene.csv") %>% as.data.frame()
it <- colnames(xx)
i=it[5]
g=gene[1]
library(CyDataPro)
for(i in it){
  gene <-intersect(pr_deg_ids$gene_short_name,xx[,i])
  if(length(gene)>0){
    dir.create(paste0("SingleCell/UCEC/////Monocle/",i,"/pseudotimePlot"),recursive = T)
    dir.create(paste0("SingleCell/UCEC/////Monocle/",i,"/cellplot"),recursive = T)
    for(g in gene){
      p<-plot_cells(cds, genes=g,
                 show_trajectory_graph=FALSE,
                 label_cell_groups=F,
                 label_leaves=T)+mytheme
      save_plot(p,filename = paste0("SingleCell/UCEC/////Monocle/",i,"/cellplot/",g),style = "g")
      p<-plot_genes_in_pseudotime(cds[g,],color_cells_by = "CellType",
                               min_expr=0.5, ncol = 1)+mytheme
      save_plot(p,filename = paste0("SingleCell/UCEC/////Monocle/",i,"/pseudotimePlot/",g),style = "g")
    }
  }else{
    print(i)
    next
  }
  
}
fwrite(ciliated_cds_pr_test_res,file = "SingleCell/UCEC/////Monocle/time_DEG_result.csv",row.names = T)

dir.create(paste0("SingleCell/BRCA/gene_exp/",it),recursive = T)




###返修##

#NCOA4与多死亡相关性
setwd("D:/WORK/Project/结直肠癌_程序性死亡_模型25/")
rm(list = ls())
load("../结直肠癌_RNAmSNP_MR_Singlecell16/data/TCGA_CRC.Rdata")
exp <- dectect_log(count_matix)
ssGSEA<- fread("ssGSEA/ssGSEAReslut/ssGSEA_result.csv",data.table = F,check.names = F)
ssGSEA<- ssGSEA %>% column_to_rownames("V1")

inn<- intersect(rownames(ssGSEA),colnames(exp))

exp_G<- exp["NCOA4",inn] %>% t()
ssGSEA<- ssGSEA[inn,]

library(psych)

re<- corr.test(as.numeric(exp_G),ssGSEA)
R<- re$r %>% as.data.frame()
P<- re$p.adj %>% as.data.frame()
result<- rbind(R,P) %>% t() %>% as.data.frame()
colnames(result)<- c("R2","Padj")
result$PDC<- rownames(result)

dir.create("返修1/NCOA4与PCD相关性",recursive = T)
fwrite(result,file = "返修1/NCOA4与PCD相关性/NCOA4_PCD_ssGSEA_corr_result.csv")

result$Gene<- 1
result$sig<- ifelse(result$Padj<0.001,"***",
                    ifelse(result$Padj<0.01,"**",
                           ifelse(result$Padj<0.05,"*"," ")))
result$label<- paste0(round(result$R2,3),"\n",result$sig)

p<-ggplot(result, aes(x=reorder(PDC,R2), y=Gene, fill=R2)) +
  geom_tile(colour="white") +
  
  scale_fill_gradient(low = "white", high = "#c85c32")+
  coord_polar(theta = "x", start = pi/3)+
  geom_text(aes(label=label,y=Gene+0.2),size=6)+
  theme_void()+labs(x="",y="")+
  theme(axis.text.y = element_blank(),
        legend.position = "bottom",
        legend.direction = "horizontal",
        axis.text = element_text(size = 15,face = "bold",colour = "black"))
p
save_plot(p,filename = "返修1/NCOA4与PCD相关性/NCOA4_PCD_ssGSEA_corr_plot",style = "g",height = 10,width = 12)

plotd<- cbind(exp_G,ssGSEA)
dir.create("返修1/NCOA4与PCD相关性/CorrPlot")
library(ggpubr)
for(i in colnames(plotd)[-1]){
  p<-ggscatter(plotd,x="NCOA4",y=i,
               size = 1.5,
               add = "reg.line",
               add.params = list(color = "#c85c32", fill = "#FEE08B", size = 1),
               conf.int = TRUE)+
    labs(x="NCOA4 Expression Level",y=paste0(i," ssGSEA Score"))+
    stat_cor(method = "pearson", label.sep = "\n",size=8)+theme_test(base_rect_size = 1.5)+mytheme
  save_plot(p,filename = paste0("返修1/NCOA4与PCD相关性/CorrPlot/",i),style = "g")
}


##功能分析
library(limma)
library(clusterProfiler)
library(org.Hs.eg.db)
gene="NCOA4"
expr=exp
Pvalue=0.05;logFC_cutoff=0.7
dir.create("返修1/功能富集")
Out_dir="返修1/功能富集/"
group<- ifelse(as.numeric(expr[gene,])> median(as.numeric(expr[gene,])),"High","Low")
group <- factor(group,levels = c("Low","High"))
exp<- normalizeBetweenArrays(expr)
design=model.matrix(~group)
fit=lmFit(exp,design)
fit=eBayes(fit)
deg=topTable(fit,coef=2,number = Inf)
K1 <- (deg$P.Value<Pvalue)&(deg$logFC < -logFC_cutoff)
K2 <- (deg$P.Value<Pvalue)&(deg$logFC > logFC_cutoff)
deg$change <- ifelse(K1,"Down",ifelse(K2,"Up","not"))
gg<- deg %>% filter(change !="not")

Erich_GO_KEGG(rownames(gg),out = paste0(Out_dir,"/",gene),top = 15,plotStyle = "bar",method = "GO",width = 10,height = 8,color = c("#2166AC","#B2182B"))
Erich_GO_KEGG(rownames(gg),out = paste0(Out_dir,"/",gene),top = 15,plotStyle = "bar",method = "KEGG",width = 10,height = 8,color = c("#2166AC","#B2182B"))
colnames(plotd)
go<- fread(paste0(Out_dir,"/",gene,"_GO_result.txt"),data.table = F)
item1<- c("Apoptosis","apoptotic","cell death","Entosis",
          "Pyroptosis","autophagy","p53","Anoikis")
enrich_all<- go %>% filter(pvalue < 0.05,qvalue<0.2)
all<- data.frame()

for(i in item1){
  tmp<-enrich_all %>% filter(str_detect(enrich_all$Description,i))
  if(nrow(tmp)==0){
    next
  }
  
  all<- rbind(tmp,all)
}
fwrite(all,file = paste0(Out_dir,"/",gene,"_GO_细胞死亡相关.csv"))
plot<- all %>% arrange(pvalue) %>% head(10)
plot<- plot %>% arrange(Count)
plot$Description<-factor(plot$Description,levels = plot$Description)
library(ggthemes)
p<-ggplot(plot,aes(Count,Description))+
  geom_point(aes(size=Count,color=-log10(pvalue)),shape=16)+
  scale_color_gradient(high = "#c85c32",low = "#A6D96A")+
  labs(x="Gene Counts",y="",title = "GO enrichment")+
  theme_clean()+
  theme(panel.grid = element_blank())+
  theme(axis.title = element_text(size = 15,face = "bold"),
        axis.text = element_text(size = 15,face = "bold",colour = "black"),
        text = element_text(size = 15,colour = "black",face = "bold"))

p
CyDataPro::save_plot(p,filename = paste0(Out_dir,"/",gene,"_GO_细胞死亡相关"),style = "g",width = 12)


item2<- unique(c("Cancer","Chemical Carcinogenesis","PI3K-Akt","Focal Adhesion",
                 "Apoptosis","apoptotic","Cell death","Entosis",
                 "Pyroptosis","autophagy","p53","Anoikis",
                 "Epithelial to Mesenchymal Transition","ECM-","Cell Cycle",
                 "DNA Replication","p53","PI3K-Akt","MAPK","Wnt","Hedgehog","Notch",
                 "TGF-beta","Apoptosis","FoxO","PI3K-AKT","mTOR","AMPK",
                 "apoptotic","FAS-JNK","TNF","NF-kB","Cell cycle",
                 "DNA repair","Cancer pathways","Pathways in cancer","Leukocyte transendothelial migration",
                 "Cytokine-cytokine","Hematopoietic cell lineage"))
kegg<- fread(paste0(Out_dir,"/",gene,"_KEGG_result.txt"),data.table = F)

enrich_all<- kegg %>% filter(pvalue < 0.05,qvalue<0.2)
all<- data.frame()
for(i in unique(item2)){
  tmp<-enrich_all %>% filter(str_detect(enrich_all$Description,i))
  if(nrow(tmp)==0){
    next
  }
  all<- rbind(tmp,all)
}

fwrite(all,file = paste0(Out_dir,"/",gene,"_KEGG_细胞死亡相关.csv"))
plot<- all %>% arrange(pvalue) %>% head(10)
plot<- plot %>% arrange(Count)
plot$Description<-factor(plot$Description,levels = plot$Description)

p<-ggplot(plot,aes(Count,Description))+
  geom_point(aes(size=Count,color=-log10(pvalue)),shape=16)+
  scale_color_gradient(high = "#B2182B",low = "#A6D96A")+
  labs(x="Gene Counts",y="",title = "KEGG enrichment")+
  theme_clean()+
  theme(panel.grid = element_blank())+
  theme(axis.title = element_text(size = 15,face = "bold"),
        axis.text = element_text(size = 15,face = "bold",colour = "black"),
        text = element_text(size = 15,colour = "black",face = "bold"))
p
CyDataPro::save_plot(p,filename = paste0(Out_dir,"/",gene,"_KEGG_细胞死亡相关"),style = "g",width = 12)
