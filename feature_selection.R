# library(tidyverse)
library(mixOmics)
library(e1071)
library(ROCR)
library(RColorBrewer)
library(ggrepel)
# library(caret)

#################################
# Set universal parameters
#################################
is.residual = TRUE

optimal_comp = FALSE
components = 5
vip_threshold = 2
foldchange_threshold = 0
is.log = TRUE
fold_cv = 10
nrepeat = 50
pred.eval.method = "BER"
cluster.method = "dist" #"dist" or "bicor" #Residual can only use "dist"
cpu = 4

dir.folder <- "C:/Users/QiYan/Dropbox/AIME/Pesticide/HILICpos/HILIC_OPs/"

#################################
# Start!
#################################

time.start <- proc.time()

# Set log and graph direction
dir.temp <- paste(dir.folder,"PANDA_output_PLSDA",sep = "")
setwd(dir.temp)

if(is.residual == TRUE){
  logfile <- paste("Res_Log_",Sys.Date(),"_vip",vip_threshold,"fc",foldchange_threshold,".log",sep = "")
  pdf_file<-paste("Res_Output_",Sys.Date(),"_vip",vip_threshold,"fc",foldchange_threshold,".pdf",sep="")
}else{
  logfile <- paste("Log_",Sys.Date(),"_vip",vip_threshold,"fc",foldchange_threshold,".log",sep = "")
  pdf_file<-paste("Output_",Sys.Date(),"_vip",vip_threshold,"fc",foldchange_threshold,".pdf",sep="")
}

logfile <- file(logfile)
sink(logfile) # sink()

pdf(file = pdf_file,width=10,height=10)

# Load data
print("Good Luck, Have Fun!")
print("Read metabolomics data")
print("After preprocess")

dir.temp <- paste(dir.folder,"PANDA_input",sep = "")
setwd(dir.temp)

{

# For fold change

load(file = "HILIC_OPs_transformed.RData")

X <- final_feature
sampleID <- final_sample
row.names(sampleID) <- sampleID$SampleID
sampleID$factorcase <- as.numeric(ifelse(sampleID$factorcase == "Exposed", 1, ifelse(sampleID$factorcase == "Unexposed", 0, 99)))
sampleID$factorcase <- as.factor(sampleID$factorcase)
Y <- sampleID$factorcase
levels(Y) = c("control","case")

X = X[ ,apply(X[,1:ncol(X)],2,var) != 0] # remove 0 variance variables

fc.X_caret <- as.data.frame(cbind(X,Y))
fc.X_caret$Y <- as.factor(fc.X_caret$Y)
levels(fc.X_caret$Y) = c("control","case")

# For the rest
load(file = "HILIC_OPs_redisual.RData")

X <- final_feature
sampleID <- final_sample
row.names(sampleID) <- sampleID$SampleID
sampleID$factorcase <- as.numeric(ifelse(sampleID$factorcase == "Exposed", 1, ifelse(sampleID$factorcase == "Unexposed", 0, 99)))
sampleID$factorcase <- as.factor(sampleID$factorcase)
Y <- sampleID$factorcase
levels(Y) = c("control","case")

X = X[ ,apply(X[,1:ncol(X)],2,var) != 0] # remove 0 variance variables

X_caret <- as.data.frame(cbind(X,Y))
X_caret$Y <- as.factor(X_caret$Y)
levels(X_caret$Y) = c("control","case")

}

###################################################
### I. PCA
###################################################

# PCA
print("PCA analysis")
pca.datExpr = pca(X, ncomp = 10, center = TRUE, scale = TRUE)
plot(pca.datExpr,main = "Variance explained by PCs using all features")
plotIndiv(pca.datExpr, group = Y, ind.names = FALSE, 
          legend = TRUE, title = 'PC score using all features after preprocessing')
print("PCA done")

###################################################
### II. PLSDA
###################################################

# PLSDA
print("PLS-DA analysis")
datExpr.plsda <- plsda(X, Y, ncomp = 10) # set ncomp to 10 for performance assessment later
plotIndiv(datExpr.plsda, comp = 1:2,
          group = Y, ind.names = FALSE, 
          ellipse = TRUE, legend = TRUE, title = 'PLS score using all features after preprocessing')

# Plot PLS1 vs. PLS2 with background
background = background.predict(datExpr.plsda, comp.predicted=2, dist = "max.dist") 
plotIndiv(datExpr.plsda, comp = 1:2,
          group = Y, ind.names = FALSE, title = "PLS score using all features after preprocessing (Maximum distance)", ellipse = TRUE,
          legend = TRUE,  background = background)

# Assess the preformance of PLSDA, tuning, select number of component
print("Optimize number of component if necessary")
if (optimal_comp == TRUE){
set.seed(1106)
system.time(
perf.plsda.datExpr <- perf(datExpr.plsda, validation = "Mfold", folds = fold_cv,   ## 5 fold CV
                  progressBar = TRUE, auc = TRUE, nrepeat = nrepeat, cpus = cpu)
)
plot(perf.plsda.datExpr, col = color.mixo(5:7), sd = TRUE, legend.position = "horizontal")

error_rate <- as.data.frame(perf.plsda.datExpr$error.rate$BER)
num_comp <- rownames(error_rate)[which.min(apply(error_rate,MARGIN=1,min))]
dist_method <- as.numeric(which.min(apply(error_rate,MARGIN=2,min)))
# if (optimal_comp == TRUE){
  opt_comp = as.numeric(perf.plsda.datExpr$choice.ncomp[2,dist_method])
}else{
  opt_comp = components
}

# select important variables
print("Select important variables")
vip.plsda.datExpr <- vip(datExpr.plsda)
good_feats<-{} 
for(c1 in 1:opt_comp){
  good_feats<-c(good_feats,which(vip.plsda.datExpr[,c1]>=vip_threshold))
}
good_feats<-unique(good_feats)

# fold change
print("Calculate fold change")

fc.empty <- data.frame()
for(i in 1:(ncol(fc.X_caret)-1)){
  case.temp <- fc.X_caret[which(as.character(fc.X_caret$Y) == "case"),i]
  case <- mean(case.temp[!(abs(case.temp - median(case.temp)) > 2*sd(case.temp))])
  control.temp <- fc.X_caret[which(as.character(fc.X_caret$Y) == "control"),i]
  control <- mean(control.temp[!(abs(control.temp - median(control.temp)) > 2*sd(control.temp))])
  fc.empty[i,1] <- case
  fc.empty[i,2] <- control
}
colnames(fc.empty) <- c("case","control")
fc.empty$foldchange = fc.empty$case-fc.empty$control
fc.updown <- fc.empty
rm(fc.empty,case,case.temp,control, control.temp)

# get vip score for all expr features and plot manhattan plots
print("Calculate vip score and plot manhattan plots")
vip.for.selection<-as.data.frame(apply(vip.plsda.datExpr[,c(1:opt_comp)],1,max))
vip.for.selection <- cbind(vip.for.selection,fc.updown$foldchange)
colnames(vip.for.selection) = c("vip","foldchange")

vip.for.selection <- cbind(linkid,vip.for.selection)
vip.for.selection$group[vip.for.selection$vip < vip_threshold] <- 0
vip.for.selection$group[vip.for.selection$vip >= vip_threshold & vip.for.selection$foldchange <= 0] <- 1 ## down
vip.for.selection$group[vip.for.selection$vip >= vip_threshold & vip.for.selection$foldchange > 0] <- 2  ## up

ggplot2::ggplot(vip.for.selection,aes(x=mz,y=vip)) +
  geom_point(aes(colour = cut(group, c(-Inf,0,1,2,Inf))),size=1,show.legend = FALSE) + 
  scale_fill_hue(c=20, l=20) + 
  scale_color_manual(values = c("#999999","springgreen3","firebrick1")) +
  geom_hline(aes(yintercept = 2),color = "red",size = 0.5,linetype = "dashed") +
  scale_x_continuous(breaks = pretty(vip.for.selection$mz, n = 10))+
  scale_y_continuous(breaks = pretty(vip.for.selection$vip, n = 10))+
  # theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme_bw(base_size = 12)+
  theme(axis.text=element_text(face="bold"),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        plot.title = element_text(size=12, face="bold"),
        axis.title.x = element_text(size=12, face="bold"),
        axis.title.y = element_text(size=12, face="bold"))+
  labs(x="M/Z",y="VIP Score") +
  ggtitle("Type 1 manhattan plot (VIP vs mz)")

ggplot2::ggplot(vip.for.selection,aes(x=time,y=vip)) +
  geom_point(aes(colour = cut(group, c(-Inf,0,1,2,Inf))),size=1,show.legend = FALSE) + 
  scale_fill_hue(c=20, l=20) + 
  scale_color_manual(values = c("#999999","springgreen3","firebrick1")) +
  geom_hline(aes(yintercept = 2),color = "red",size = 0.5,linetype = "dashed") +
  scale_x_continuous(breaks = pretty(vip.for.selection$time, n = 10))+
  scale_y_continuous(breaks = pretty(vip.for.selection$vip, n = 10))+
  # theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme_bw(base_size = 12)+
  theme(axis.text=element_text(face="bold"),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        plot.title = element_text(size=12, face="bold"),
        axis.title.x = element_text(size=12, face="bold"),
        axis.title.y = element_text(size=12, face="bold"))+
  labs(x="Retention Time",y="VIP Score") +
  ggtitle("Type 2 manhattan plot (VIP vs retention time)")

volcano <- vip.for.selection[,c(1:4)]
volcano$group <- ifelse(volcano$vip>=vip_threshold&volcano$foldchange>=1,3,
                        ifelse(volcano$vip>=vip_threshold&volcano$foldchange<=-1,2,
                               ifelse(volcano$vip>=vip_threshold&abs(volcano$foldchange)<=1,1,0)))
volcano$size <- ifelse(volcano$group == 0,0,1)

ggplot2::ggplot(volcano,aes(x=foldchange,y=vip)) +
  geom_point(aes(colour=cut(group, c(-Inf,0,1,2,Inf))),show.legend = FALSE) + 
  xlim(-3.8,3.8) +
  scale_fill_hue(c=20, l=20) + 
  scale_color_manual(values = c("#999999","goldenrod3","springgreen3","firebrick1")) + 
  scale_size_continuous(range = c(1,3)) +
  geom_hline(aes(yintercept = 2),color = "black",size = 0.5,linetype = "dashed") +
  geom_vline(aes(xintercept = 1),color = "black",size = 0.5,linetype = "dashed") +
  geom_vline(aes(xintercept = -1),color = "black",size = 0.5,linetype = "dashed") +
  labs(x="Log2(Fold Change)",y="VIP Score") +
  # theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme_bw(base_size = 12)+
  theme(axis.text=element_text(face="bold"),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        plot.title = element_text(size=12, face="bold"),
        axis.title.x = element_text(size=12, face="bold"),
        axis.title.y = element_text(size=12, face="bold"))+
  geom_label_repel(data = subset(volcano,volcano$vip>=3&abs(volcano$foldchange)>=1),
                   aes(label = paste("mz:",round(mz,4),sep = "")))+
  ggtitle("Volcano Plot")

print("PLS-DA phase 1 done")

# Using significant features
print(paste("Assess performance of PLS-DA using significant features, threshold vip score >=",vip_threshold,sep = " "))
good_feats_ordered <- vip.for.selection[order(-vip.for.selection$vip,-abs(vip.for.selection$foldchange)),]
good_feats_ordered <- subset(good_feats_ordered,vip>=vip_threshold & abs(foldchange)>=foldchange_threshold)
good_feats_ordered.name <- {}
for (i in 1:nrow(good_feats_ordered)) {
  good_feats_ordered.name <- c(good_feats_ordered.name,which(colnames(X)==row.names(good_feats_ordered[i,])))
}
sig.X <- X[,c(good_feats_ordered.name)]

#pca using significant features
print("PCA using significant features")
sig.pca.datExpr = pca(sig.X, ncomp = 10, center = TRUE, scale = TRUE)
plot(sig.pca.datExpr,main = "Variance explained by PCs using significant features")
plotIndiv(sig.pca.datExpr, group = Y, ind.names = FALSE, 
          legend = TRUE, title = 'PC score using significant features')

# PLSDA using significant features
print("PLS-DA using significant features")
sig.datExpr.plsda <- plsda(sig.X, Y, ncomp = 10) # set ncomp to 10 for performance assessment later
plotIndiv(sig.datExpr.plsda, comp = 1:2,
          group = Y, ind.names = FALSE, 
          ellipse = TRUE, legend = TRUE, title = 'PLS score using significant features')

# Plot PLS1 vs. PLS2 with background using significant features
background = background.predict(sig.datExpr.plsda, comp.predicted=2, dist = "max.dist")
plotIndiv(sig.datExpr.plsda, comp = 1:2,
          group = Y, ind.names = FALSE, title = "PLS score using significant features (Maximum distance)", ellipse = TRUE,
          legend = TRUE,  background = background)

# Assess the preformance of PLSDA using significant features
print("Assess the preformance")
set.seed(1106)
system.time(
  sig.perf.plsda.datExpr <- perf(sig.datExpr.plsda, validation = "Mfold", folds = fold_cv,   ## 5 fold CV
                           progressBar = TRUE, auc = TRUE, nrepeat = nrepeat, cpus = cpu)
)

# get ROC curve using significant features
print("Generating ROC curve using top features on training set")
source("C:/Users/QiYan/Dropbox/AIME/Archive/get_roc.R")
roc.dataA <- t(sig.X)
get_roc(dataA=roc.dataA,classlabels=Y,classifier="svm",kname="radial",
        rocfeatlist=seq(2,10,1),rocfeatincrement=TRUE,mainlabel="Training set ROC curve using top features")
print("ROC done")

# get CV accuracy using significant features via SVM
print("get k-fold CV accuracy")
source("C:/Users/QiYan/Dropbox/AIME/Archive/svm_cv.R")
xvec<-{}
yvec<-{}
best_acc<-0
for(i in 2:nrow(roc.dataA)){
  subdata<-t(roc.dataA[1:i,])
  svm_model<-try(svm_cv(v=fold_cv,x=subdata,y=Y,kname="radial",
                        errortype=pred.eval.method,conflevel=95))
  if(is(svm_model,"try-error")){
    svm_model<-NA
    }else{
      xvec<-c(xvec,i)
      yvec<-c(yvec,svm_model$avg_acc)
      if(svm_model$avg_acc>best_acc){
      best_acc<-svm_model$avg_acc
      best_subset<-seq(1,i)
      }
      if(svm_model$avg_acc<best_acc){
        diff_acc<-best_acc-svm_model$avg_acc
        if(diff_acc>50){
          break;
        }
      }
    }
  }

if(pred.eval.method=="CV"){
  ylab_text=paste(pred.eval.method," accuracy (%)",sep="")
  }else{
    if(pred.eval.method=="BER"){
      ylab_text=paste("Balanced accuracy"," (%)",sep="")
      }else{
        ylab_text=paste("AUC"," (%)",sep="")
      }
    }
print(length(yvec))
print(length(xvec))
print(xvec)
print(yvec)
if(length(yvec)>0){
  plot(x=xvec,y=yvec,main="k-fold CV classification accuracy based on forward selection of top features",xlab="Feature index",ylab=ylab_text,type="b",col="brown")
  
  cv_mat<-cbind(xvec,yvec)
  colnames(cv_mat)<-c("Feature Index",ylab_text)
  
  # write.table(cv_mat,file="kfold_cv_mat.txt",sep="\t")
}

cv.acc.sigfeats <- svm_model$avg_acc
print(paste(paste("K-fold CV accuracy is",cv.acc.sigfeats,sep = " "),"%",sep = ""))

#permutation test
print("Calculating permuted CV accuracy")
cv.acc.permut<-{}
subdata <- t(roc.dataA)
cv.acc.permut<-lapply(1:100,function(j){
  rand_order<-sample(1:dim(as.data.frame(Y))[1],size=dim(as.data.frame(Y))[1])
  classlabels_permut<-as.data.frame(Y)[rand_order,]
  classlabels_permut<-as.data.frame(classlabels_permut)
  svm_permut_res<-svm_cv(v=fold_cv,x=subdata,y=classlabels_permut,kname="radial",errortype=pred.eval.method,conflevel=95)
  return(svm_permut_res$avg_acc)
})

cv.acc.permut<-unlist(cv.acc.permut)
cv.acc.permut<-mean(cv.acc.permut,na.rm=TRUE)
cv.acc.permut<-round(cv.acc.permut,2)

print(paste(paste("mean Permuted accuracy is:",cv.acc.permut,sep = " "),"%",sep = ""))

# Plot var plot for permutation test
bar <- as.data.frame(c(cv.acc.sigfeats,cv.acc.permut))
bar$name <- c("Actual","Permutated")
colnames(bar) <- c("value","name")
bar$value <- round(bar$value,2)
ggplot2::ggplot(bar, aes(x=name, y=value)) +
  geom_bar(stat = "identity",aes(fill = as.factor(name)),width = 0.5,show.legend = FALSE) +
  geom_text(aes(label=value),vjust=-0.3) +
  labs(title = "K-fold CV accuracy and permutation test",x = element_blank(),y = "Accuracy (%)")

print("PLS-DA phase 2 done")

# Two-way HCA
print("Two-way HCA")
source("C:/Users/QiYan/Dropbox/AIME/Archive/get_hca.R")
get_hca(data_m = sig.X,classlabels = Y,is.data.znorm = FALSE,clu.method = cluster.method)

# Save files
save.plsresults.allfeatures <- cbind(vip.for.selection,t(X))
good_feats_ordered$rank <- 1:nrow(good_feats_ordered)
save.plsresults.sigfeatures <- cbind(good_feats_ordered,roc.dataA)
save.cv.accuracy <- as.data.frame(cv_mat)

print(paste("Number of significant features:",nrow(save.plsresults.sigfeatures),sep = " "))

###################################################
### III. Prepare for mummichog
###################################################

save.mummichog_PLSDA_VIP2 <- save.plsresults.allfeatures[,1:4]
save.mummichog_PLSDA_VIP2$"p-value" = 0.051
save.mummichog_PLSDA_VIP2$`p-value`[save.mummichog_PLSDA_VIP2$vip>=vip_threshold&abs(save.mummichog_PLSDA_VIP2$foldchange)>=foldchange_threshold] <- 0.04
save.mummichog_PLSDA_VIP2 <- save.mummichog_PLSDA_VIP2[order(save.mummichog_PLSDA_VIP2$`p-value`,-save.mummichog_PLSDA_VIP2$vip),]
save.mummichog_PLSDA_VIP2 <- save.mummichog_PLSDA_VIP2[,c(1,2,5,4)]

# ###################################################
# ### IV. Annotation
# ###################################################
# 
# #PLSDA_vip2
# annotations_filename<-"C:/Users/QiYan/Dropbox/AIME/PNS_Ritz/HILICpos_ThermoHFQE_85to1275_mz_range/DBmatches_KEGG.txt"
# annotations_filename_2<-"C:/Users/QiYan/Dropbox/AIME/PNS_Ritz/HILICpos_ThermoHFQE_85to1275_mz_range/DBmatches_HMDB.txt"
# annotations_filename_3<-"C:/Users/QiYan/Dropbox/AIME/PNS_Ritz/HILICpos_ThermoHFQE_85to1275_mz_range/DBmatches_LipidMaps.txt"
# 
# number_significant_digits_rounding<-4
# 
# d_sigfeature<-save.plsresults.sigfeatures
# dim(d_sigfeature)
# 
# a_KEGG<-read.table(annotations_filename,sep="\t",header=TRUE)
# a_HMDB<-read.table(annotations_filename_2,sep="\t",header=TRUE)
# a_LipidMaps<-read.table(annotations_filename_3,sep="\t",header=TRUE)
# 
# a_KEGG$mz<-round(a_KEGG$mz,number_significant_digits_rounding)
# a_KEGG<-a_KEGG[,-c(20:2032)]
# a_HMDB$mz<-round(a_HMDB$mz,number_significant_digits_rounding)
# a_HMDB<-a_HMDB[,-c(20:2032)]
# a_LipidMaps$mz<-round(a_LipidMaps$mz,number_significant_digits_rounding)
# a_LipidMaps<-a_LipidMaps[,-c(20:2032)]
# 
# a_KEGG$time<-round(a_KEGG$time,2)
# a_KEGG<-a_KEGG[,-c(20:2032)]
# a_HMDB$time<-round(a_HMDB$time,2)
# a_HMDB<-a_HMDB[,-c(20:2032)]
# a_LipidMaps$time<-round(a_LipidMaps$time,2)
# a_LipidMaps<-a_LipidMaps[,-c(20:2032)]
# 
# d_sigfeature$mz<-round(d_sigfeature$mz,number_significant_digits_rounding)
# d_sigfeature$time<-round(d_sigfeature$time,2)
# 
# m_KEGG<-merge(a_KEGG,d_sigfeature,by.x=c("mz","time"),by.y=c("mz","time"))      ## Merged data
# m_KEGG<-m_KEGG[order(m_KEGG$rank),]
# m_HMDB<-merge(a_HMDB,d_sigfeature,by.x=c("mz","time"),by.y=c("mz","time"))
# m_HMDB<-m_HMDB[order(m_HMDB$rank),]
# m_LipidMaps<-merge(a_LipidMaps,d_sigfeature,by.x=c("mz","time"),by.y=c("mz","time"))
# m_LipidMaps<-m_LipidMaps[order(m_LipidMaps$rank),]
# 
# save.KEGG<-subset(m_KEGG,select=c("mz","time","KEGGID","Name","rank"))
# names(save.KEGG)[names(save.KEGG)=="Name"] <- "KEGG_name"
# save.HMDB<-subset(m_HMDB,select=c("mz","time","HMDBID","Name","rank"))
# names(save.HMDB)[names(save.HMDB)=="Name"] <- "HMDB_name"
# save.LipidMaps<-subset(m_LipidMaps,select=c("mz","time","LM_ID","Name","rank"))
# names(save.LipidMaps)[names(save.LipidMaps)=="Name"] <- "LipidMaps_name"
# 
# # save.all<-dplyr::full_join(save.HMDB,save.KEGG,by=c('mz','rank'))
# # save.all<-dplyr::full_join(save.all,save.LipidMaps,by=c('mz','rank'))
# 
# rm(a_HMDB,a_KEGG,a_LipidMaps,m_HMDB,m_KEGG,m_LipidMaps)

###################################################
### IV. Save files
###################################################

dir.temp <- paste(dir.folder,"PANDA_output_PLSDA",sep = "")
setwd(dir.temp)

filename <- paste("Res_PLSDA_result_",Sys.Date(),"_vip",vip_threshold,"fc",foldchange_threshold,".RData",sep="")

save(save.cv.accuracy,save.mummichog_PLSDA_VIP2,save.plsresults.allfeatures,
     save.plsresults.sigfeatures,file = filename)

time.end <- proc.time()
time.run <- time.end - time.start

print("The running time is")
print(time.run)

dev.off()
sink()

# Create mummichog files

load(filename)
mummichog_name <- paste("C:/Users/QiYan/Dropbox/AIME/Pesticide/HILICpos/HILIC_OPs/HILIC_mummichog/mummichog_input","_vip",vip_threshold,"fc",foldchange_threshold,"_",Sys.Date(),".txt",sep="")
write.table(save.mummichog_PLSDA_VIP2,file=mummichog_name,sep = "\t",row.names = F,quote = F)