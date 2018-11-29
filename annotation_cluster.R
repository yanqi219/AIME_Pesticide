###################################################
### VIII. Re-annotate significant features with loose adducts
###################################################

# Install packages

# source("http://bioconductor.org/biocLite.R")
# 
# install.packages("snow",repos="http://cran.r-project.org")
# install.packages("doSNOW",repos="http://cran.r-project.org")
# install.packages("parallel",repos="http://cran.r-project.org")
# install.packages("e1071",repos="http://cran.r-project.org")
# install.packages("XML",repos="http://cran.r-project.org")
# install.packages("R2HTML",repos="http://cran.r-project.org")
# install.packages("RCurl",repos="http://cran.r-project.org")
# source("http://bioconductor.org/biocLite.R")
# biocLite("Rdisop",suppressUpdates=TRUE)
# biocLite("KEGGREST",suppressUpdates=TRUE)
# biocLite("pcaMethods",suppressUpdates=TRUE)
# install.packages("flashClust",repos="http://cran.r-project.org")
# install.packages("plyr",repos="http://cran.r-project.org")
# install.packages("png",repos="http://cran.r-project.org")
# install.packages("rjson",repos="http://cran.r-project.org")
# 
# install.packages("/u/home/q/qyan/AIME/r_package/SSOAP_0.9-0.tar.gz", repos = NULL, type = "source")
# install.packages("/u/home/q/qyan/AIME/r_package/XMLSchema_0.7-0.tar.gz", repos = NULL, type = "source")
# install.packages("/u/home/q/qyan/AIME/r_package/xMSannotator_1.3.2.tar.gz", repos = NULL, type = "source")

library(xMSannotator)

annotation <- function(wd, data, outloc, database, mode){
  # setwd("C:/Users/QiYanYan/Dropbox/AIME/Panda_HILICpos/HILIC_Controls_ExpoUnexpo/PANDA_output_PLSDA")
  # load(file = "Res_PLSDA_result_2018-05-31_vip2fc0.RData")
  setwd(wd)
  # setwd("/u/home/q/qyan/AIME/Annotation")
  # load(file = "HILIC_control_expo_unexpo_residual_nonorm_WGCNA.RData")
  load(file = data)
  
  data(adduct_table)
  data(adduct_weights)
  
  normalization <- function(x){
    return((x-mean(x))/(max(x)-min(x)))
  }
  
  ###########Parameters to change##############
  
  # input_data <- save.plsresults.sigfeatures[,-c(3:6)]
  # input_data.norm <- cbind(input_data[,c(1,2)],apply(input_data[,-c(1,2)],2,normalization))
  # input_data <- wide_save_residual
  input_data <- as.data.frame(t(after.prepro.feature))
  input_data <- cbind(after.prepro.linkid,input_data)
  dataA<-input_data
  
  outloc<-outloc
  # outloc<-"/u/home/q/qyan/AIME/Annotation/HILIC"
  
  max.mz.diff<-10  #mass search tolerance for DB matching in ppm
  max.rt.diff<-10 #retention time tolerance between adducts/isotopes
  corthresh<-0.7 #correlation threshold between adducts/isotopes
  max_isp=5 #maximum number of isotopes to search for
  mass_defect_window=0.01 #mass defect window for isotope search
  
  num_nodes<-4   #number of cores to be used; 2 is recommended for desktop computers due to high memory consumption
  
  db_name=database #other options: HMDB,Custom,KEGG, LipidMaps, T3DB
  status="Detected and Quantified" #other options: "Detected", NA, "Detected and Quantified", "Expected and Not Quantified"
  num_sets<-300 #number of sets into which the total number of database entries should be split into;
  
  mode<-mode #ionization mode
  # queryadductlist=c("M+H","M+2H","M+H+NH4","M+ACN+2H","M+2ACN+2H","M+NH4","M+Na","M+ACN+H","M+ACN+Na","M+2ACN+H","2M+H","2M+Na","2M+ACN+H","M+2Na-H","M+H-H2O","M+H-2H2O") #other options: c("M-H","M-H2O-H","M+Na-2H","M+Cl","M+FA-H"); c("positive"); c("negative"); c("all");see data(adduct_table) for complete list
  queryadductlist=c("all")
  if(mode == "pos"){
    filter = "M+H"
  }else{
    filter = "M-H"
  }
  
  #provide list of database IDs (depending upon selected database) for annotating only specific metabolites
  customIDs<-NA #c("HMDB15352","HMDB60549","HMDB00159","HMDB00222"); read.csv("/Users/mzmatch_95stdmx_HMDBIDs.csv")
  customDB<-NA 
  
  #########################
  
  dataA<-unique(dataA)
  print(dim(dataA))
  print(format(Sys.time(), "%a %b %d %X %Y"))
  
  system.time(annotres<-multilevelannotation(dataA=dataA,max.mz.diff=max.mz.diff,max.rt.diff=max.rt.diff,cormethod="pearson",
                                             num_nodes=num_nodes,queryadductlist=queryadductlist,
                                             mode=mode,outloc=outloc,db_name=db_name, adduct_weights=adduct_weights,num_sets=num_sets,allsteps=TRUE,
                                             corthresh=corthresh,NOPS_check=TRUE,customIDs=customIDs,missing.value=NA,
                                             deepsplit=2,networktype="unsigned",minclustsize=10,module.merge.dissimilarity=0.2,filter.by=c(filter),
                                             biofluid.location="Blood",origin=NA,status=status,boostIDs=NA,max_isp=max_isp,
                                             customDB=customDB,MplusH.abundance.ratio.check = FALSE, min_ions_perchem = 1,
                                             HMDBselect="all",mass_defect_window=mass_defect_window,pathwaycheckmode="pm",mass_defect_mode="both")
  )
  
  
  print(format(Sys.time(), "%a %b %d %X %Y"))
  
  # pkg <- "package:xMSannotator"
  # detach(pkg, character.only = TRUE)
}

annotation(wd = "C:/Users/QiYanYan/Dropbox/AIME/Panda_C18neg/C18_Controls_ExpoUnexpo/PANDA_input",
           data = "C18_control_expo_unexpo_classification_nonorm.RData",
           outloc = "C:/Users/QiYanYan/Dropbox/AIME/Panda_C18neg/C18_Controls_ExpoUnexpo/C18_Annotation/sigfeature_annotation/KEGG",
           database = "KEGG",
           mode = "neg")

annotation(wd = "C:/Users/QiYanYan/Dropbox/AIME/Panda_C18neg/C18_Non_Exposed_CasesControls/PANDA_input",
           data = "C18_case_control_noexposure_classification_nonorm.RData",
           outloc = "C:/Users/QiYanYan/Dropbox/AIME/Panda_C18neg/C18_Non_Exposed_CasesControls/C18_Annotation/sigfeature_annotation/HMDB",
           database = "HMDB",
           mode = "neg")
annotation(wd = "C:/Users/QiYanYan/Dropbox/AIME/Panda_C18neg/C18_Non_Exposed_CasesControls/PANDA_input",
           data = "C18_case_control_noexposure_classification_nonorm.RData",
           outloc = "C:/Users/QiYanYan/Dropbox/AIME/Panda_C18neg/C18_Non_Exposed_CasesControls/C18_Annotation/sigfeature_annotation/LipidMaps",
           database = "LipidMaps",
           mode = "neg")
annotation(wd = "C:/Users/QiYanYan/Dropbox/AIME/Panda_C18neg/C18_Non_Exposed_CasesControls/PANDA_input",
           data = "C18_case_control_noexposure_classification_nonorm.RData",
           outloc = "C:/Users/QiYanYan/Dropbox/AIME/Panda_C18neg/C18_Non_Exposed_CasesControls/C18_Annotation/sigfeature_annotation/KEGG",
           database = "KEGG",
           mode = "neg")

annotation(wd = "C:/Users/QiYanYan/Dropbox/AIME/Panda_HILICpos/HILIC_Non_Exposed_CasesControls/PANDA_input",
           data = "HILIC_case_control_noexposure_classification_nonorm.RData",
           outloc = "C:/Users/QiYanYan/Dropbox/AIME/Panda_HILICpos/HILIC_Non_Exposed_CasesControls/HILIC_Annotation/sigfeature_annotation/LipidMaps",
           database = "LipidMaps",
           mode = "pos")
annotation(wd = "C:/Users/QiYanYan/Dropbox/AIME/Panda_HILICpos/HILIC_Controls_ExpoUnexpo/PANDA_input",
           data = "HILIC_control_expo_unexpo_classification_nonorm.RData",
           outloc = "C:/Users/QiYanYan/Dropbox/AIME/Panda_HILICpos/HILIC_Controls_ExpoUnexpo/HILIC_Annotation/sigfeature_annotation/LipidMaps",
           database = "LipidMaps",
           mode = "pos")

############################
# Annotate all features
############################
library(xmsPANDA)
setwd("C:/Users/QiYanYan/Dropbox/AIME/Panda_C18neg/")

# # class <- read.csv(file = "HILIC_classlabels_for_panda_all.txt", sep = '\t', header = T)
# # feature <- read.csv(file = "HILIC_ftrsmzcalib_combat_ordered_all.txt", sep = '\t', header = T)
# class <- "C:/Users/QiYanYan/Dropbox/AIME/Panda_C18neg/C18_classlabels_for_panda_all.txt"
# feature <- "C:/Users/QiYanYan/Dropbox/AIME/Panda_C18neg/C18_ftrsmzcalib_combat_ordered_all.txt"
# outloc <- "C:/Users/QiYanYan/Dropbox/AIME/Panda_C18neg"
# 
# ready_for_regression<-data_preprocess(Xmat=NA,Ymat=NA,feature_table_file=feature,parentoutput_dir=outloc,class_labels_file=class,num_replicates=3,feat.filt.thresh=NA,
#                                       summarize.replicates=TRUE,summary.method="median",all.missing.thresh=0,group.missing.thresh=0,
#                                       log2transform=FALSE,medcenter=FALSE,znormtransform=FALSE,quantile_norm=FALSE,lowess_norm=FALSE,madscaling=FALSE,missing.val=0,
#                                       samplermindex=NA,rep.max.missing.thresh=0.5,summary.na.replacement="halfdatamin",featselmethod=NA)
# feature <- as.data.frame(ready_for_regression$data_matrix_afternorm_scaling)
# na_count <-sapply(feature, function(y) sum(is.na(y)))
# summary(na_count)
# 
# row.names(feature) <- c(paste("met_",1:nrow(feature),sep = ""))
# after.prepro.linkid <- feature[,1:2]
# after.prepro.feature <- t(feature[,-c(1:2)])
# after.prepro.feature <- after.prepro.feature[order(row.names(after.prepro.feature)), ]
# 
# save(after.prepro.feature, after.prepro.linkid, file = "C18_annotation_input.RData")
load("C18_annotation_input.RData")

input_data <- as.data.frame(t(after.prepro.feature))
input_data <- cbind(after.prepro.linkid,input_data)

data(adduct_table)
data(adduct_weights)

###########Parameters to change##############

dataA<-input_data

outloc<-"C:/Users/QiYanYan/Dropbox/AIME/Panda_C18neg/C18_Annotation/sigfeature_annotation/HMDB"

max.mz.diff<-10  #mass search tolerance for DB matching in ppm
max.rt.diff<-10 #retention time tolerance between adducts/isotopes
corthresh<-0.7 #correlation threshold between adducts/isotopes
max_isp=5 #maximum number of isotopes to search for
mass_defect_window=0.01 #mass defect window for isotope search

num_nodes<-4   #number of cores to be used; 2 is recommended for desktop computers due to high memory consumption

db_name="HMDB" #other options: HMDB,Custom,KEGG, LipidMaps, T3DB
status=NA #other options: "Detected", NA, "Detected and Quantified", "Expected and Not Quantified"
num_sets<-300 #number of sets into which the total number of database entries should be split into;

mode<-"neg" #ionization mode
# queryadductlist=c("M+H","M+2H","M+H+NH4","M+ACN+2H","M+2ACN+2H","M+NH4","M+Na","M+ACN+H","M+ACN+Na","M+2ACN+H","2M+H","2M+Na","2M+ACN+H","M+2Na-H","M+H-H2O","M+H-2H2O") #other options: c("M-H","M-H2O-H","M+Na-2H","M+Cl","M+FA-H"); c("positive"); c("negative"); c("all");see data(adduct_table) for complete list
queryadductlist=c("all")

#provide list of database IDs (depending upon selected database) for annotating only specific metabolites
customIDs<-NA #c("HMDB15352","HMDB60549","HMDB00159","HMDB00222"); read.csv("/Users/mzmatch_95stdmx_HMDBIDs.csv")
customDB<-NA 

#########################

dataA<-unique(dataA)
print(dim(dataA))
print(format(Sys.time(), "%a %b %d %X %Y"))

system.time(annotres<-multilevelannotation(dataA=dataA,max.mz.diff=max.mz.diff,max.rt.diff=max.rt.diff,cormethod="pearson",
                                           num_nodes=num_nodes,queryadductlist=queryadductlist,
                                           mode=mode,outloc=outloc,db_name=db_name, adduct_weights=adduct_weights,num_sets=num_sets,allsteps=TRUE,
                                           corthresh=corthresh,NOPS_check=TRUE,customIDs=customIDs,missing.value=NA,
                                           deepsplit=2,networktype="unsigned",minclustsize=10,module.merge.dissimilarity=0.2,filter.by=c("M-H"),
                                           biofluid.location=NA,origin=NA,status=status,boostIDs=NA,max_isp=max_isp,
                                           customDB=customDB,MplusH.abundance.ratio.check = FALSE, min_ions_perchem = 1,
                                           HMDBselect="all",mass_defect_window=mass_defect_window,pathwaycheckmode="pm",mass_defect_mode="both")
)


print(format(Sys.time(), "%a %b %d %X %Y"))

################################
# Verify annotation using in-house library
################################

# mz tolerance = 5ppm, time tolerance = 15s

library(fuzzyjoin)

tolerance = 5

HILIC_library <- read.table(file = "C:/Users/QiYan/Dropbox/AIME/PNS_Ritz/IROA_input_file_hilicpos_vjune282018.txt", sep = "\t", header = T)
C18_library <- read.table(file = "C:/Users/QiYan/Dropbox/AIME/PNS_Ritz/IROA_input_file_c18neg_vjune282018.txt", sep = "\t", header = T)

load("C:/Users/QiYan/Dropbox/AIME/Panda_HILICpos/HILIC_annotation_input.RData")
after.prepro.linkid$delta_mz = (tolerance) * (after.prepro.linkid$mz/1e+06)
after.prepro.linkid$min_mz = after.prepro.linkid$mz - after.prepro.linkid$delta_mz
after.prepro.linkid$max_mz = after.prepro.linkid$mz + after.prepro.linkid$delta_mz

HILIC_verify <- data.frame()
for(i in 1:nrow(after.prepro.linkid)){
  rownum <- which(HILIC_library$mz >= after.prepro.linkid[i,4] & HILIC_library$mz <= after.prepro.linkid[i,5] & abs(HILIC_library$time-after.prepro.linkid[i,2])<=15)
  if (length(rownum)!=0){
    temp <- cbind(after.prepro.linkid[i,],HILIC_library[rownum,])
    HILIC_verify <- rbind(HILIC_verify, temp)
  }
  rm(temp)
}

HILIC_verify <- HILIC_verify[,c(1,2,6:10)]
colnames(HILIC_verify) <- c("mz.x","time.x","mz.y","time.y","ID","KEGGID","HMDBID")

# HILIC_verify <- difference_left_join(after.prepro.linkid, HILIC_library, by = "mz", max_dist = 0.005)
# HILIC_verify <- HILIC_verify[-which(is.na(HILIC_verify$mz.y)),]
# timeout <- which(abs(HILIC_verify$time.x - HILIC_verify$time.y)>=15)
# HILIC_verify <- HILIC_verify[-timeout,]
write.table(HILIC_verify, file = "C:/Users/QiYan/Dropbox/AIME/Panda_HILICpos/HILIC_Annotation/HILIC_annotation_verified.txt", sep = "\t", row.names = F,quote = F)

load("C:/Users/QiYan/Dropbox/AIME/Panda_C18neg/C18_annotation_input.RData")
# C18_verify <- difference_left_join(after.prepro.linkid, C18_library, by = "mz", max_dist = 0.005)
# C18_verify <- C18_verify[-which(is.na(C18_verify$mz.y)),]
# timeout <- which(abs(C18_verify$time.x - C18_verify$time.y)>=15)
# C18_verify <- C18_verify[-timeout,]
after.prepro.linkid$delta_mz = (tolerance) * (after.prepro.linkid$mz/1e+06)
after.prepro.linkid$min_mz = after.prepro.linkid$mz - after.prepro.linkid$delta_mz
after.prepro.linkid$max_mz = after.prepro.linkid$mz + after.prepro.linkid$delta_mz

C18_verify <- data.frame()
for(i in 1:nrow(after.prepro.linkid)){
  rownum <- which(C18_library$mz >= after.prepro.linkid[i,4] & C18_library$mz <= after.prepro.linkid[i,5] & abs(C18_library$time-after.prepro.linkid[i,2])<=15)
  if (length(rownum)!=0){
    temp <- cbind(after.prepro.linkid[i,],C18_library[rownum,])
    C18_verify <- rbind(C18_verify, temp)
  }
  rm(temp)
}

C18_verify <- C18_verify[,c(1,2,6:10)]
colnames(C18_verify) <- c("mz.x","time.x","mz.y","time.y","ID","KEGGID","HMDBID")

write.table(C18_verify, file = "C:/Users/QiYan/Dropbox/AIME/Panda_C18neg/C18_Annotation/C18_annotation_verified.txt", sep = "\t", row.names = F,quote = F)
