library(tidyverse)
library(zoo)
library(xlsx)

# Note: if we use mummichog version 1, then we use MTBNK+PLSDAvip2fc0.58 features, otherwise we use PLSDAvip2fc0 features

final_table <- function(wd_annotation, wd_feature, data_feature_pls, data_feature_mummichogInput, wd_pathway,
                        data_pathway, data_verified, whole_list=FALSE){
  if(whole_list == FALSE){
    ####################
    # annotation
    ####################
    setwd(wd_annotation)
    annotation_HMDB <- read.csv(file = "HMDB/Stage4.csv", sep = ",", header = T)
    annotation_KEGG <- read.csv(file = "KEGG/Stage4.csv", sep = ",", header = T)
    annotation_LipidMaps <- read.csv(file = "LipidMaps/Stage4.csv", sep = ",", header = T)
    
    ####################
    # feature selection - mummichog server
    ####################
    # setwd("C:/Users/QiYan/Dropbox/AIME/Panda_HILICpos/HILIC_Non_Exposed_CasesControls/PANDA_output_PLSDA/")
    # load("Res_PLSDA_result_2018-05-01_vip2fc0.58.RData")
    # feature_PLSDA <- save.plsresults.allfeatures[which(save.plsresults.allfeatures$vip>=2),]
    # feature_PLSDA <- feature_PLSDA[order(-feature_PLSDA$vip),]
    # rm(save.cv.accuracy, save.HMDB, save.KEGG, save.LipidMaps, save.mummichog_PLSDA_VIP2, save.plsresults.allfeatures, save.plsresults.sigfeatures)
    
    ####################
    # pathway analysis - mummichog server
    ####################
    # setwd("C:/Users/QiYan/Dropbox/AIME/Panda_HILICpos/HILIC_Non_Exposed_CasesControls/HILIC_mummichog/HILIC_mummichogOL_vip2fc0/tables/")
    # mummichog_input <- read.table(file = "userInputData.txt", sep = "\t", header = T)
    # mummichog_pathway <- read.xlsx(file = "mcg_pathwayanalysis_.xlsx", 1, header = T)
    # mummichog_empirical <- read.table(file = "ListOfEmpiricalCompounds.tsv", sep = "\t", header = T)
    # mummichog_link <- read.xlsx(file = "mcg_pathwayanalysis_.xlsx", 2, header = T)    # second sheet copied from the web, is the link between feature and name
    # 
    # # rearrange the link table
    # flag <- which(grepl(pattern = "E", x = as.character(mummichog_link[,1])))
    # temp <- mummichog_link[flag,]
    # temp <- temp[,c(1,2,4)]
    # colnames(temp) <- c("EmpiricalCompound", "CompoundName", "KEGGID")
    # mummichog_link$EmpiricalCompound <- na.locf(mummichog_link$EmpiricalCompound)
    # mummichog_link <- mummichog_link[-flag,]
    # mummichog_link <- merge(x = mummichog_link, y = temp, by = "EmpiricalCompound", all.x = T)
    # rm(temp,flag)
    # 
    # # rearrange pathway data from wide to long
    # mummichog_pathway <- mummichog_pathway[which(mummichog_pathway$p.value<=0.05),]  #add overlap size or not
    # mummichog_pathway_expand <- data.frame()
    # for(i in 1:nrow(mummichog_pathway)){
    #   temp <- mummichog_pathway[i,]
    #   temp.compound <- as.character(temp$overlap_EmpiricalCompounds..id.)
    #   temp.compound <- as.data.frame(unlist(strsplit(temp.compound, ",")))
    #   temp.rest <- temp[,c(1:4)]
    #   temp.merge <- merge(x = temp.rest, y = temp.compound, all.y = T)
    #   mummichog_pathway_expand <- rbind(mummichog_pathway_expand, temp.merge)
    # }
    # colnames(mummichog_pathway_expand) <- c("pathway", "overlap_size", "pathway_size", "p.value", "EmpiricalCompound")
    # rm(temp, temp.compound, temp.merge, temp.rest)
    # 
    # # combine pathway with link
    # mummichog_pathway_complete <- merge(mummichog_link, mummichog_pathway_expand, by = "EmpiricalCompound", all = T)
    # mummichog_pathway_complete <- mummichog_pathway_complete[-which(is.na(mummichog_pathway_complete$pathway)),]
    # mummichog_pathway_complete <- mummichog_pathway_complete[order(mummichog_pathway_complete$pathway),]
    # 
    # names(mummichog_pathway_complete)[names(mummichog_pathway_complete) == "Input.m.z"] <- "mz"
    # names(mummichog_pathway_complete)[names(mummichog_pathway_complete) == "Retention.time"] <- "time"
    # 
    # mummichog_pathway_complete$time <- round(mummichog_pathway_complete$time, 1)
    # mummichog_pathway_complete$mz <- round(as.numeric(as.character(mummichog_pathway_complete$mz)), 4)
    
    ####################
    # feature selection - mummichog version 1
    ####################
    setwd(wd_feature)
    load(data_feature_pls)
    feature_MTBNK <- read.csv(file = data_feature_mummichogInput,
                              sep = "\t", header = T)
    feature_MTBNK <- feature_MTBNK[which(feature_MTBNK$p.value <= 0.04), c(1:2)]
    feature_PLSDA <- merge(feature_MTBNK, save.plsresults.allfeatures, by = c("mz","time"), all.X = T)
    feature_PLSDA <- feature_PLSDA[order(-feature_PLSDA$vip),]
    rm(save.cv.accuracy, save.HMDB, save.KEGG, save.LipidMaps, save.mummichog_PLSDA_VIP2, save.plsresults.allfeatures, save.plsresults.sigfeatures, feature_MTBNK)
    
    ####################
    # pathway analysis - mummichog version 1
    ####################
    setwd(wd_pathway)
    mummichog_pathway <- read.xlsx(file = data_pathway, 2, header = T)
    
    mummichog_pathway_complete <- mummichog_pathway
    names(mummichog_pathway_complete)[names(mummichog_pathway_complete) == "m.z"] <- "mz"
    mummichog_pathway_complete$mz <- round(as.numeric(as.character(mummichog_pathway_complete$mz)), 4)
    mummichog_pathway_complete <- mummichog_pathway_complete[,-4]
    
    ####################
    # combine three datasets
    ####################
    feature_PLSDA$time <- round(feature_PLSDA$time, 1)
    annotation_HMDB$time <- round(annotation_HMDB$time, 1)
    annotation_KEGG$time <- round(annotation_KEGG$time, 1)
    annotation_LipidMaps$time <- round(annotation_LipidMaps$time, 1)
    
    feature_PLSDA$mz <- round(feature_PLSDA$mz, 4)
    annotation_HMDB$mz <- round(annotation_HMDB$mz, 4)
    annotation_KEGG$mz <- round(annotation_KEGG$mz, 4)
    annotation_LipidMaps$mz <- round(annotation_LipidMaps$mz, 4)
    
    # merge only confidence score >= 2
    
    # KEGG first
    annotation_KEGG_highconf <- annotation_KEGG[which(annotation_KEGG$Confidence >= 1),]
    feature_w_annotation <- merge(feature_PLSDA[,c(1:4)], annotation_KEGG_highconf[,c(1,2,5,6,10)], by = c("mz", "time"), all = T)
    # feature_w_annotation <- merge(feature_w_annotation[,c(1:4)], annotation_KEGG_highconf[,c(1,2,5,6,10)], by = c("mz", "time"), all = T)
    feature_w_annotation <- feature_w_annotation[-which(is.na(feature_w_annotation$vip)),]
    flag <- which(!is.na(feature_w_annotation$chemical_ID))
    temp.KEGG <- feature_w_annotation[flag,]
    feature_w_annotation <- feature_w_annotation[-flag,]   # the rest pass to HMDB
    
    # HMDB
    annotation_HMDB_highconf <- annotation_HMDB[which(annotation_HMDB$Confidence >= 1),]
    # feature_w_annotation <- merge(feature_PLSDA[,c(1:4)], annotation_HMDB_highconf[,c(1,2,5,6,10)], by = c("mz", "time"), all = T)
    feature_w_annotation <- merge(feature_w_annotation[,c(1:4)], annotation_HMDB_highconf[,c(1,2,5,6,10)], by = c("mz", "time"), all = T)
    feature_w_annotation <- feature_w_annotation[-which(is.na(feature_w_annotation$vip)),]
    flag <- which(!is.na(feature_w_annotation$chemical_ID))
    temp.HMDB <- feature_w_annotation[flag,]
    feature_w_annotation <- feature_w_annotation[-flag,]   
    
    # LipidMaps
    annotation_LipidMaps_highconf <- annotation_LipidMaps[which(annotation_LipidMaps$Confidence >= 1),]
    # feature_w_annotation <- merge(feature_PLSDA[,c(1:4)], annotation_LipidMaps_highconf[,c(1,2,5,6,10)], by = c("mz", "time"), all = T)
    feature_w_annotation <- merge(feature_w_annotation[,c(1:4)], annotation_LipidMaps_highconf[,c(1,2,5,6,10)], by = c("mz", "time"), all = T)
    feature_w_annotation <- feature_w_annotation[-which(is.na(feature_w_annotation$vip)),]
    flag <- which(!is.na(feature_w_annotation$chemical_ID))
    temp.LipidMaps <- feature_w_annotation[flag,]
    feature_w_annotation <- feature_w_annotation[-flag,]
    # add back
    feature_w_annotation <- rbind(feature_w_annotation, temp.HMDB, temp.KEGG, temp.LipidMaps)
    feature_w_annotation <- feature_w_annotation[order(-feature_w_annotation$vip),]
    rm(temp.HMDB, temp.KEGG, temp.LipidMaps)
    
    # add pathway
    
    # merge_all <- merge(feature_w_annotation, mummichog_pathway_complete, by = c("mz", "time"), all = T)  #this is for version 2!!!
    # merge_all <- merge_all[,-c(8,10,16:18)]
    merge_all <- merge(feature_w_annotation, mummichog_pathway_complete, by = c("mz"), all = T)
    # remove those non-significant features picked by mummichog
    merge_all <- merge_all[-which(is.na(merge_all$vip)),]
    merge_all <- merge_all[order(-merge_all$vip),]
    
    if(is.na(data_verified)==TRUE){
      return(merge_all)
    }else{
      verify_all <- read.table(file = data_verified, sep = "\t", header = T)
      verify_all <- verify_all[,c(1,2,6,7)]
      colnames(verify_all) <- c("mz","time","Verified_KEGGID","Verified_HMDBID")
      verify_all$mz <- round(verify_all$mz,4)
      verify_all$time <- round(verify_all$time,1)
      merge_all_verify <- merge(merge_all,verify_all,by = c("mz","time"),all.x = T)
      return(merge_all_verify)
    }
  }else{
    #####################################
    # Generate whole list of metabolites for searching in the future
    #####################################
    
    setwd(wd_annotation)
    annotation_HMDB <- read.csv(file = "HMDB/Stage4.csv", sep = ",", header = T) # round before R!!!!!
    annotation_KEGG <- read.csv(file = "KEGG/Stage4.csv", sep = ",", header = T)
    annotation_LipidMaps <- read.csv(file = "LipidMaps/Stage4.csv", sep = ",", header = T)
    
    setwd(wd_feature)
    load(data_feature_pls)
    feature_PLSDA <- save.plsresults.allfeatures[order(save.plsresults.allfeatures$mz),c(1:4)]
    rm(save.cv.accuracy, save.HMDB, save.KEGG, save.LipidMaps, save.mummichog_PLSDA_VIP2, save.plsresults.allfeatures, save.plsresults.sigfeatures)
    
    # merge only confidence score >= 2
    
    feature_PLSDA$time <- round(feature_PLSDA$time, 1)
    annotation_HMDB$time <- round(annotation_HMDB$time, 1)
    annotation_KEGG$time <- round(annotation_KEGG$time, 1)
    annotation_LipidMaps$time <- round(annotation_LipidMaps$time, 1)
    
    feature_PLSDA$mz <- round(feature_PLSDA$mz, 4)
    annotation_HMDB$mz <- round(annotation_HMDB$mz, 4)
    annotation_KEGG$mz <- round(annotation_KEGG$mz, 4)
    annotation_LipidMaps$mz <- round(annotation_LipidMaps$mz, 4)
    
    # HMDB
    annotation_HMDB_highconf <- annotation_HMDB[which(annotation_HMDB$Confidence >= 1),]
    feature_w_annotation <- merge(feature_PLSDA[,c(1:4)], annotation_HMDB_highconf[,c(1,2,5,6,10,13)], by = c("mz", "time"), all = T)
    # feature_w_annotation <- merge(feature_w_annotation[,c(1:4)], annotation_HMDB_highconf[,c(1,2,5,6,10)], by = c("mz", "time"), all = T)
    # feature_w_annotation <- feature_w_annotation[-which(is.na(feature_w_annotation$vip)),]
    flag <- which(!is.na(feature_w_annotation$chemical_ID))
    temp.HMDB <- feature_w_annotation[flag,]
    feature_w_annotation <- feature_w_annotation[-flag,]
    
    # KEGG first
    annotation_KEGG_highconf <- annotation_KEGG[which(annotation_KEGG$Confidence >= 1),]
    # feature_w_annotation <- merge(feature_PLSDA[,c(1:4)], annotation_KEGG_highconf[,c(1,2,5,6,10)], by = c("mz", "time"), all = T)
    feature_w_annotation <- merge(feature_w_annotation[,c(1:4)], annotation_KEGG_highconf[,c(1,2,5,6,10,13)], by = c("mz", "time"), all = T)
    feature_w_annotation <- feature_w_annotation[-which(is.na(feature_w_annotation$vip)),]
    flag <- which(!is.na(feature_w_annotation$chemical_ID))
    temp.KEGG <- feature_w_annotation[flag,]
    feature_w_annotation <- feature_w_annotation[-flag,]   # the rest pass to HMDB
    
    # LipidMaps
    annotation_LipidMaps_highconf <- annotation_LipidMaps[which(annotation_LipidMaps$Confidence >= 1),]
    # feature_w_annotation <- merge(feature_PLSDA[,c(1:4)], annotation_LipidMaps_highconf[,c(1,2,5,6,10)], by = c("mz", "time"), all = T)
    feature_w_annotation <- merge(feature_w_annotation[,c(1:4)], annotation_LipidMaps_highconf[,c(1,2,5,6,10,13)], by = c("mz", "time"), all = T)
    feature_w_annotation <- feature_w_annotation[-which(is.na(feature_w_annotation$vip)),]
    flag <- which(!is.na(feature_w_annotation$chemical_ID))
    temp.LipidMaps <- feature_w_annotation[flag,]
    feature_w_annotation <- feature_w_annotation[-flag,]
    # add back
    feature_w_annotation <- rbind(feature_w_annotation, temp.HMDB, temp.KEGG, temp.LipidMaps)
    feature_w_annotation <- feature_w_annotation[order(feature_w_annotation$mz),]
    rm(temp.HMDB, temp.KEGG, temp.LipidMaps)
    
    merge_all <- feature_w_annotation[-which(is.na(feature_w_annotation$vip)),]
    
    if(is.na(data_verified)==TRUE){
      return(merge_all)
    }else{
      verify_all <- read.table(file = data_verified, sep = "\t", header = T)
      verify_all <- verify_all[,c(1,2,6,7)]
      colnames(verify_all) <- c("mz","time","Verified_KEGGID","Verified_HMDBID")
      verify_all$mz <- round(verify_all$mz,4)
      verify_all$time <- round(verify_all$time,1)
      merge_all_verify <- merge(merge_all,verify_all,by = c("mz","time"),all.x = T)
      return(merge_all_verify)
    }
  }
  
  
}

## Now the problem is that we need to use mummichog 1 instead of mummichog 2!!

HILIC_autism <- final_table(wd_annotation = "C:/Users/QiYan/Dropbox/AIME/Panda_HILICpos/HILIC_Annotation/sigfeature_annotation/",
                            wd_feature = "C:/Users/QiYan/Dropbox/AIME/Panda_HILICpos/HILIC_Non_Exposed_CasesControls/PANDA_output_PLSDA/",
                            data_feature_pls = "Res_PLSDA_result_2018-05-01_vip2fc0.58.RData",
                            data_feature_mummichogInput = "C:/Users/QiYan/Dropbox/AIME/Panda_HILICpos/HILIC_Non_Exposed_CasesControls/HILIC_mummichog/mummichogMTBNK_input_vip2fc0.58.txt",
                            wd_pathway = "C:/Users/QiYan/Dropbox/AIME/Panda_HILICpos/HILIC_Non_Exposed_CasesControls/HILIC_mummichog/HILIC_mummichogMTBNK_vip2fc0.58/tsv/",
                            data_pathway = "mcg_pathwayanalysis_HILIC_mummichogMTBNK_vip2fc0.58.xlsx",
                            data_verified = "C:/Users/QiYan/Dropbox/AIME/Panda_HILICpos/HILIC_Annotation/HILIC_annotation_verified.txt")

C18_autism <- final_table(wd_annotation = "C:/Users/QiYan/Dropbox/AIME/Panda_C18neg/C18_Annotation/sigfeature_annotation/",
                          wd_feature = "C:/Users/QiYan/Dropbox/AIME/Panda_C18neg/C18_Non_Exposed_CasesControls/PANDA_output_PLSDA/",
                          data_feature_pls = "Res_PLSDA_result_2018-05-01_vip2fc0.58.RData",
                          data_feature_mummichogInput = "C:/Users/QiYan/Dropbox/AIME/Panda_C18neg/C18_Non_Exposed_CasesControls/C18_mummichog/mummichogMTBNK_input_vip2fc0.58.txt",
                          wd_pathway = "C:/Users/QiYan/Dropbox/AIME/Panda_C18neg/C18_Non_Exposed_CasesControls/C18_mummichog/C18_mummichogMTBNK_vip2fc0.58/tsv/",
                          data_pathway = "mcg_pathwayanalysis_C18_mummichogMTBNK_vip2fc0.58.xlsx",
                          data_verified = "C:/Users/QiYan/Dropbox/AIME/Panda_C18neg/C18_Annotation/C18_annotation_verified.txt")

HILIC_autism <- HILIC_autism[-which(is.na(HILIC_autism$pathway)&is.na(HILIC_autism$Verified_KEGGID)&is.na(HILIC_autism$Verified_HMDBID)),]
C18_autism <- C18_autism[-which(is.na(C18_autism$pathway)&is.na(C18_autism$Verified_KEGGID)&is.na(C18_autism$Verified_HMDBID)),]

write.table(HILIC_autism, file = "C:/Users/QiYan/Dropbox/AIME/Results/Autism/HILIC_final_table.txt", sep = "\t", row.names = F,quote = F)
write.table(C18_autism, file = "C:/Users/QiYan/Dropbox/AIME/Results/Autism/C18_final_table.txt", sep = "\t", row.names = F,quote = F)

# Air pollution
HILIC_ap <- final_table(wd_annotation = "C:/Users/QiYan/Dropbox/AIME/Panda_HILICpos/HILIC_Annotation/sigfeature_annotation/",
                        wd_feature = "C:/Users/QiYan/Dropbox/AIME/Panda_HILICpos/HILIC_Controls_ExpoUnexpo/PANDA_output_PLSDA/",
                        data_feature_pls = "Res_PLSDA_result_2018-08-19_vip2fc0.RData",
                        data_feature_mummichogInput = "C:/Users/QiYan/Dropbox/AIME/Panda_HILICpos/HILIC_Controls_ExpoUnexpo/HILIC_mummichog/mummichog_input_vip2fc0_2018-08-19.txt",
                        wd_pathway = "C:/Users/QiYan/Dropbox/AIME/Panda_HILICpos/HILIC_Controls_ExpoUnexpo/HILIC_mummichog/HILIC_mummichog_vip2fc0_0820_major/tsv/",
                        data_pathway = "mcg_pathwayanalysis_HILIC_mummichog_vip2fc0_0820.xlsx",
                        data_verified = "C:/Users/QiYan/Dropbox/AIME/Panda_HILICpos/HILIC_Annotation/HILIC_annotation_verified.txt"
)

C18_ap <- final_table(wd_annotation = "C:/Users/QiYan/Dropbox/AIME/Panda_C18neg/C18_Annotation/sigfeature_annotation/",
                      wd_feature = "C:/Users/QiYan/Dropbox/AIME/Panda_C18neg/C18_Controls_ExpoUnexpo/PANDA_output_PLSDA/",
                      data_feature_pls = "Res_PLSDA_result_2018-08-19_vip2fc0.RData",
                      data_feature_mummichogInput = "C:/Users/QiYan/Dropbox/AIME/Panda_C18neg/C18_Controls_ExpoUnexpo/C18_mummichog/mummichog_input_vip2fc0_2018-08-19.txt",
                      wd_pathway = "C:/Users/QiYan/Dropbox/AIME/Panda_C18neg/C18_Controls_ExpoUnexpo/C18_mummichog/C18_mummichog_vip2fc0_0820_major/tsv/",
                      data_pathway = "mcg_pathwayanalysis_C18_mummichog_vip2fc0_0820.xlsx",
                      data_verified = "C:/Users/QiYan/Dropbox/AIME/Panda_C18neg/C18_Annotation/C18_annotation_verified.txt")

HILIC_ap <- HILIC_ap[-which(is.na(HILIC_ap$chemical_ID)&is.na(HILIC_ap$pathway)&is.na(HILIC_ap$Verified_KEGGID)&is.na(HILIC_ap$Verified_HMDBID)),]
C18_ap <- C18_ap[-which(is.na(C18_ap$chemical_ID)&is.na(C18_ap$pathway)&is.na(C18_ap$Verified_KEGGID)&is.na(C18_ap$Verified_HMDBID)),]

write.table(HILIC_ap, file = "C:/Users/QiYan/Dropbox/AIME/Results/Air pollution/HILIC_final_table_0920.txt", sep = "\t", row.names = F,quote = F)
write.table(C18_ap, file = "C:/Users/QiYan/Dropbox/AIME/Results/Air pollution/C18_final_table_0920.txt", sep = "\t", row.names = F,quote = F)


# Air pollution Whole list
HILIC_ap <- final_table(wd_annotation = "C:/Users/QiYan/Dropbox/AIME/Panda_HILICpos/HILIC_Annotation/sigfeature_annotation/",
                        wd_feature = "C:/Users/QiYan/Dropbox/AIME/Panda_HILICpos/HILIC_Controls_ExpoUnexpo/PANDA_output_PLSDA/",
                        data_feature_pls = "Res_PLSDA_result_2018-08-19_vip2fc0.RData",
                        data_feature_mummichogInput = "C:/Users/QiYan/Dropbox/AIME/Panda_HILICpos/HILIC_Controls_ExpoUnexpo/HILIC_mummichog/mummichog_input_vip2fc0_2018-08-19.txt",
                        wd_pathway = "C:/Users/QiYan/Dropbox/AIME/Panda_HILICpos/HILIC_Controls_ExpoUnexpo/HILIC_mummichog/HILIC_mummichog_vip2fc0_0820_major/tsv/",
                        data_pathway = "mcg_pathwayanalysis_HILIC_mummichog_vip2fc0_0820.xlsx",
                        data_verified = "C:/Users/QiYan/Dropbox/AIME/Panda_HILICpos/HILIC_Annotation/HILIC_annotation_verified.txt",
                        whole_list = TRUE)

C18_ap <- final_table(wd_annotation = "C:/Users/QiYan/Dropbox/AIME/Panda_C18neg/C18_Annotation/sigfeature_annotation/",
                      wd_feature = "C:/Users/QiYan/Dropbox/AIME/Panda_C18neg/C18_Controls_ExpoUnexpo/PANDA_output_PLSDA/",
                      data_feature_pls = "Res_PLSDA_result_2018-08-19_vip2fc0.RData",
                      data_feature_mummichogInput = "C:/Users/QiYan/Dropbox/AIME/Panda_C18neg/C18_Controls_ExpoUnexpo/C18_mummichog/mummichog_input_vip2fc0_2018-08-19.txt",
                      wd_pathway = "C:/Users/QiYan/Dropbox/AIME/Panda_C18neg/C18_Controls_ExpoUnexpo/C18_mummichog/C18_mummichog_vip2fc0_0820_major/tsv/",
                      data_pathway = "mcg_pathwayanalysis_C18_mummichog_vip2fc0_0820.xlsx",
                      data_verified = "C:/Users/QiYan/Dropbox/AIME/Panda_C18neg/C18_Annotation/C18_annotation_verified.txt",
                      whole_list = TRUE)

HILIC_ap <- HILIC_ap[-which(is.na(HILIC_ap$chemical_ID)&is.na(HILIC_ap$Verified_KEGGID)&is.na(HILIC_ap$Verified_HMDBID)),]
C18_ap <- C18_ap[-which(is.na(C18_ap$chemical_ID)&is.na(C18_ap$Verified_KEGGID)&is.na(C18_ap$Verified_HMDBID)),]

write.table(HILIC_ap, file = "C:/Users/QiYan/Dropbox/AIME/Results/Air pollution/HILIC_whole_table_0920.txt", sep = "\t", row.names = F,quote = F)
write.table(C18_ap, file = "C:/Users/QiYan/Dropbox/AIME/Results/Air pollution/C18_whole_table_0920.txt", sep = "\t", row.names = F,quote = F)

# Autism whole list
HILIC_autism <- final_table(wd_annotation = "C:/Users/QiYan/Dropbox/AIME/Panda_HILICpos/HILIC_Annotation/sigfeature_annotation/",
                            wd_feature = "C:/Users/QiYan/Dropbox/AIME/Panda_HILICpos/HILIC_Non_Exposed_CasesControls/PANDA_output_PLSDA/",
                            data_feature_pls = "Res_PLSDA_result_2018-05-01_vip2fc0.58.RData",
                            data_feature_mummichogInput = "C:/Users/QiYan/Dropbox/AIME/Panda_HILICpos/HILIC_Non_Exposed_CasesControls/HILIC_mummichog/mummichogMTBNK_input_vip2fc0.58.txt",
                            wd_pathway = "C:/Users/QiYan/Dropbox/AIME/Panda_HILICpos/HILIC_Non_Exposed_CasesControls/HILIC_mummichog/HILIC_mummichogMTBNK_vip2fc0.58/tsv/",
                            data_pathway = "mcg_pathwayanalysis_HILIC_mummichogMTBNK_vip2fc0.58.xlsx",
                            data_verified = "C:/Users/QiYan/Dropbox/AIME/Panda_HILICpos/HILIC_Annotation/HILIC_annotation_verified.txt",
                            whole_list = T)

C18_autism <- final_table(wd_annotation = "C:/Users/QiYan/Dropbox/AIME/Panda_C18neg/C18_Annotation/sigfeature_annotation/",
                          wd_feature = "C:/Users/QiYan/Dropbox/AIME/Panda_C18neg/C18_Non_Exposed_CasesControls/PANDA_output_PLSDA/",
                          data_feature_pls = "Res_PLSDA_result_2018-05-01_vip2fc0.58.RData",
                          data_feature_mummichogInput = "C:/Users/QiYan/Dropbox/AIME/Panda_C18neg/C18_Non_Exposed_CasesControls/C18_mummichog/mummichogMTBNK_input_vip2fc0.58.txt",
                          wd_pathway = "C:/Users/QiYan/Dropbox/AIME/Panda_C18neg/C18_Non_Exposed_CasesControls/C18_mummichog/C18_mummichogMTBNK_vip2fc0.58/tsv/",
                          data_pathway = "mcg_pathwayanalysis_C18_mummichogMTBNK_vip2fc0.58.xlsx",
                          data_verified = "C:/Users/QiYan/Dropbox/AIME/Panda_C18neg/C18_Annotation/C18_annotation_verified.txt",
                          whole_list = T)

HILIC_autism <- HILIC_autism[-which(is.na(HILIC_autism$chemical_ID)),]
C18_autism <- C18_autism[-which(is.na(C18_autism$chemical_ID)),]

write.table(HILIC_autism, file = "C:/Users/QiYan/Dropbox/AIME/Results/Autism/HILIC_whole_table.txt", sep = "\t", row.names = F,quote = F)
write.table(C18_autism, file = "C:/Users/QiYan/Dropbox/AIME/Results/Autism/C18_whole_table.txt", sep = "\t", row.names = F,quote = F)
