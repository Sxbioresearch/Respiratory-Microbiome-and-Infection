##### Co-occurrence network analysis

# 1、correlation calculation
library(dplyr)
library(data.table)
meta_new <- fread("DNA_OP_group_27.csv")
meta_new <- cbind(meta_new$sample,meta_new$group1) %>% as.data.frame()
names(meta_new) <- c("#NAME","Disease")
meta_new$Disease[meta_new$Disease == "1"] <- "Time1"
meta_new$Disease[meta_new$Disease == "2"] <- "Time2"
meta_new$Disease[meta_new$Disease == "3"] <- "Time3"
meta_new$Disease[meta_new$Disease == "4"] <- "Time4"
unique(meta_new$Disease)
combined.dat_new <- read.table("DNA_genus_OP_relative_filter_27.csv",header = T,sep = ",",row.names = 1,stringsAsFactors = FALSE)
combined.dat_new <- as.data.frame(t(combined.dat_new))
library(tibble)
combined.dat_new <- rownames_to_column(combined.dat_new,var="#NAME")
combined.dat_new <- merge(combined.dat_new,meta_new,by="#NAME")
names(combined.dat_new)[1] <- "Row.names"
names(combined.dat_new)[85] <- "Group"
# abundance tables were organized separately for each disease state
data <- read.table("DNA_genus_OP_relative_filter_27.csv",header = T,sep = ",",row.names = 1,stringsAsFactors = FALSE)
group <- read.table("DNA_OP_group_27.csv",header = T,sep = ",",stringsAsFactors = FALSE)
T1 <- group[group$group1 == "1",]
T1 <- T1$sample
data_T1 <- data[,T1]
write.csv(data_T1,"abund_T1.csv",quote = F)
T2 <- group[group$group1 == "2",]
T2 <- T2$sample
data_T2 <- data[,T2]
write.csv(data_T2,"abund_T2.csv",quote = F)
T3 <- group[group$group1 == "3",]
T3 <- T3$sample
data_T3 <- data[,T3]
write.csv(data_T3,"abund_T3.csv",quote = F)
T4 <- group[group$group1 == "4",]
T4 <- T4$sample
data_T4 <- data[,T4]
write.csv(data_T4,"abund_T4.csv",quote = F)
# net edges were organized into correlation matrix
##Time1
dat.tmp_new <- combined.dat_new %>%    
  filter(Group == "Time1") %>%       
  tibble::column_to_rownames("Row.names") %>%
  select(-Group)
Corr_new <- rcorr(as.matrix(dat.tmp_new) , type="spearman")
#occor.r_new <- Corr_new$r
# occor.p_new <- Corr_new$P
Corr_new$P[is.na(Corr_new$P)] <- 0
Corr_new$P[Corr_new$P >= 0.05] <- -1
Corr_new$P[Corr_new$P < 0.05 & Corr_new$P >= 0] <- 1
Corr_new$P[Corr_new$P == -1] <- 0
Corr_new$r[abs(Corr_new$r) <= 0.6] <- 0
Corr_new$r <- Corr_new$r * Corr_new$P
occor.r_new <- Corr_new$r
write.csv(occor.r_new,"net_Time1.csv",quote = F)
##Time2
dat.tmp_new <- combined.dat_new %>%    
  filter(Group == "Time2") %>%       
  tibble::column_to_rownames("Row.names") %>%
  select(-Group)
##计算相关矩阵
Corr_new <- rcorr(as.matrix(dat.tmp_new) , type="spearman") 
#occor.r_new <- Corr_new$r
# occor.p_new <- Corr_new$P
Corr_new$P[is.na(Corr_new$P)] <- 0
Corr_new$P[Corr_new$P >= 0.05] <- -1
Corr_new$P[Corr_new$P < 0.05 & Corr_new$P >= 0] <- 1
Corr_new$P[Corr_new$P == -1] <- 0
Corr_new$r[abs(Corr_new$r) <= 0.6] <- 0
Corr_new$r <- Corr_new$r * Corr_new$P
occor.r_new <- Corr_new$r
write.csv(occor.r_new,"net_Time2.csv",quote = F)
##Time3
dat.tmp_new <- combined.dat_new %>%    
  filter(Group == "Time3") %>%      
  tibble::column_to_rownames("Row.names") %>%
  select(-Group)
Corr_new <- rcorr(as.matrix(dat.tmp_new) , type="spearman")   
#occor.r_new <- Corr_new$r
# occor.p_new <- Corr_new$P
Corr_new$P[is.na(Corr_new$P)] <- 0
Corr_new$P[Corr_new$P >= 0.05] <- -1
Corr_new$P[Corr_new$P < 0.05 & Corr_new$P >= 0] <- 1
Corr_new$P[Corr_new$P == -1] <- 0
Corr_new$r[abs(Corr_new$r) <= 0.6] <- 0
Corr_new$r <- Corr_new$r * Corr_new$P
occor.r_new <- Corr_new$r
write.csv(occor.r_new,"net_Time3.csv",quote = F)
##Time4
dat.tmp_new <- combined.dat_new %>%    
  filter(Group == "Time4") %>%      
  tibble::column_to_rownames("Row.names") %>%
  select(-Group)
##计算相关矩阵
Corr_new <- rcorr(as.matrix(dat.tmp_new) , type="spearman")  
#occor.r_new <- Corr_new$r
# occor.p_new <- Corr_new$P
Corr_new$P[is.na(Corr_new$P)] <- 0
Corr_new$P[Corr_new$P >= 0.05] <- -1
Corr_new$P[Corr_new$P < 0.05 & Corr_new$P >= 0] <- 1
Corr_new$P[Corr_new$P == -1] <- 0
Corr_new$r[abs(Corr_new$r) <= 0.6] <- 0
Corr_new$r <- Corr_new$r * Corr_new$P
occor.r_new <- Corr_new$r
write.csv(occor.r_new,"net_Time4.csv",quote = F)

#2、NetMoss
###Time1 - Time2
abund_T1 <- data_T1
abund_T1 <- rownames_to_column(abund_T1,var = "x")
names(abund_T1)[1] <- ""
abund_T2 <- data_T2
abund_T2 <- rownames_to_column(abund_T2,var = "x")
names(abund_T2)[1] <- ""
abund_T3 <- data_T3
abund_T3 <- rownames_to_column(abund_T3,var = "x")
names(abund_T3)[1] <- ""
abund_T4 <- data_T4
abund_T4 <- rownames_to_column(abund_T4,var = "x")
names(abund_T4)[1] <- ""

net_T1 <- read.table("net_Time1.csv",header = T,sep = ",",row.names = 1,stringsAsFactors = FALSE)
net_T2 <- read.table("net_Time2.csv",header = T,sep = ",",row.names = 1,stringsAsFactors = FALSE)
net_T3 <- read.table("net_Time3.csv",header = T,sep = ",",row.names = 1,stringsAsFactors = FALSE)
net_T4 <- read.table("net_Time4.csv",header = T,sep = ",",row.names = 1,stringsAsFactors = FALSE)

nodes_result_T1.vs.T2 = 
  NetMoss(case_dir = abund_T2,
          control_dir = abund_T1,
          net_case_dir = net_T2,
          net_control_dir = net_T1)

nodes_result_T2.vs.T3 = 
  NetMoss(case_dir = abund_T3,
          control_dir = abund_T2,
          net_case_dir = net_T3,
          net_control_dir = net_T2)

nodes_result_T2.vs.T3_new = 
  NetMoss(case_dir = abund_T2,
          control_dir = abund_T3,
          net_case_dir = net_T2,
          net_control_dir = net_T3)

nodes_result_T3.vs.T4 = 
  NetMoss(case_dir = abund_T4,
          control_dir = abund_T3,
          net_case_dir = net_T4,
          net_control_dir = net_T3)

nodes_result_T1.vs.T3 = 
  NetMoss(case_dir = abund_T3,
          control_dir = abund_T1,
          net_case_dir = net_T3,
          net_control_dir = net_T1)

nodes_result_T1.vs.T4 = 
  NetMoss(case_dir = abund_T4,
          control_dir = abund_T1,
          net_case_dir = net_T4,
          net_control_dir = net_T1)

microbeAbb <- read.table("microbeAbb.csv",header = T,sep = ",",stringsAsFactors = FALSE)
NMSS_1.2 = nodes_result_T1.vs.T2[[1]] %>% 
  mutate(`1-2` = NetMoss_Score)
write.csv(NMSS_1.2,"NMSS_1.2.csv",quote = F,row.names = F)
NMSS_2.3 = nodes_result_T2.vs.T3[[1]]  %>% 
  mutate(`2-3` = NetMoss_Score)
write.csv(NMSS_2.3,"NMSS_2.3.csv",quote = F,row.names = F)
NMSS_3.4 = nodes_result_T3.vs.T4[[1]]  %>% 
  mutate(`3-4` = NetMoss_Score)
write.csv(NMSS_3.4,"NMSS_3.4.csv",quote = F,row.names = F)
NMSS_1.3 = nodes_result_T1.vs.T3[[1]]  %>% 
  mutate(`1-3` = NetMoss_Score)
write.csv(NMSS_1.3,"NMSS_1.3.csv",quote = F,row.names = F)
NMSS_1.4 = nodes_result_T1.vs.T4[[1]]  %>% 
  mutate(`1-4` = NetMoss_Score)
write.csv(NMSS_1.4,"NMSS_1.4.csv",quote = F,row.names = F)
res <- merge(NMSS_1.2, 
             NMSS_2.3, 
             by=c("taxon_names")) %>%
  select(taxon_names, `1-2`, `2-3`)
res <- merge(res, 
             NMSS_3.4, 
             by=c("taxon_names")) %>%
  select(taxon_names, `1-2`, `2-3`,`3-4`)
res <- merge(res, 
             NMSS_1.3, 
             by=c("taxon_names")) %>%
  select(taxon_names, `1-2`, `2-3`,`3-4`,`1-3`)
res <- merge(res, 
             NMSS_1.4, 
             by=c("taxon_names")) %>%
  select(taxon_names, `1-2`, `2-3`,`3-4`,`1-3`,`1-4`)
colnames(res)[1] <- "microbeAbb"

#3、calculate noise NetMoss score by randomly shuffle sample IDs 100 times
library(Hmisc)
library(NetMoss2)
abund_T1 <- read.table("abund_T1.csv",header = T,sep = ",",stringsAsFactors = FALSE)
abund_T2 <- read.table("abund_T2.csv",header = T,sep = ",",stringsAsFactors = FALSE)
abund_T3 <- read.table("abund_T3.csv",header = T,sep = ",",stringsAsFactors = FALSE)
abund_T4 <- read.table("abund_T4.csv",header = T,sep = ",",stringsAsFactors = FALSE)
abund_all_new <- cbind.data.frame(abund_T1, abund_T2 %>% select(-X) , abund_T3 %>% select(-X),abund_T4 %>% select(-X))
# calculate noise NetMoss score by randomly shuffle sample IDs 100 times
AllShuff_Res_new <- NULL
for(shuf_new in 1:100){
  shuffled_new  <- abund_all_new %>% tibble::column_to_rownames("X")
  set.seed(shuf_new)
  num_columns <- ncol(shuffled_new)
  shuffled_new <- shuffled_new[,sample(num_columns)] ##Randomly shuffle the order of all columns
  
  abund_T1 <- shuffled_new[,1:27]
  abund_T2 <- shuffled_new[,28:54]
  abund_T3 <- shuffled_new[,55:81]
  abund_T4 <- shuffled_new[,82:108]
  # remove all 0 features
  abund_T1 <- abund_T1[!rownames(abund_T1) %in% names(which(apply(abund_T1, 1, function(x) sum(x != 0)) <=1)), ]
  abund_T2 <- abund_T2[!rownames(abund_T2) %in% names(which(apply(abund_T2, 1,function(x) sum(x != 0)) <=1)),]
  abund_T3 <- abund_T3[!rownames(abund_T3) %in% names(which(apply(abund_T3, 1, function(x) sum(x != 0)) <=1)),]
  abund_T4 <- abund_T4[!rownames(abund_T4) %in% names(which(apply(abund_T4, 1, function(x) sum(x != 0)) <=1)),]
  # T1
  Corr <- rcorr(abund_T1 %>% t(), type="spearman")
  net_T1 <- Corr$r
  p_T1 <- Corr$P
  net_T1_revised <- net_T1
  net_T1_revised[which(p_T1 > 0.05)] = 0
  net_T1_revised[which(abs(net_T1_revised) < 0.6)] = 0
  # T2
  Corr <- rcorr(abund_T2 %>% t(), type="spearman")
  net_T2 <- Corr$r
  p_T2 <- Corr$P
  net_T2_revised <- net_T2
  net_T2_revised[which(p_T2 > 0.05)] = 0
  net_T2_revised[which(abs(net_T2_revised) < 0.6)] = 0
  # T3
  Corr <- rcorr(abund_T3 %>% t(), type="spearman")
  net_T3 <- Corr$r
  p_T3 <- Corr$P
  net_T3_revised <- net_T3
  net_T3_revised[which(p_T3 > 0.05)] = 0
  net_T3_revised[which(abs(net_T3_revised) < 0.6)] = 0
  # T4
  Corr <- rcorr(abund_T4 %>% t(), type="spearman")
  net_T4 <- Corr$r
  p_T4 <- Corr$P
  net_T4_revised <- net_T4
  net_T4_revised[which(p_T4 > 0.05)] = 0
  net_T4_revised[which(abs(net_T4_revised) < 0.6)] = 0
  
  abund_T1 <- rownames_to_column(abund_T1,var = "x")
  names(abund_T1)[1] <- ""
  abund_T2 <- rownames_to_column(abund_T2,var = "x")
  names(abund_T2)[1] <- ""
  abund_T3 <- rownames_to_column(abund_T3,var = "x")
  names(abund_T3)[1] <- ""
  abund_T4 <- rownames_to_column(abund_T4,var = "x")
  names(abund_T4)[1] <- ""
  # netMoss ----
  # T1 vs T2
  netMossRes_T1.vs.T2 = try(NetMoss(case_dir = abund_T2,
                                    control_dir = abund_T1 ,
                                    net_case_dir = net_T2,
                                    net_control_dir = net_T1))
  
  netMossRes_T1.vs.T2_revisedNet = try(NetMoss(case_dir = abund_T2,
                                               control_dir = abund_T1 ,
                                               net_case_dir = net_T2_revised,
                                               net_control_dir = net_T1_revised))
  
  # T2 vs T3
  netMossRes_T2.vs.T3 = try(NetMoss(case_dir = abund_T3,
                                    control_dir = abund_T2 ,
                                    net_case_dir = net_T3,
                                    net_control_dir = net_T2))
  
  netMossRes_T2.vs.T3_revisedNet = try(NetMoss(case_dir = abund_T3,
                                               control_dir = abund_T2 ,
                                               net_case_dir = net_T3_revised,
                                               net_control_dir = net_T2_revised))
  
  # T3 vs T4
  netMossRes_T3.vs.T4 = try(NetMoss(case_dir = abund_T4,
                                    control_dir = abund_T3,
                                    net_case_dir = net_T4,
                                    net_control_dir = net_T3))
  
  netMossRes_T3.vs.T4_revisedNet = try(NetMoss(case_dir = abund_T4,
                                               control_dir = abund_T3 ,
                                               net_case_dir = net_T4_revised,
                                               net_control_dir = net_T3_revised))
  
  # T1 vs T3
  netMossRes_T1.vs.T3 = try(NetMoss(case_dir = abund_T3,
                                    control_dir = abund_T1,
                                    net_case_dir = net_T3,
                                    net_control_dir = net_T1))
  
  netMossRes_T1.vs.T3_revisedNet = try(NetMoss(case_dir = abund_T3,
                                               control_dir = abund_T1 ,
                                               net_case_dir = net_T3_revised,
                                               net_control_dir = net_T1_revised))
  
  # T1 vs T4
  netMossRes_T1.vs.T4 = try(NetMoss(case_dir = abund_T4,
                                    control_dir = abund_T1,
                                    net_case_dir = net_T4,
                                    net_control_dir = net_T1))
  
  netMossRes_T1.vs.T4_revisedNet = try(NetMoss(case_dir = abund_T4,
                                               control_dir = abund_T1 ,
                                               net_case_dir = net_T4_revised,
                                               net_control_dir = net_T1_revised))
  
  # summarise results
  ##T1 VS T2
  if('try-error' %in% class(netMossRes_T1.vs.T2) ){
    NNMS.1_2_new <- cbind.data.frame(taxon_names = abund_all_new$X,
                                     NMSS.1_2_new = NA, 
                                     stringsAsFactors=F)
  }else{
    NNMS.1_2_new <- netMossRes_T1.vs.T2[[1]] %>% mutate(NMSS.1_2_new = NetMoss_Score) %>% select(taxon_names, NMSS.1_2_new)
  }
  
  if('try-error' %in% class(netMossRes_T1.vs.T2_revisedNet) ){
    NNMS.1_2.revisedNet_new <- cbind.data.frame(taxon_names = abund_all_new$X,
                                                NMSS.1_2.revised_new = NA, 
                                                stringsAsFactors=F)
  }else{
    NNMS.1_2.revisedNet_new <- netMossRes_T1.vs.T2_revisedNet[[1]] %>% 
      mutate(NMSS.1_2.revised_new = NetMoss_Score) %>% 
      select(taxon_names, NMSS.1_2.revised_new)
  }
  
  ##T2 VS T3
  if('try-error' %in% class(netMossRes_T2.vs.T3) ){
    NNMS.2_3_new <- cbind.data.frame(taxon_names = abund_all_new$X,
                                     NMSS.2_3_new = NA, 
                                     stringsAsFactors=F)
  }else{
    NNMS.2_3_new <-  netMossRes_T2.vs.T3[[1]] %>% mutate(NMSS.2_3_new = NetMoss_Score) %>% select(taxon_names, NMSS.2_3_new)
  }
  
  if('try-error' %in% class(netMossRes_T2.vs.T3_revisedNet) ){
    NNMS.2_3.revisedNet_new <- cbind.data.frame(taxon_names = abund_all_new$X,
                                                NMSS.2_3.revised_new = NA, 
                                                stringsAsFactors=F)
  }else{
    NNMS.2_3.revisedNet_new <-  netMossRes_T2.vs.T3_revisedNet[[1]] %>% mutate(NMSS.2_3.revised_new = NetMoss_Score) %>% select(taxon_names, NMSS.2_3.revised_new)
    
  }
  
  #T3 VS T4
  if('try-error' %in% class(netMossRes_T3.vs.T4) ){
    NNMS.3_4_new <- cbind.data.frame(taxon_names = abund_all_new$X,
                                     NMSS.3_4_new = NA, 
                                     stringsAsFactors=F)
  }else{
    NNMS.3_4_new <- netMossRes_T3.vs.T4[[1]] %>% mutate(NMSS.3_4_new = NetMoss_Score) %>% select(taxon_names, NMSS.3_4_new)
    
  }
  
  if('try-error' %in% class(netMossRes_T3.vs.T4_revisedNet) ){
    NNMS.3_4.revisedNet_new <- cbind.data.frame(taxon_names = abund_all_new$X,
                                                NMSS.3_4.revised_new = NA, 
                                                stringsAsFactors=F)
  }else{
    NNMS.3_4.revisedNet_new <-  netMossRes_T3.vs.T4_revisedNet[[1]] %>% mutate(NMSS.3_4.revised_new = NetMoss_Score) %>% select(taxon_names, NMSS.3_4.revised_new)
    
  }
  
  ##T1 VS T3
  
  
  if ('try-error' %in% class(netMossRes_T1.vs.T3) || 
      is.null(netMossRes_T1.vs.T3[[1]]) || 
      !is.data.frame(netMossRes_T1.vs.T3[[1]])) {
    
    NNMS.1_3_new <- cbind.data.frame(taxon_names = abund_all_new$X,
                                     NMSS.1_3_new = NA, 
                                     stringsAsFactors = FALSE)
  } else {
    NNMS.1_3_new <- netMossRes_T1.vs.T3[[1]] %>% 
      mutate(NMSS.1_3_new = NetMoss_Score) %>% 
      select(taxon_names, NMSS.1_3_new)
  }
  
  
  
  if('try-error' %in% class(netMossRes_T1.vs.T3_revisedNet) ){
    NNMS.1_3.revisedNet_new <- cbind.data.frame(taxon_names = abund_all_new$X,
                                                NMSS.1_3.revised_new = NA, 
                                                stringsAsFactors=F)
  }else{
    NNMS.1_3.revisedNet_new <- netMossRes_T1.vs.T3_revisedNet[[1]] %>% 
      mutate(NMSS.1_3.revised_new = NetMoss_Score) %>% 
      select(taxon_names, NMSS.1_3.revised_new)
  }
  
  
  ##T1 VS T4
  if('try-error' %in% class(netMossRes_T1.vs.T4) ){
    NNMS.1_4_new <- cbind.data.frame(taxon_names = abund_all_new$X,
                                     NMSS.1_4_new = NA, 
                                     stringsAsFactors=F)
  }else{
    NNMS.1_4_new <- netMossRes_T1.vs.T4[[1]] %>% mutate(NMSS.1_4_new = NetMoss_Score) %>% select(taxon_names, NMSS.1_4_new)
  }
  
  if('try-error' %in% class(netMossRes_T1.vs.T4_revisedNet) ){
    NNMS.1_4.revisedNet_new <- cbind.data.frame(taxon_names = abund_all_new$X,
                                                NMSS.1_4.revised_new = NA, 
                                                stringsAsFactors=F)
  }else{
    NNMS.1_4.revisedNet_new <- netMossRes_T1.vs.T4_revisedNet[[1]] %>% 
      mutate(NMSS.1_4.revised_new = NetMoss_Score) %>% 
      select(taxon_names, NMSS.1_4.revised_new)
  }
  
  shuffRes_new <- merge(merge(merge(merge(merge(NNMS.1_2_new, NNMS.1_2.revisedNet_new, by = "taxon_names", all=T),
                                          merge(NNMS.2_3_new, NNMS.2_3.revisedNet_new, by="taxon_names", all=T),
                                          by="taxon_names", all=T),
                                    merge(NNMS.3_4_new, NNMS.3_4.revisedNet_new, by="taxon_names", all=T),
                                    by="taxon_names", all=T),
                              merge(NNMS.1_3_new, NNMS.1_3.revisedNet_new, by="taxon_names", all=T),
                              by="taxon_names", all=T),
                        merge(NNMS.1_4_new, NNMS.1_4.revisedNet_new, by="taxon_names", all=T),
                        by="taxon_names", all=T)%>%
    mutate(shuffled_new = shuf_new)
  
  AllShuff_Res_new <- bind_rows(AllShuff_Res_new, shuffRes_new)
}


# one sample t test to calculate whether the NetMoss score for each node are larger than noise NetMoss score  -------------------------
library(dplyr)

dat.shuf_new <- AllShuff_Res_new


dat.1sp_new <- data.table::fread("NetMoss_new.txt", data.table = F)

results_new <- NULL

for(nd in dat.1sp_new$microbeAbb) {
  
  # 1-2 ============================
  value.1_2 = dat.1sp_new$`1-2`[which(dat.1sp_new$microbeAbb == nd)]
  shuff.1_2 = dat.shuf_new$NMSS.1_2.revised[which(dat.shuf_new$taxon_names == nd)]
  
  if (length(na.omit(shuff.1_2)) > 1) {
    wilcoxP.1_2 <- wilcox.test(shuff.1_2, mu = value.1_2, alternative = "less")$p.value
    ttestP.1_2 <- t.test(shuff.1_2, mu = value.1_2, alternative = "less")$p.value
  } else {
    wilcoxP.1_2 <- NA
    ttestP.1_2 <- NA
  }
  
  # 2-3 ============================
  value.2_3 = dat.1sp_new$`2-3`[which(dat.1sp_new$microbeAbb == nd)]
  shuff.2_3 = dat.shuf_new$NMSS.2_3.revised[which(dat.shuf_new$taxon_names == nd)]
  
  if (length(na.omit(shuff.2_3)) > 1) {
    wilcoxP.2_3 <- wilcox.test(shuff.2_3, mu = value.2_3, alternative = "less")$p.value
    ttestP.2_3 <- t.test(shuff.2_3, mu = value.2_3, alternative = "less")$p.value
  } else {
    wilcoxP.2_3 <- NA
    ttestP.2_3 <- NA
  }
  
  
  # 3-4 ============================
  value.3_4 = dat.1sp_new$`3-4`[which(dat.1sp_new$microbeAbb == nd)]
  shuff.3_4 = dat.shuf_new$NMSS.3_4.revised[which(dat.shuf_new$taxon_names == nd)]
  
  if (length(na.omit(shuff.3_4)) > 1) {
    wilcoxP.3_4 <- wilcox.test(shuff.3_4, mu = value.3_4, alternative = "less")$p.value
    ttestP.3_4 <- t.test(shuff.3_4, mu = value.3_4, alternative = "less")$p.value
  } else {
    wilcoxP.3_4 <- NA
    ttestP.3_4 <- NA
  }
  
  # 1-3 ============================
  value.1_3 = dat.1sp_new$`1-3`[which(dat.1sp_new$microbeAbb == nd)]
  shuff.1_3 = dat.shuf_new$NMSS.1_3.revised[which(dat.shuf_new$taxon_names == nd)]
  
  if (length(na.omit(shuff.1_3)) > 1) {
    wilcoxP.1_3 <- wilcox.test(shuff.1_3, mu = value.1_3, alternative = "less")$p.value
    ttestP.1_3 <- t.test(shuff.1_3, mu = value.1_3, alternative = "less")$p.value
  } else {
    wilcoxP.1_3 <- NA
    ttestP.1_3 <- NA
  }
  
  # 1-4 ============================
  value.1_4 = dat.1sp_new$`1-4`[which(dat.1sp_new$microbeAbb == nd)]
  shuff.1_4 = dat.shuf_new$NMSS.1_4.revised[which(dat.shuf_new$taxon_names == nd)]
  
  if (length(na.omit(shuff.1_4)) > 1) {
    wilcoxP.1_4 <- wilcox.test(shuff.1_4, mu = value.1_4, alternative = "less")$p.value
    ttestP.1_4 <- t.test(shuff.1_4, mu = value.1_4, alternative = "less")$p.value
  } else {
    wilcoxP.1_4 <- NA
    ttestP.1_4 <- NA
  }
  
  
  res_c_new <- data.frame(
    node = nd,
    ttestPval.1_2 = ttestP.1_2, 
    wilcoxPval.1_2 = wilcoxP.1_2, 
    ttestPval.2_3 = ttestP.2_3, 
    wilcoxPval.2_3 = wilcoxP.2_3,
    ttestPval.3_4 = ttestP.3_4, 
    wilcoxPval.3_4 = wilcoxP.3_4,
    ttestPval.1_3 = ttestP.1_3, 
    wilcoxPval.1_3 = wilcoxP.1_3,
    ttestPval.1_4 = ttestP.1_4, 
    wilcoxPval.1_4 = wilcoxP.1_4
    
  )
  
  results_new <- bind_rows(results_new, res_c_new)
}
write.table(results_new, file = 'NetMoss_oneSampleTest.w.100shuffle_new.txt', quote = F, sep = "\t", row.names = F)

#####Fig_paint
remove(list = ls())
library(data.table)
dat <- fread("NetMoss_new.txt", data.table = F)
head(dat)
# color by phylum and node significance
# see "shuffleSample_NetMoss_noise.r" for how nodes significance were calcualted 
noise <- fread("NetMoss_oneSampleTest.w.100shuffle_new.txt",data.table = F)
insig.nodes_T12 <- noise$node[which(noise$ttestPval.1_2 > 0.05) ]
T12_sig <- subset(dat,!dat$microbeAbb %in% insig.nodes_T12)
write.csv(T12_sig,"NetMoss_T12_sig_new.csv",quote = F,row.names = F)
insig.nodes_T23 <- noise$node[which(noise$ttestPval.2_3 > 0.05) ]
T23_sig <- subset(dat,!dat$microbeAbb %in% insig.nodes_T23)
write.csv(T23_sig,"NetMoss_T23_sig_new.csv",quote = F,row.names = F)
insig.nodes_T34 <- noise$node[which(noise$ttestPval.3_4 > 0.05) ]
T34_sig <- subset(dat,!dat$microbeAbb %in% insig.nodes_T34)
write.csv(T34_sig,"NetMoss_T34_sig_new.csv",quote = F,row.names = F)
insig.nodes_T13 <- noise$node[which(noise$ttestPval.1_3 > 0.05) ]
T13_sig <- subset(dat,!dat$microbeAbb %in% insig.nodes_T13)
write.csv(T13_sig,"NetMoss_T13_sig_new.csv",quote = F,row.names = F)
insig.nodes_T14 <- noise$node[which(noise$ttestPval.1_4 > 0.05) ]
T14_sig <- subset(dat,!dat$microbeAbb %in% insig.nodes_T14)
write.csv(T14_sig,"NetMoss_T14_sig_new.csv",quote = F,row.names = F)
T12 <- read.table("NetMoss_T12_sig_new.csv",header = T,sep = ",",stringsAsFactors = FALSE,check.names = F)
T23 <- read.table("NetMoss_T23_sig_new.csv",header = T,sep = ",",stringsAsFactors = FALSE,check.names = F)
T34 <- read.table("NetMoss_T34_sig_new.csv",header = T,sep = ",",stringsAsFactors = FALSE,check.names = F)
T13 <- read.table("NetMoss_T13_sig_new.csv",header = T,sep = ",",stringsAsFactors = FALSE,check.names = F)
T14 <- read.table("NetMoss_T14_sig_new.csv",header = T,sep = ",",stringsAsFactors = FALSE,check.names = F)
merge <- merge(T12,T23,by="microbeAbb",all = TRUE)
merge <- merge(merge,T34,by="microbeAbb",all = TRUE)
merge <- merge(merge,T13,by="microbeAbb",all = TRUE)
merge <- merge(merge,T14,by="microbeAbb",all = TRUE)
###paint
merge_paint <- merge
rownames(merge_paint) <- merge_paint$microbeAbb
merge_paint <- select(merge_paint,-microbeAbb)
merge_paint[is.na(merge_paint)] <- 0
merge_paint_new <- select(merge_paint,-`T2-T3`,-`T3-T4`)
####only T1-T2,T1-T3,T1-T4
merge_paint_order <- merge_paint_new
#merge_paint_order$max_value <- apply(merge_paint_order[, 1:5], 1, max)
merge_paint_order <- merge_paint_order %>%
  arrange(`T1-T2`,`T1-T3`,`T1-T4`)
merge_paint_order <- merge_paint_order[rowSums(merge_paint_order != 0) > 0, ]
rownames(mean_abundance) <- mean_abundance$taxa
merge_paint_order2 <- merge(merge_paint_order,mean_abundance,by = "row.names")
merge_paint_order2 <- merge_paint_order2 %>%
  arrange(`T1-T2`,mean_abundance,`T1-T3`,`T1-T4`)
order <- merge_paint_order2$taxa
###FigS2D
merge_paint <- merge_paint_order2
rownames(merge_paint) <- merge_paint$taxa
merge_paint <- select(merge_paint_order2,-Row.names,-taxa,-mean_abundance,-log)
rownames(merge_paint) <- merge_paint_order2$taxa
color_map <- colorRampPalette(c( "white","#FDDBC7","#F4A582","#D6604D", "#B2182B","#67001F"))(100)
order <- as.character(order)
merge_paint <- merge_paint[rev(order), , drop = FALSE]
corrplot::corrplot(merge_paint %>% as.matrix(), is.corr = F,method = "color",col=color_map,
                   outline = T,
                   tl.pos="lt", tl.cex=0.7, tl.col="black",
                   cl.pos = "r",cl.length = 11,
                   cl.ratio = 0.35,cl.offset=1,cl.cex = 1)
###FigS2E
netPlot(result = nodes_result_T1.vs.T2,
        num.top = 5,
        num.score = 20,
        e.th = 0.6,
        my.layout = layout_components,
        my.label = TRUE)