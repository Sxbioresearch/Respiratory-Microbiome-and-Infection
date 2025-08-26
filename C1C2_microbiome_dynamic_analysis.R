###microbiome  dynamic in C1 and C2

###1、Fig5C
library(ape)
data <- read.table("DNA_OP_filter_jsd_27.csv",header = T,sep = ",",row.names = 1,stringsAsFactors = FALSE)
mdata <- read.table("DNA_OP_group_27.csv",header = T,sep = ",",stringsAsFactors = FALSE)
mdata <- mdata[mdata$group1 %in% c(1,2),]
group <- read.table("C1C2_group.csv",header = T,sep = ",",stringsAsFactors = FALSE)
mdata <- merge(mdata,group,by="people")
rownames(mdata) <- mdata$sample
##PCOA
mdata <- mdata
mdata$type1 <- paste0(mdata$cluster,"_",mdata$group1)
a <- rownames(mdata)
data <- data[,a]
data_trans=t(data)
jsd <- read.table("DNA_OP_filter_jsd_27.csv",header = T,sep = ",",row.names = 1,stringsAsFactors = FALSE)
jsd <- jsd[,a]
jsd <- jsd[a,]
pcoa <- cmdscale(jsd, k = (nrow(data_trans) - 1), eig = TRUE)
eig_values <- pcoa$eig
pcoa_eig <- pcoa$eig
barplot(pcoa_eig)
pcoa_exp <- pcoa$eig/sum(pcoa$eig)
# pcoa_exp <- pcoa_eig / sum(pcoa_eig)
site <- pcoa$points
#positive_indices <- which(pcoa_eig > 0)
#site <- site[, positive_indices]

site <- data.frame(site)[1:2]
site$name <- rownames(site)
rownames(mdata)->mdata$name
#merge(shannon_index_s1,s_meta,by='SampleID')->s1_merge
merge(site,mdata,by='name')->merge
merge[order(merge$name,decreasing = T),]->merge
site[order(site$name,decreasing = T),]->site
colnames(merge)
site$group <- merge$type1
#site$type<-merge$type
# site$group<-factor(site$group,levels=c("S2","S3","S4","S5"))
pcoa1 <- paste('PCoA axis1 :', round(100*pcoa_exp[1], 2), '%')
pcoa2 <- paste('PCoA axis2 :', round(100*pcoa_exp[2], 2), '%')
aggregate(X1~group,site,mean)->mean.x
aggregate(X2~group,site,mean)->mean.y
df1<-merge(site,mean.x,by="group")
df1<-merge(df1,mean.y,by="group")
colnames(df1)<-c("group","X1","X2","name","mean.x","mean.y")
#install.packages("ggExtra")
library(ggExtra)  
library(ggsci)
library(ggplot2)
df1$group <- as.character(df1$group)
ggplot(df1, aes(X1,X2,color=group))+
  geom_point(size=4,alpha=0.8)+ 
  geom_segment(aes(x=mean.x, y=mean.y, xend=X1, yend=X2),size=0.05,alpha=1)+ 
  geom_point(aes(x=mean.x,y=mean.y),size=7,alpha=2)+
  #    scale_color_jama()+
  scale_color_manual(values = c("#868686FF","#374E55FF","#EFC000FF","#DF8F44FF"))+
  scale_x_continuous(limits = c(-0.3, 0.4))+
  scale_y_continuous(limits = c(-0.3, 0.4))+
  #scale_linetype_manual(values = "dashed")+
  # ggrepel::geom_label_repel(data=unique(select(df1,1,5,6)),aes(mean.x,mean.y,color=group),label=unique(df1$group),fontface="bold",show.legend = F,box.padding = 0,size=5)+
  theme(axis.title = element_text(size = rel(1.5), color = "black"))+
  #scale_color_manual(values =c('#91D1C2FF','#F39B7FFF','#E64B35FF','#00A087FF') )+
  # guides(fill=guide_legend(title = NULL))+
  # guides (fill = "none")+
  # theme(legend.position = "none")+
  xlab(paste("PCo1"," (",round(pcoa_exp[1]*100,1),"%",")",sep = ""))+
  ylab(paste("PCo2"," (",round(pcoa_exp[2]*100,1),"%",")",sep = ""))+
  theme_bw()+
  theme(panel.grid = element_blank())+
  theme(aspect.ratio = 1)

###2、Fig5D-E
######T1(C1 vs C2 ), T2(C1 vs C2), C1(T1 vs T2), C2(T1 vs T2)
data <- read.table('DNA_OP_filter_jsd_unique_group_27_DR.csv',header = T, sep = ',', stringsAsFactors = FALSE, check.names = FALSE)
mdata <- read.table("DNA_OP_group_27.csv",header = T,sep = ",",stringsAsFactors = FALSE)
mdata <- mdata[mdata$group1 %in% c(1,2),]
group <- read.table("C1C2_group.csv",header = T,sep = ",",stringsAsFactors = FALSE)
mdata <- merge(mdata,group,by="people")
group <- mdata
groupA <- cbind(group$sample,group$cluster) %>% as.data.frame()
names(groupA) <- c("sampleA","typeA")
data1 <- merge(data,groupA,by="sampleA")
group <- mdata
groupB <- cbind(group$sample,group$cluster) %>% as.data.frame()
names(groupB) <- c("sampleB","typeB")
data2 <- merge(data1,groupB,by="sampleB")

##C1: T2 VS T1 
data3 <-  data2[data2$typeA == data2$typeB, ]
data4 <- data3[data3$typeA == "C1",]
data5 <- data4[data4$timeA != data4$timeB,]
data5 <- data5[data5$peopleA == data5$peopleB,]

#C2: T2 VS T1
data6 <-  data2[data2$typeA == data2$typeB, ]
data7 <- data6[data6$typeA == "C2",]
data8 <- data7[data7$timeA != data7$timeB,]
data8 <- data8[data8$peopleA == data8$peopleB,]

#ref
ref <- read.table("ref_species_jsd_group.csv",header = T,sep = ",",stringsAsFactors = FALSE)
ref1 <- ref[ref$peopleA == ref$peopleB,]
ref2 <- ref1[ref1$timeA == "2020/12/14" | ref1$timeB == "2020/12/14",]
ref2$typeA <- "Reference"
ref2$typeB <- "Reference"

#T2: C1 vs C2
data9 <- data2[data2$typeA != data2$typeB, ]
data10 <- data9[data9$timeA == data9$timeB,]
data11 <- data10[data10$timeA == 2,]

#T1: C1 vs C2
data12 <- data2[data2$typeA != data2$typeB, ]
data13 <- data12[data12$timeA == data12$timeB,]
data14 <- data13[data13$timeA == 1,]
data5$group <- "C1"
data8$group <- "C2"
data11$group <- "T2"
data14$group <- "T1"
ref2$group <- "Reference"
my_comparisons <- list( c("C1", "C2"), 
                        c("C1", "T2"), 
                        c("C1", "T1"),
                        c("C2","T2"),
                        c("C2","T1"),
                        c("T2","T1")
)
data15 <- rbind(data5,data8,data11,data14,ref2)
data15 <-  data15 %>%
  mutate(type = ifelse(group %in% c("T1", "T2"), "Cluster compare", "Time compare"))
my_comparisons <- list( c("C1", "C2"), 
                        c("C1", "Reference"),
                        c("C2", "Reference"),
                        c("T1","T2")
)
library(ggpubr)
ggplot(data15, aes(x = group, y = value, color = group)) + 
  geom_violin(trim = FALSE, width=0.8) +       
  geom_jitter(width = 0.1, alpha = 0.3) + 
  scale_color_manual(values=c("C1"="#374E55FF","T2"="#B24745FF","C2"="#DF8F44FF","T1"="#374E55FF","Reference"="#79AF97FF"))+ 
  stat_compare_means(comparisons = my_comparisons,method = "wilcox.test", hide.ns = FALSE,label= "p.signif")+ 
  theme_bw() +
  xlab("Group") +
  ylab("Distance between samples") +
  facet_wrap(.~type,scales = "free_x",nrow =1)+
  theme(
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.ticks = element_line(colour = "black"),
    axis.text = element_text(colour = "black"),
    axis.title = element_text(colour = "black")
  )+
  theme(aspect.ratio = 1)

###3、Fig5F
###
group <- read.table('C1C2_group.csv',header = T,sep = ',', stringsAsFactors = FALSE, check.names = FALSE)
mdata <- read.table("DNA_OP_group_27.csv",header = T,sep = ",",stringsAsFactors = FALSE)
mdata <- merge(mdata,group,by="people")
mdata$type <- paste0(mdata$cluster,"_",mdata$group1)
rownames(mdata) <- mdata$sample
data <- read.table("DNA_species_OP_relative_filter_27.csv",header = T,sep = ",",row.names = 1,stringsAsFactors = FALSE)
data_trans=t(data)
jsd <- read.table("DNA_OP_filter_jsd_27.csv",header = T,sep = ",",row.names = 1,stringsAsFactors = FALSE)
pcoa <- cmdscale(jsd, k = (nrow(data_trans) - 1), eig = TRUE)
eig_values <- pcoa$eig
pcoa_eig <- pcoa$eig
barplot(pcoa_eig)
pcoa_exp <- pcoa$eig/sum(pcoa$eig)
# pcoa_exp <- pcoa_eig / sum(pcoa_eig)
site <- pcoa$points
#positive_indices <- which(pcoa_eig > 0)
#site <- site[, positive_indices]

site <- data.frame(site)[1:2]
site$name <- rownames(site)
rownames(mdata)->mdata$name
#merge(shannon_index_s1,s_meta,by='SampleID')->s1_merge
merge(site,mdata,by='name')->merge
merge[order(merge$name,decreasing = T),]->merge
# site[order(site$name,decreasing = T),]->site
#colnames(merge)
#site$group <- merge$type
#site$type<-merge$type
# site$group<-factor(site$group,levels=c("S2","S3","S4","S5"))
pcoa1 <- paste('PCoA axis1 :', round(100*pcoa_exp[1], 2), '%')
pcoa2 <- paste('PCoA axis2 :', round(100*pcoa_exp[2], 2), '%')
aggregate(X1~type,merge,mean)->mean.x
names(mean.x)[2] <- "mean.x"
aggregate(X2~type,merge,mean)->mean.y
names(mean.y)[2] <- "mean.y"
df1<-merge(merge,mean.x,by="type")
df1<-merge(df1,mean.y,by="type")
# colnames(df1)<-c("group","X1","X2","name","mean.x","mean.y")
#install.packages("ggExtra")
library(ggExtra)  
library(ggsci)
library(ggplot2)
df1$group1 <- as.character(df1$group1)

ggplot(df1, aes(mean.x,mean.y,color=group1,shape=cluster))+
  #   geom_point(size=4,alpha=0.7)+  
  #   geom_segment(aes(x=mean.x, y=mean.y, xend=X1, yend=X2),size=0.05,alpha=0.8)+ 
  geom_point(aes(x=mean.x,y=mean.y),size=10,alpha=2)+  
  scale_shape_manual(
    values = c(C1 = 17, C2 = 19),  
    name = "Type"               
  ) +
  geom_vline(xintercept = 0, linetype = "dashed")+
  geom_hline(yintercept = 0, linetype = "dashed")+
  #    scale_color_jama()+
  scale_color_manual(values = c("#374E55FF","#B24745FF","#DF8F44FF","#00A1D5FF"))+
  scale_x_continuous(limits = c(-0.15, 0.15))+
  scale_y_continuous(limits = c(-0.15, 0.15))+
  #scale_linetype_manual(values = "dashed")+
  # ggrepel::geom_label_repel(data=unique(select(df1,1,5,6)),aes(mean.x,mean.y,color=group),label=unique(df1$group),fontface="bold",show.legend = F,box.padding = 0,size=5)+
  theme(axis.title = element_text(size = rel(1.5), color = "black"))+
  #scale_color_manual(values =c('#91D1C2FF','#F39B7FFF','#E64B35FF','#00A087FF') )+
  # guides(fill=guide_legend(title = NULL))+
  # guides (fill = "none")+
  # theme(legend.position = "none")+
  xlab(paste("PCo1"," (",round(pcoa_exp[1]*100,1),"%",")",sep = ""))+
  ylab(paste("PCo2"," (",round(pcoa_exp[2]*100,1),"%",")",sep = ""))+
  theme_bw()+
  theme(panel.grid = element_blank())+
  theme(aspect.ratio = 1)


###4、Fig5G
data <- read.table('DNA_species_OP_relative_filter_27.csv',header = T, row.names = 1, sep = ',', stringsAsFactors = FALSE, check.names = FALSE)
new_rownames <- sub(".*\\|s__", "", rownames(data))
rownames(data) <- new_rownames
taxa <- read.table("DR_taxa_inter.csv",header = T,sep = ",",stringsAsFactors = FALSE,check.names = F)
taxa <- taxa$taxa_inter
data <- data[taxa,]
data <- as.data.frame(t(data))
minval <- min(data[data > 0])
library(tibble)
data <- rownames_to_column(data,var = "sample")
group <- read.table('DNA_OP_group_27_PCAPCB.csv',header = T,sep = ',', stringsAsFactors = FALSE, check.names = FALSE)
group_new <- read.table('C1C2_group.csv',header = T,sep = ',', stringsAsFactors = FALSE, check.names = FALSE)
group <- merge(group_new,group,by="people")
#group <- select(group,-group2,-group3,-group4,-institution,-Infection,-DMM_clusters)
group <- group %>%
  mutate(time = recode(time, `1` = "T1", `2` = "T2", `3` = "T3", `4` = "T4"))
data <- merge(group,data,by="sample")
data$people <- as.character(data$people)
data$time <- as.factor(data$time)
data$cluster <- as.factor(data$cluster)
for (i in 18:ncol(data)) {
  data[,i] <- as.numeric(data[,i])
}
data$time <- relevel(data$time, ref = "T1")
levels(data$time)
microbesofinterest <- colnames(data)[18:ncol(data)]
library(tidyverse)
library(lme4)
library(broom)
#install.packages("broom.mixed")
library(broom.mixed)
library(reshape2)
#install.packages("lmerTest")
library(lmerTest)
####C1 group
data_C1 <- data[data$cluster == "C1",]
# 回归分析
regression_output_C1 <- list()
for(m in microbesofinterest) {
  print(m)
  regression_output_C1[[m]] <- try(
    lmer(data = data_C1, log(data_C1[[m]] + minval) ~ time + (1 | people)) %>% 
      tidy() %>% 
      mutate(yvar = m) %>% 
      filter(term != '(Intercept)'),
    silent = TRUE
  )
}
regression_output_C1[[1]]
regression_output_use_C1 = bind_rows(regression_output_C1)
rregression_output_use_C1 = regression_output_use_C1 %>% mutate(BH_adjusted = p.adjust(p.value,method='BH'))
write.csv(rregression_output_use_C1,"regression_output_lme_C1.csv",quote=F,row.names = F)
####C2 group
data_C2 <- data[data$cluster == "C2",]
regression_output_C2 <- list()
for(m in microbesofinterest) {
  print(m)
  regression_output_C2[[m]] <- try(
    lmer(data = data_C2, log(data_C2[[m]] + minval) ~ time + (1 | people)) %>% 
      tidy() %>% 
      mutate(yvar = m) %>% 
      filter(term != '(Intercept)'),
    silent = TRUE
  )
}

regression_output_C2[[1]]
regression_output_use_C2 = bind_rows(regression_output_C2)
rregression_output_use_C2 = regression_output_use_C2 %>% mutate(BH_adjusted = p.adjust(p.value,method='BH'))
write.csv(rregression_output_use_C2,"regression_output_lme_C2.csv",quote=F,row.names = F)
####diff species
C1 <- read.table("regression_output_lme_C1.csv",header = T,sep = ",",stringsAsFactors = FALSE)
C1$type <- "C1"
C2 <- read.table("regression_output_lme_C2.csv",header = T,sep = ",",stringsAsFactors = FALSE)
C2$type <- "C2"
C1C2 <- as.data.frame(rbind(C1,C2))
C1C2 <- C1C2[!is.na(C1C2$BH_adjusted), ]
C1C2_sig <- C1C2[C1C2$BH_adjusted < 0.05,]
C1C2_sig$time <- paste0(C1C2_sig$type,"_",C1C2_sig$term)
library(dplyr)
C1C2_sig_paint  <- C1C2_sig %>%
  arrange(estimate)
a <- C1C2_sig_paint$yvar
a <- unique(a)
C1C2_sig_paint$yvar <- factor(C1C2_sig_paint$yvar,levels = rev(a))
library(ggplot2)
ggplot(C1C2_sig_paint, aes(x = time, y = yvar, fill = estimate)) +
  geom_tile(color = "white") + 
  # scale_fill_viridis() + 
  scale_fill_gradientn(colours = c("#053061", "#2166AC", "#4393C3", "#92C5DE","#D1E5F0","#FFFFFF","#FDDBC7","#F4A582","#D6604D", "#B2182B","#67001F"),
                       values = scales::rescale(c(-2,-1.5,-1,-0.5,-0.1,0,0.1,0.5,1,1.5,2,2.5))
  )+
  geom_text(aes(label = ifelse(BH_adjusted < 0.001, "***", 
                               ifelse(BH_adjusted < 0.01, "**", 
                                      ifelse(BH_adjusted < 0.05, "*", 
                                             ifelse(BH_adjusted < 0.1, ".", ""))))), size = 6, color = "white",vjust = 0.8
  ) + 
  theme_minimal() +  
  labs(x = "Time", y = "Phenotype", fill = "Coefficient")+ 
  theme(
    # panel.border = element_blank(),   
    axis.ticks = element_blank(),      
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 0.5)
  )

