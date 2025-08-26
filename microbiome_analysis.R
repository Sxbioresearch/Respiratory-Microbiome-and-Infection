###1、decontamination
###NP decontamination (use NC + NC-NP together)
library(phyloseq)
library(ggplot2)
library(decontam)
DNA <- read.table("dna_S_all.csv",header = T,sep = ",",row.names = 1)
DNA_t <- as.data.frame(t(DNA))
library(tibble)
DNA_t <- rownames_to_column(DNA_t, var = "sample")
group <- read.table("DNA_group_all.csv",header = T,sep = ",")
DNA_merge <- merge(DNA_t,group,by = "sample",all = T)
DNA_NP <- subset(DNA_merge, group3 %in% c("NP", "NC","NC-NP"))
#NP <- subset(DNA_merge, group3 %in% c("NP"))
#NP <- select(NP,-group1,-group2,-group3,-group4,-people,-institution,-Infection)
#rownames(NP) <- NP$sample
#NP <- select(NP,-sample)
#NP <- as.data.frame(t(NP))
#df <- NP[rowSums(NP != 0) > 0, ]
library(dplyr)
DNA_NP <- select(DNA_NP,-group1,-group2,-group3,-group4,-people,-institution,-Infection)
DNA_NP2 <- as.data.frame(t(DNA_NP))
write.table(DNA_NP2,"DNA_NP_raw_all.csv",quote = F,sep=",",col.names = F)
load("lab_metadata.RData")
lab_info_dna$X16S.2.ng.ul. <- as.numeric(lab_info_dna$X16S.2.ng.ul.)
lab_info_dna$X16S.2.ng.ul.[is.na(lab_info_dna$X16S.2.ng.ul.)] <- 0
lab_info_dna$type <- lab_info_dna$样本类型2
lab_info_dna$type[lab_info_dna$样本类型1 == "NC"] <- "NC"
metadecontam <- as.data.frame(cbind(lab_info_dna$seqID,lab_info_dna$X16S.2.ng.ul.,lab_info_dna$样本体积,lab_info_dna$type))
names(metadecontam) <- c("sample","X16S","volume","type")
metadecontam$X16S <- as.numeric(metadecontam$X16S)
metadecontam$volume <- as.numeric(metadecontam$volume)
metadecontam$conc <- metadecontam$X16S*metadecontam$volume
write.csv(metadecontam,"metadecontam_NP.csv",quote = F,row.names = F)
metadecontam <- read.table("metadecontam_NP.csv",header = T,sep = ",",quote = "")
library(dplyr)
metadecontam <- select(metadecontam,-X16S,-volume)
data <- merge(DNA_NP,metadecontam,by = "sample",all.x = TRUE)
DNA_st <- select(data,-type,-conc)
rownames(DNA_st) <- DNA_st$sample
DNA_st <- select(DNA_st,-sample)
DNA_st <- as.matrix(DNA_st)
data$group <- data$type
data$group[data$group == "NC"] <- TRUE
data$group[data$group != "TRUE"] <- FALSE
DNA_conc <- data$conc
DNA_neg <- data$group
DNA_conc_new5 <- log(1 + DNA_conc) / log(1 + max(DNA_conc))
DNA_conc_new6 <- DNA_conc_new5+1
DNA_neg <- as.logical(DNA_neg)
decontam_DNA <- isContaminant(DNA_st, conc=DNA_conc_new6, neg=DNA_neg, method="both",threshold=c(0.5,0.500001),normalize =TRUE)
write.csv(decontam_DNA,"decontam_DNA_NP_all.csv",quote = F)

###quantile test
data <- read.table("DNA_NP_raw_all.csv",header = T,sep = ",",row.names = 1)
library(tidyr)
data=data %>% apply(.,2,function(x){x/sum(x)})
dna_meta1 <- read.table("DNA_group_all.csv",header = T,sep = ",")
dna_meta1 <- dna_meta1[dna_meta1$sample != "D19",]
dna_meta1 <- subset(dna_meta1,dna_meta1$group3 %in% c("NC","NP","NC-NP"))
nc_species_rela= data[,dna_meta1$sample[dna_meta1$group2=="NC"]] %>% .[rowSums(.)!=0,colSums(.)!=0]
sample_species_rela=data[,dna_meta1$sample[dna_meta1$group2!="NC"]] %>% .[rowSums(.)!=0,colSums(.)!=0]
species_test=intersect(rownames(nc_species_rela),rownames(sample_species_rela))
two_sample_quantail_test=function(x){
  s1=sample_species_rela[x,] %>% as.numeric() %>% na.omit()
  s2=nc_species_rela[x,] %>% as.numeric() %>% na.omit()
  test=EnvStats::quantileTest(x = s1,y = s2,
                              target.quantile = 0.95,
                              exact.p = T,
                              alternative = "greater")
  return(test$p.value)
}
p_value=lapply(species_test,two_sample_quantail_test)
hist(p_value %>% unlist)
species_contami=data.frame(species=species_test,p_value=unlist(p_value))

####2、top_genus
###物种组成(genus)——corplot
##目标物种（每种里面的top-5放在一起画heatmap）
#DNA-NP
data <- read.table("DNA_genus_NP_relative_filter_all.csv",header = T,sep = ",",row.names = 1,stringsAsFactors = FALSE) #data=species_decontam
#data <- data[rowSums(data != 0) > (ncol(data) / 2), ]
library(tidyr)
tmp1=order(apply(data,1,mean),decreasing = T)[1:5] %>% rownames(data)[.]
main_taxa=apply(data,1,max)
#main_taxa=names(main_taxa[main_taxa>0.15])
#tmp1=unique(c(tmp1,main_taxa))
tmp=data[tmp1,] %>% as.data.frame()
tmp$row_mean <- rowMeans(tmp, na.rm = TRUE)
DNA_NP <- rownames(tmp)
#DNA-OP
data <- read.table("DNA_genus_OP_relative_filter.csv",header = T,sep = ",",row.names = 1,stringsAsFactors = FALSE) #data=species_decontam
library(tidyr)
tmp1=order(apply(data,1,mean),decreasing = T)[1:5] %>% rownames(data)[.]
#main_taxa=apply(data,1,max)
#main_taxa=names(main_taxa[main_taxa>0.15])
#tmp1=unique(c(tmp1,main_taxa))
tmp=data[tmp1,] %>% as.data.frame()
tmp$row_mean <- rowMeans(tmp, na.rm = TRUE)
DNA_OP <- rownames(tmp)
#RNA-NP
data <- read.table("RNA_genus_NP_relative_filter.csv",header = T,sep = ",",row.names = 1,stringsAsFactors = FALSE) #data=species_decontam
library(tidyr)
tmp1=order(apply(data,1,mean),decreasing = T)[1:5] %>% rownames(data)[.]
tmp=data[tmp1,] %>% as.data.frame()
tmp$row_mean <- rowMeans(tmp, na.rm = TRUE)
RNA_NP <- rownames(tmp)
#RNA-OP
data <- read.table("RNA_genus_OP_relative_filter_new.csv",header = T,sep = ",",row.names = 1,stringsAsFactors = FALSE) #data=species_decontam
library(tidyr)
tmp1=order(apply(data,1,mean),decreasing = T)[1:5] %>% rownames(data)[.]
#main_taxa=apply(data,1,max)
#main_taxa=names(main_taxa[main_taxa>0.15])
#tmp1=unique(c(tmp1,main_taxa))
tmp=data[tmp1,] %>% as.data.frame()
tmp$row_mean <- rowMeans(tmp, na.rm = TRUE)
RNA_OP <- rownames(tmp)
###target genus
genus_target <- c(DNA_NP,DNA_OP,RNA_NP,RNA_OP)
genus_target <- unique(genus_target)
###DNA_NP
DNA_NP_data <- read.table("DNA_NP_genus_relative_decontam.csv",header = T,sep = ",",row.names = 1,stringsAsFactors = FALSE)
DNA_NP_data_use <- DNA_NP_data[genus_target,]
DNA_NP_data_use <- DNA_NP_data_use[complete.cases(DNA_NP_data_use), ]
DNA_NP_paint <- rowMeans(DNA_NP_data_use) %>% as.data.frame()
names(DNA_NP_paint) <- "DNA_NP"
###DNA_OP
DNA_OP_data <- read.table("DNA_OP_genus_relative_decontam.csv",header = T,sep = ",",row.names = 1,stringsAsFactors = FALSE)
DNA_OP_data_use <- DNA_OP_data[genus_target,]
DNA_OP_data_use <- DNA_OP_data_use[complete.cases(DNA_OP_data_use), ]
DNA_OP_paint <- rowMeans(DNA_OP_data_use) %>% as.data.frame()
names(DNA_OP_paint) <- "DNA_OP"
###RNA_NP
RNA_NP_data <- read.table("RNA_NP_genus_relative_decontam.csv",header = T,sep = ",",row.names = 1,stringsAsFactors = FALSE)
RNA_NP_data_use <- RNA_NP_data[genus_target,]
RNA_NP_data_use <- RNA_NP_data_use[complete.cases(RNA_NP_data_use), ]
RNA_NP_paint <- rowMeans(RNA_NP_data_use) %>% as.data.frame()
names(RNA_NP_paint) <- "RNA_NP"
###RNA_OP
RNA_OP_data <- read.table("RNA_OP_genus_relative_decontam.csv",header = T,sep = ",",row.names = 1,stringsAsFactors = FALSE)
RNA_OP_data_use <- RNA_OP_data[genus_target,]
RNA_OP_data_use <- RNA_OP_data_use[complete.cases(RNA_OP_data_use), ]
RNA_OP_paint <- rowMeans(RNA_OP_data_use) %>% as.data.frame()
names(RNA_OP_paint) <- "RNA_OP"
merge1 <- merge(DNA_NP_paint,DNA_OP_paint,by="row.names",all=TRUE)
merge2 <- merge(RNA_NP_paint,RNA_OP_paint,by="row.names",all=TRUE)
merge_paint <- merge(merge1,merge2,by="Row.names",all=TRUE)
merge_paint[is.na(merge_paint)] <- 0
merge_paint1 <- merge_paint
rownames(merge_paint1) <- merge_paint1$Row.names
merge_paint1 <- select(merge_paint1,-Row.names)
merge_paint2 <- merge_paint1
merge_paint2$row_mean <- rowMeans(merge_paint2, na.rm = TRUE)
merge_paint2 <- merge_paint2 %>%
  arrange(row_mean)
order <- rownames(merge_paint2)
order <- rev(order)
merge_paint3 <- merge_paint1[order,]
color_map <- colorRampPalette(c("#FFFFFF","#FDDBC7","#F4A582","#D6604D", "#B2182B","#67001F"))(100)
###Fig1C
corrplot::corrplot(merge_paint3 %>% as.matrix(), is.corr = F,method = "circle",col=color_map,
                   
                   outline = T,
                   #  addCoef.col = 'black',
                   tl.pos="lt", tl.cex=1, tl.col="black",
                   cl.pos = "b",cl.length = 5,
                   cl.ratio = 0.05,cl.offset=0.3,cl.cex = 0.8)


###3、alpha diversity
library(vegan)
alpha_diversity <- function(x) {
  observed_species <- estimateR(x)[1, ]
  Chao1 <- estimateR(x)[2, ]
  ACE <- estimateR(x)[4, ]
  Shannon <- diversity(x, index = 'shannon',base = 2)
  Simpson <- diversity(x, index = 'simpson')    
  goods_Coverage <- 1 - rowSums(x == 1) / rowSums(x)
  result <- data.frame(observed_species, ACE,Chao1, Shannon, Simpson, goods_Coverage)
  result
}
data1 <- read.table("DNA_OP_genus_filter_raw.csv",header = T,sep = ",",row.names = 1)
data2 <- t(data1)
alpha_OP <- alpha_diversity(data2)

###4、beta diversity
library(phyloseq)
#jsd distance
data <- read.table("DNA_genus_OP_relative_filter_27.csv",header = T,sep = ",",row.names = 1,stringsAsFactors = FALSE)
data1 <- t(data)
data2 = otu_table(data1, taxa_are_rows = F)
data3 = phyloseq(data2)
jsd=phyloseq::distance(data3, method = "jsd")
jsd1 <- sqrt(jsd)
jsd2 <- as.matrix(jsd1)
###Fig2C
mdata <- read.table("DNA_OP_group_27.csv",header = T,sep = ",",row.names = 1,stringsAsFactors = FALSE)
data <- read.table("DNA_genus_OP_relative_filter_27.csv",header = T,sep = ",",row.names = 1,stringsAsFactors = FALSE)
data_trans=t(data)
jsd <- read.table("DNA_OP_genus_filter_jsd_27.csv",header = T,sep = ",",row.names = 1,stringsAsFactors = FALSE)
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
site$group <- merge$group1
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
  geom_point(size=4,alpha=0.5)+  
  geom_segment(aes(x=mean.x, y=mean.y, xend=X1, yend=X2),size=0.05,alpha=0.7)+ 
  geom_point(aes(x=mean.x,y=mean.y),size=7,alpha=2)+ 
  #    scale_color_jama()+
  scale_color_manual(values = c("#374E55FF","#B24745FF","#DF8F44FF","#00A1D5FF"))+
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

###5、differential analysis
###zicoseq method
### T1 vs T2
data <- read.table("DNA_genus_OP_relative_filter_27.csv",header = T,sep = ",",row.names = 1,stringsAsFactors = FALSE)
group <- read.table("DNA_OP_group_27.csv",header = T,sep = ",",stringsAsFactors = FALSE)
group <- group[group$group1 != 3,]
group <- group[group$group1 != 4,]
a <- group$sample
input1 <- data[,a]
meta1 <- group
library(GUniFrac)
meta1$group1 <- as.character(meta1$group1)
ZicoSeq1 <- ZicoSeq(meta.dat = meta1, feature.dat = as.matrix(input1), 
                    grp.name = 'group1',feature.dat.type = "proportion",
                    # Winsorization to replace outliers
                    is.winsor = TRUE, outlier.pct=0.01,winsor.end = 'top',
                    # Posterior sampling 
                    is.post.sample = TRUE,
                    # Use the square-root transformation
                    link.func = list(function (x) x^0.5),stats.combine.func = max,
                    # Permutation-based multiple testing correction
                    perm.no = 999,  strata =  meta1$people , 
                    # Reference-based multiple stage normalization
                    ref.pct = 0.5, stage.no = 6, excl.pct = 0.2,
                    # Family-wise error rate control
                    is.fwer = TRUE, verbose = TRUE, return.feature.dat = TRUE)
ZicoSeq.plot(ZicoSeq1,
             pvalue.type = c('p.adj.fdr'),
             cutoff = 0.1)
##foldchange
Fd_an1 <- as.data.frame(colnames(input1))
names(Fd_an1) <- "sample"
library(dplyr)
Fd_an1 <- left_join(Fd_an1,meta1)
Fd_an1 <- Fd_an1 %>%
  arrange(people,group1)
b <- Fd_an1$sample
input1 <- input1[,b]
Fd_an1 <- Fd_an1$group1

min = min(input1[input1!=min(input1)]) 
for (i in 1:nrow(input1)){
  for (j in 1:ncol(input1)) {
    if(input1[i,j]==0)
      input1[i,j] <- runif(1,0.1*min,min)
  }
  
}

foldChange <- function(inData, classLabel) {
  sampleIdsCase <- which(classLabel == "2")  
  sampleIdsControl <- which(classLabel == "1") 
  if (length(sampleIdsCase) != length(sampleIdsControl)) {
    stop("配对样本数量不相等。")
  }
  probeFC <- rep(0, nrow(inData))
  for (i in 1:nrow(inData)) {
    ratios <- numeric(length(sampleIdsCase))  
    for (j in 1:length(sampleIdsCase)) {
      caseValue <- as.numeric(inData[i, sampleIdsCase[j]])
      controlValue <- as.numeric(inData[i, sampleIdsControl[j]])
      ratios[j] <- caseValue / controlValue
    }
    probeFC[i] <- median(ratios)  
  }
  probeFC <- log(probeFC, base = 2)
  return(probeFC)
}
Fd_r1 <- as.data.frame(foldChange(input1,Fd_an1))
rownames(Fd_r1) <- rownames(input1)
names(Fd_r1) <- "Log2FC"
library(tibble)
Fd_r1 <- rownames_to_column(Fd_r1,var="taxa")
names(Fd_r1)[2] <- "Log2FC(T2/T1)"
ZicoSeq1_p <- as.data.frame(ZicoSeq1$p.adj.fdr)
ZicoSeq1_p <- tibble::rownames_to_column(ZicoSeq1_p, var = "taxa")
ZicoSeq1_p$R2 <- ZicoSeq1$R2
ZicoSeq1_p <- left_join(ZicoSeq1_p,Fd_r1)
names(ZicoSeq1_p)[1:3] <- c("taxa","p.adj.fdr","R2")
ZicoSeq1_p$log10P <- log10(ZicoSeq1_p$p.adj.fdr)
#average relative abundance
input1 <- data[,a]
mean_abundance <- apply(input1, 1, mean) %>% as.data.frame()
mean_abundance <- tibble::rownames_to_column(mean_abundance, var = "taxa")
names(mean_abundance) <- c("taxa","mean_abundance")
#prevalence
prevalence <- apply(input1, 1, function(row) sum(row != 0)) %>% as.data.frame() 
prevalence <- tibble::rownames_to_column(prevalence, var = "taxa")
names(prevalence) <- c("taxa","prevalence")
prevalence$prevalence <- prevalence$prevalence/54
ZicoSeq1_p <- left_join(ZicoSeq1_p,mean_abundance) 
ZicoSeq1_p <- left_join(ZicoSeq1_p,prevalence) 
write.csv(ZicoSeq1_p,"DNA_OP_genus_zicoseq_T12.csv",quote = F,row.names = F)

###lme method
data <-  read.table('DNA_genus_OP_relative_filter_27.csv',header = T, row.names = 1, sep = ',', stringsAsFactors = FALSE, check.names = FALSE)
data1 <- as.data.frame(t(data))
data1 <- rownames_to_column(data1,var = "sample")
group <- read.table("DNA_OP_group_27.csv",header = T,sep = ",",stringsAsFactors = FALSE)
data2 <- merge(group,data1,by="sample")
data2 <- data2 %>%
  mutate(group1 = recode(group1, `1` = "T1", `2` = "T2", `3` = "T3", `4` = "T4"))
#data3 <- data2[data2$group1 != "T4",]
data2 <- select(data2,-group2,-group3,-group4,-institution,-Infection)
data3 <- select(data2,-sample,-group1,-people)
microbesofinterest <- names(data3)
minval <- min(data3[data3 > 0])
names(data2)[2:3] <- c("Time","People")
DNA_OP <- data2
DNA_OP$Time <- as.factor(DNA_OP$Time)
DNA_OP$Time <- relevel(DNA_OP$Time, ref = "T1")
levels(DNA_OP$Time)
library(tidyverse)
library(lme4)
library(broom)
#install.packages("broom.mixed")
library(broom.mixed)
library(reshape2)
#install.packages("lmerTest")
library(lmerTest)
regression_output_OP <- list()
for(m in microbesofinterest) {
  print(m)
  regression_output_OP[[m]] <- try(
    lmer(data = DNA_OP, log(DNA_OP[[m]] + minval) ~ Time + (1 | People)) %>% 
      tidy() %>% 
      mutate(yvar = m) %>% 
      filter(term != '(Intercept)'),
    silent = TRUE
  )
}
regression_output_OP[[1]]
regression_output_OP_use = bind_rows(regression_output_OP)
rregression_output_OP_use = regression_output_OP_use %>% mutate(BH_adjusted = p.adjust(p.value,method='BH'))
write.csv(rregression_output_OP_use,"regression_output_OP.csv",quote=F)
###Fig2G
high_paint <- read.table("zicoseq_lme_inter_high.csv",header = T,sep = ",",stringsAsFactors = FALSE)
high_paint$group4 <- ifelse(high_paint$Log2FC > 0, "increase", "decrease")
high_paint$abs <- abs(high_paint$Log2FC)
high_paint$taxa <- factor(high_paint$taxa,levels = high_target)
ggplot(high_paint,aes(x = group1, 
                      y = taxa,
                      size = abs,
                      #alpha = p,
                      shape=group4
)) +
  geom_point(aes(fill = group2), colour='#4C4C4C',stroke = 0.25)+
  # scale_fill_nejm()+
  scale_fill_manual(values=c(zicoseq='#C16622FF',lme='#FFA319FF',all='#800000FF',no='grey'))+
  #  scale_alpha_continuous(range = c(1, 0), limits = c(0, 0.23)) +  # 设置透明度的范围和限制
  scale_shape_manual(values = c(25, 24))+
  labs(x = "Time", y = "Taxa")+           
  # scale_radius(                               
  #   range=c(2,6),                             
  #  name="log(sum)")+                             
  #scale_size_continuous(range = c(1, 5))+
  theme_bw() +                               
  theme(axis.text.x = element_text(angle = 30,vjust = 0.85,hjust = 0.75))+
  theme(panel.grid = element_blank())+
  theme(aspect.ratio = 3)
####fitted curve
data <- read.table('DNA_genus_OP_relative_filter_27.csv',header = T,sep = ',',row.names = 1, stringsAsFactors = FALSE, check.names = FALSE)
group <- read.table("DNA_OP_group_27.csv",header = T,sep = ",",stringsAsFactors = FALSE)
data <- as.data.frame(t(data))
library(tibble)
data <- rownames_to_column(data,var = "sample")
merge <- merge(group,data,by="sample")
#merge <- merge %>%
# mutate(group1 = recode(group1, `1` = "T1", `2` = "T2", `3` = "T3", `4` = "T4"))
merge <- select(merge,-sample,-group2,-group3,-group4,-institution,-Infection,-Infection)
library(reshape2)
merge$group1 <- as.character(merge$group1)
merge_long <- melt(merge,by="group1")
merge_c3 <- subset(merge_long,merge_long$variable %in% high_target)
merge_c3$variable <- factor(merge_c3$variable,levels = rev(high_target))
library(ggplot2)
ggplot(merge_c3, aes(x=group1, y=value,group=1)) +
  # geom_point() +
  # geom_line()
  geom_smooth(method = "loess")+
  facet_wrap(.~variable,scales = "free_y")+
  ylab("Relative abundance")+
  theme_bw()
###GAMM fitted curve
library(mgcv)
merge_c3$group1 <- as.numeric(merge_c3$group1)
unique_values <- length(unique(merge_c3$group1))
print(paste("Unique values in group1:", unique_values))
k_value <- min(10, unique_values)
split_data <- split(merge_c3, merge_c3$variable)
library(mgcv) 
library(dplyr)
dt <- NULL
fit <- NULL
ano <- NULL
ano_p <- list()
for (i in 1:7) {
  dt[[i]] <- split_data[[i]]
  dt[[i]]$group1 <- as.numeric(dt[[i]]$group1)
  fit[[i]]<-gamm(value~s(group1, k = k_value),random=list(people=~1),data=dt[[i]],family=gaussian)  ###拟合gamm模型
  ano[[i]] <- anova(fit[[i]]$gam)
  ano_p[[i]] <- ano[[i]]$s.pv
}
ano_p <- unlist(ano_p) %>% as.data.frame()
rownames(ano_p) <- names(split_data)
names(ano_p) <- "p"
ano_p$BH <- p.adjust(ano_p$p,method = "BH")
for (i in 1:7) {
  plot(fit[[i]]$gam,xaxt = "n", main =paste(names(split_data)[i]))
  axis(1, at = 1:4, labels = 1:4)  
}
par(mfrow=c(3,3));for (i in 1:7) {plot(fit[[i]]$gam,xaxt = "n", main =paste(names(split_data)[i])) 
  axis(1, at = 1:4, labels = 1:4)}


