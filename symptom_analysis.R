###Comparison of symptoms between C1 and C2 group

###1、overall symptom
###number
data <- read.table("overall_symptom_27_all_new.csv",header = T,sep = ",",stringsAsFactors = FALSE)
data_C1 <- data[data$cluster == "C1",]
data_C2 <- data[data$cluster == "C2",] 
##c1
counts <- numeric(15)
for (k in 0:14) {
  counts[k + 1] <- sum(data_C1$number > k, na.rm = TRUE)
}
result_C1 <- data.frame(k = 0:14, count = counts)
##C2
counts <- numeric(15)
for (k in 0:14) {
  counts[k + 1] <- sum(data_C2$number > k, na.rm = TRUE)
}
result_C2 <- data.frame(k = 0:14, count = counts)
result <- merge(result_C1,result_C2,by="k")
names(result) <- c("number","C1_number","C2_number")
result$C1_ratio <- result$C1_number/14
result$C2_ratio <- result$C2_number/13
data <- result
data1 <- select(data,-C1_number,-C2_number)
data2 <- melt(data1,id="number")
library(ggalt)
###Fig6A
ggplot(data2,aes(x=number,y=value,group=variable))+
  geom_smooth(aes(color=variable),se=FALSE, size=1.5) +
  geom_point(aes(fill=variable),size=3,shape=21)+
  # geom_xspline(spline_shape = 1)+
  scale_color_manual(values=c("#374E55FF","#DF8F44FF"))+
  scale_fill_manual(values=c("#374E55FF","#DF8F44FF"))+
  xlab("Symptom number")+ylab("The proportion of individual")+
  theme_bw()+
  theme(aspect.ratio = 1)
ks_test_result <- ks.test(data$C1_ratio, data$C2_ratio,alternative =  "greater")###0.1832


### score
data <- read.table("overall_symptom_27_all_new.csv",header = T,sep = ",",stringsAsFactors = FALSE)
data_C1 <- data[data$cluster == "C1",]
data_C2 <- data[data$cluster == "C2",] 
counts <- numeric(36)
for (k in 0:35) {
  counts[k + 1] <- sum(data_C1$score > k, na.rm = TRUE)
}
result_C1 <- data.frame(k = 0:35, count = counts)
counts <- numeric(36)
for (k in 0:35) {
  counts[k + 1] <- sum(data_C2$score > k, na.rm = TRUE)
}
result_C2 <- data.frame(k = 0:35, count = counts)
result <- merge(result_C1,result_C2,by="k")
names(result) <- c("score","C1_number","C2_number")
result$C1_ratio <- result$C1_number/14
result$C2_ratio <- result$C2_number/13
data <- result_new
data1 <- select(data,-C1_number,-C2_number)
data2 <- melt(data1,id="score")
library(ggalt)
###Fig6B
ggplot(data2,aes(x=score,y=value,group=variable))+
  geom_smooth(aes(color=variable),se=FALSE, size=1.5) +
  geom_point(aes(fill=variable),size=3,shape=21)+
  scale_color_manual(values=c("#374E55FF","#DF8F44FF"))+
  scale_fill_manual(values=c("#374E55FF","#DF8F44FF"))+
  xlab("Symptom score")+ylab("The proportion of individual")+
  theme_bw()+
  theme(aspect.ratio = 1)
ks_test_result <- ks.test(data$C1_ratio, data$C2_ratio,alternative =  "greater")###0.02659

###duration
duration <- read.table("duration_40_last.csv",header = T,sep = ",",stringsAsFactors = FALSE)
meta <- read.table("symptom_meta_27.csv",header = T,sep = ",",stringsAsFactors = FALSE)
people <- meta$ID
data <- duration[duration$ID %in% people,]
data1 <- melt(data,id="ID")
data2 <- merge(meta,data1,by="ID")
data_C1 <- data2[data2$cluster == "C1",]
data_C2 <- data2[data2$cluster == "C2",] 
counts <- numeric(27)
for (k in 0:26) {
  counts[k + 1] <- sum(data_C1$value > k, na.rm = TRUE)
}
result_C1 <- data.frame(k = 0:26, count = counts)
counts <- numeric(27)
for (k in 0:26) {
  counts[k + 1] <- sum(data_C2$value > k, na.rm = TRUE)
}
result_C2 <- data.frame(k = 0:26, count = counts)
result <- merge(result_C1,result_C2,by="k")
names(result) <- c("time","C1_number","C2_number")
result$C1_ratio <- result$C1_number/266
result$C2_ratio <- result$C2_number/247
data <- result
data1 <- select(data,-C1_number,-C2_number)
data2 <- melt(data1,id="time")
library(ggalt)
##FigS8A
ggplot(data2,aes(x=time,y=value,group=variable))+
  geom_smooth(aes(color=variable),se=FALSE, size=1.5) +
  geom_point(aes(fill=variable),size=3,shape=21)+
  # geom_xspline(spline_shape = 1)+
  scale_color_manual(values=c("#374E55FF","#DF8F44FF"))+
  scale_fill_manual(values=c("#374E55FF","#DF8F44FF"))+
  xlab("Symptom_duration")+ylab("The proportion of individual")+
  theme_bw()+
  theme(aspect.ratio = 1)
ks_test_result <- ks.test(data$C1_ratio, data$C2_ratio,alternative =  "less") ##0.0778

###2、single symptom
##number
##fisher test
C1C2_phenotype <- read.table("dichotomy_number.csv",header = T,sep = ",",stringsAsFactors = FALSE)
results <- data.frame(p_value = numeric(nrow(C1C2_phenotype)), OR = numeric(nrow(C1C2_phenotype)))
data <- NULL
for (i in 1:nrow(C1C2_phenotype)) {
  data[[i]] <- matrix(c(C1C2_phenotype[i, 2], C1C2_phenotype[i, 3],(14 - C1C2_phenotype[i, 2]),(13 - C1C2_phenotype[i, 3])), nrow = 2)
  rownames(data[[i]]) <- c("C1", "C2")
  colnames(data[[i]]) <- c("exist", "absence")
  fisher_test <- fisher.test(data[[i]])
  results$p_value[i] <- fisher_test$p.value
  results$OR[i] <- fisher_test$estimate
}
results$BH <- p.adjust(results$p_value, method = "BH")
rownames(results) <- C1C2_phenotype$X
write.csv(results,"dichotomy_fisher.csv",quote = F)
###Fig6C
dichotomy <- read.table("dichotomy_number.csv",header = T,sep = ",",stringsAsFactors = FALSE)
dichotomy$order <- dichotomy$C2_ratio-dichotomy$C1_ratio
dichotomy <- dichotomy %>%
  arrange(order)
order <- dichotomy$X
paint <- select(dichotomy,-C1,-C2,-order)
paint1 <- melt(paint,id="X")
paint1$X <- factor(paint1$X,levels = order)
paint1$variable <- factor(paint1$variable,levels = c("C2_ratio","C1_ratio"))
pdf("single_number.pdf",height = 7,width = 5.3)
ggplot(paint1, aes(x =value , y = X, fill = variable)) 
  geom_bar(stat = "identity", position = position_dodge(),size=0.5) +
  scale_fill_manual(values = c("C2_ratio"="#DF8F44FF","C1_ratio"="#374E55FF"))+
  scale_x_continuous(limits = c(0,1))+
  labs(x = "Symptom Prevalence (%)", y = "") + 
  theme_bw()
###C1 VS C2
data <- matrix(c(3,16,19-3,19-16), nrow = 2)
rownames(data) <- c("C1", "C2")
colnames(data) <- c("more", "less")
fisher_test <- fisher.test(data)


##score（>2 vs <=2）
##fisher test
C1C2_phenotype <- read.table("severity_3.csv",header = T,sep = ",",stringsAsFactors = FALSE)
results <- data.frame(p_value = numeric(nrow(C1C2_phenotype)), OR = numeric(nrow(C1C2_phenotype)))
data <- NULL
for (i in 1:nrow(C1C2_phenotype)) {
  data[[i]] <- matrix(c(C1C2_phenotype[i, 2], C1C2_phenotype[i, 3],(14 - C1C2_phenotype[i, 2]),(13 - C1C2_phenotype[i, 3])), nrow = 2)
  rownames(data[[i]]) <- c("C1", "C2")
  colnames(data[[i]]) <- c("exist", "absence")
  fisher_test <- fisher.test(data[[i]])
  results$p_value[i] <- fisher_test$p.value
  results$OR[i] <- fisher_test$estimate
}
results$BH <- p.adjust(results$p_value, method = "BH")
rownames(results) <- C1C2_phenotype$X
write.csv(results,"severity_fisher.csv",quote = F)
C1C2_phenotype <- read.table("severity_3.csv",header = T,sep = ",",stringsAsFactors = FALSE)
results <- data.frame(p_value = numeric(nrow(C1C2_phenotype)), OR = numeric(nrow(C1C2_phenotype)))
data <- NULL
for (i in 1:nrow(C1C2_phenotype)) {
  data[[i]] <- matrix(c(C1C2_phenotype[i, 2], C1C2_phenotype[i, 4],(C1C2_phenotype[i, 3] - C1C2_phenotype[i, 2]),(C1C2_phenotype[i, 5] - C1C2_phenotype[i, 4])), nrow = 2)
  rownames(data[[i]]) <- c("C1", "C2")
  colnames(data[[i]]) <- c("exist", "absence")
  fisher_test <- fisher.test(data[[i]])
  results$p_value[i] <- fisher_test$p.value
  results$OR[i] <- fisher_test$estimate
}
results$BH <- p.adjust(results$p_value, method = "BH")
rownames(results) <- C1C2_phenotype$X
write.csv(results,"severity_fisher.csv",quote = F)
###Fig6D
severity <- read.table("severity_3.csv",header = T,sep = ",",stringsAsFactors = FALSE)
severity <- subset(severity, severity[, 2] != 0 & severity[, 3] != 0)
severity$order <- severity$C2_ratio-severity$C1_ratio
#severity$order1 <- severity$C2_ratio_new-severity$C1_ratio_new
severity <- severity %>%
  arrange(order)
order <- severity$X
paint <- select(severity,-C1,-C2,-C1_all,-C2_all,-C1_ratio_new,-C2_ratio_new,-order)
paint1 <- melt(paint,id="X")
paint1$X <- factor(paint1$X,levels = order)
paint1$variable <- factor(paint1$variable,levels = c("C2_ratio","C1_ratio"))
ggplot(paint1, aes(x =value , y = X, fill = variable)) +
  geom_bar(stat = "identity", position = position_dodge(),size=0.5) +
  scale_fill_manual(values = c("C2_ratio"="#DF8F44FF","C1_ratio"="#374E55FF"))+
  scale_x_continuous(limits = c(0,1))+
  labs(x = "Severe Symptom Prevalence (%)", y = "") + 
  theme_bw()
### C1 vs C2
data <- matrix(c(2,8,10-2,10-8), nrow = 2)
rownames(data) <- c("C1", "C2")
colnames(data) <- c("more", "less")
fisher_test <- fisher.test(data)

###duration
duration <- read.table("duration_40_last.csv",header = T,sep = ",",stringsAsFactors = FALSE)
meta <- read.table("symptom_meta_27.csv",header = T,sep = ",",stringsAsFactors = FALSE)
meta <- select(meta,-people,-sex,-age,-bmi)
duration1 <- merge(meta,duration,by="ID")
duration2 <- select(duration1,-ID)
duration3 <- melt(duration2,id="cluster")
result <- duration2 %>%
  group_by(cluster) %>%
  summarise(across(1:19, list(median = median,min=min,max=max,mean=mean), na.rm = TRUE))
result <- as.data.frame(t(result))
write.csv(result,"duration_mean.csv",quote = F)
result <- duration2 %>%
  group_by(cluster) %>%  
  summarise(across(
    .cols = 1:19,        
    .fns = list(
      max = ~ max(., na.rm = TRUE),   
      min = ~ min(., na.rm = TRUE),    
      median = ~ median(., na.rm = TRUE),  
      mean = ~ mean(., na.rm = TRUE)  
    )
  ))
result <- as.data.frame(t(result))
write.csv(result,"duration_paint.csv",quote = F)
duration_paint <- read.table("duration_paint.csv",header = T,sep = ",",stringsAsFactors = FALSE)
library(tidyr)
duration_paint1 <-  duration_paint %>%
  extract(
    col = X,
    into = c("symptom", "type"),
    regex = "(.*)_([^_]+)$",
    remove = FALSE
  )
duration_paint2 <- select(duration_paint1,-X)
duration_paint3 <- duration_paint2 %>%
  pivot_longer(
    cols = c(C1, C2),
    names_to = "cluster",
    values_to = "value"
  )
duration_paint4 <- duration_paint3 %>%
  pivot_wider(
    names_from = type,
    values_from = value
  )
write.csv(duration_paint4,"duration_paint_new.csv",quote = F,row.names = F)
order <- duration_paint2[duration_paint2$type == "mean",]
order$order <- order$C2 - order$C1
order <- order %>%
  arrange(order)
order <- order$symptom
duration_paint4$symptom <- factor(duration_paint4$symptom,levels = order)
duration_paint4$cluster <- factor(duration_paint4$cluster,levels = c("C2","C1"))
###FigS8B
ggplot(duration_paint4, aes(y = symptom, x = mean )) +
  geom_errorbar(aes(xmin = min, xmax = max,color = cluster), size=1,position = position_dodge(width = 0.5), width = 0.3) +  # 添加误差线，避免重叠
  geom_point(aes(fill=cluster),position = position_dodge(width = 0.5),shape=21,size=3) +  
  scale_color_manual(values=c("C2"="#DF8F44FF","C1"="#374E55FF"))+
  scale_fill_manual(values=c("C2"="#DF8F44FF","C1"="#374E55FF"))+
  labs(y = "", x = "Average Symptom Duration", color = "cluster")+  
  theme_bw()
###lm
duration <- read.table("duration_40_last.csv",header = T,sep = ",",stringsAsFactors = FALSE)
meta <- read.table("symptom_meta_27.csv",header = T,sep = ",",stringsAsFactors = FALSE)
duration_new1 <- merge(meta,duration,by="ID")
duration_new1$cluster <- as.factor(duration_new1$cluster)
duration_new1$sex <- as.factor(duration_new1$sex) 
duration_new1$age <- as.numeric(duration_new1$age)
duration_new1$bmi <- as.numeric(duration_new1$bmi)
for (i in 7:ncol(duration_new1)) {
  duration_new1[,i] <- as.numeric(duration_new1[,i])
}
taxa <- colnames(duration_new1)[7:ncol(duration_new1)]
library(broom)
results <- lapply(taxa, function(taxa_col) {
  formula <- as.formula(paste( taxa_col, " ~ cluster + sex + age + bmi", sep = ""))
  model <- lm(formula, data = duration_new1)
  tidy(model) 
})
results_df <- do.call(rbind, results)
results_df$taxa <- rep(taxa, each = nrow(results_df) / length(taxa))
results_df$BH <- p.adjust(results_df$p.value,method = "BH")
write.csv(results_df,"duration_lm.csv",quote = F)
###C1 VS C2
data <- matrix(c(9,10,19-9,19-10), nrow = 2)
rownames(data) <- c("C1", "C2")
colnames(data) <- c("more", "less")
fisher_test <- fisher.test(data)


