###########################  loading environment ############################
library("dplyr")
library("tidyr")
library(Rmisc)
library(agricolae)
library("Hmisc")
library(vegan)
library(tidyverse)
library(data.table)
library(relaimpo)
library(igraph)
library("xlsx")
library(piecewiseSEM)
library(relaimpo)

format_percent <- function(input){
  out = input
  for (x in 1:length(input)){
    out[x] = input[x]/sum(input)
  }
  return(out)
}

anova_detail<- function(phage.hostspecies.percent){
  phage.0.otu = t(cbind(phage.hostspecies.percent[,1:4],phage.hostspecies.percent[,21:24]))
  phage.3.otu = t(cbind(phage.hostspecies.percent[,5:8],phage.hostspecies.percent[,25:28]))
  phage.4.otu = t(cbind(phage.hostspecies.percent[,9:12],phage.hostspecies.percent[,29:32]))
  phage.5.otu = t(cbind(phage.hostspecies.percent[,13:16],phage.hostspecies.percent[,33:36]))
  phage.6.otu = t(cbind(phage.hostspecies.percent[,17:20],phage.hostspecies.percent[,37:40]))
  
  phage.all = t(phage.hostspecies.percent)
  group.phage = c(rep("H",20),rep("D",20))
  input = cbind.data.frame(group.phage,phage.all)
  input$group.phage= as.factor(input$group.phage)
  anova_results = matrix(nrow = ncol(input)-1,ncol = 8)
  rownames(anova_results) = colnames(input)[2:ncol(input)]
  colnames(anova_results) = c("Df_factor_total.time","Df_residuals_total.time","SumSq_factor_total.time","SumSq_residuals_total.time","MeanSq_factor_total.time","MeanSq_residuals_total.time","p_value_total.time","F_value_total.time")
  for(i in seq(2,ncol(input))){
    anova.summary <- summary(aov(as.numeric(input[,i])~input[,1],data=input))
    anova_results[i-1,1] = anova.summary[[1]][1]$Df[1]
    anova_results[i-1,2] = anova.summary[[1]][1]$Df[2]
    anova_results[i-1,3] = anova.summary[[1]][2]$`Sum Sq`[1]
    anova_results[i-1,4] = anova.summary[[1]][2]$`Sum Sq`[2]
    anova_results[i-1,5] = anova.summary[[1]][3]$`Mean Sq`[1]
    anova_results[i-1,6] = anova.summary[[1]][3]$`Mean Sq`[2]
    anova_results[i-1,7] = anova.summary[[1]][5]$`Pr(>F)`[1]
    anova_results[i-1,8] = anova.summary[[1]][4]$`F value`[1]
  }
  
  anova_results_all = as.data.frame(anova_results)
  
  group.phage = c("H","H","H","H","D","D","D","D")
  
  input = cbind.data.frame(group.phage,phage.0.otu)
  input$group.phage= as.factor(input$group.phage)
  anova_results = matrix(nrow = ncol(input)-1,ncol = 8)
  rownames(anova_results) = colnames(input)[2:ncol(input)]
  colnames(anova_results) = c("Df_factor_0","Df_residuals_0","SumSq_factor_0","SumSq_residuals_0","MeanSq_factor_0","MeanSq_residuals_0","p_value_0","F_value_0")
  for(i in seq(2,ncol(input))){
    anova.summary <- summary(aov(as.numeric(input[,i])~input[,1],data=input))
    anova_results[i-1,1] = anova.summary[[1]][1]$Df[1]
    anova_results[i-1,2] = anova.summary[[1]][1]$Df[2]
    anova_results[i-1,3] = anova.summary[[1]][2]$`Sum Sq`[1]
    anova_results[i-1,4] = anova.summary[[1]][2]$`Sum Sq`[2]
    anova_results[i-1,5] = anova.summary[[1]][3]$`Mean Sq`[1]
    anova_results[i-1,6] = anova.summary[[1]][3]$`Mean Sq`[2]
    anova_results[i-1,7] = anova.summary[[1]][5]$`Pr(>F)`[1]
    anova_results[i-1,8] = anova.summary[[1]][4]$`F value`[1]
  }
  
  anova_results0 = as.data.frame(anova_results)
  
  input = cbind.data.frame(group.phage,phage.3.otu)
  input$group.phage= as.factor(input$group.phage)
  anova_results = matrix(nrow = ncol(input)-1,ncol = 8)
  rownames(anova_results) = colnames(input)[2:ncol(input)]
  colnames(anova_results) = c("Df_factor_3","Df_residuals_3","SumSq_factor_3","SumSq_residuals_3","MeanSq_factor_3","MeanSq_residuals_3","p_value_3","F_value_3")
  for(i in seq(2,ncol(input))){
    anova.summary <- summary(aov(as.numeric(input[,i])~input[,1],data=input))
    anova_results[i-1,1] = anova.summary[[1]][1]$Df[1]
    anova_results[i-1,2] = anova.summary[[1]][1]$Df[2]
    anova_results[i-1,3] = anova.summary[[1]][2]$`Sum Sq`[1]
    anova_results[i-1,4] = anova.summary[[1]][2]$`Sum Sq`[2]
    anova_results[i-1,5] = anova.summary[[1]][3]$`Mean Sq`[1]
    anova_results[i-1,6] = anova.summary[[1]][3]$`Mean Sq`[2]
    anova_results[i-1,7] = anova.summary[[1]][5]$`Pr(>F)`[1]
    anova_results[i-1,8] = anova.summary[[1]][4]$`F value`[1]
  }
  
  anova_results3 = as.data.frame(anova_results)
  
  input = cbind.data.frame(group.phage,phage.4.otu)
  input$group.phage= as.factor(input$group.phage)
  anova_results = matrix(nrow = ncol(input)-1,ncol = 8)
  rownames(anova_results) = colnames(input)[2:ncol(input)]
  colnames(anova_results) = c("Df_factor_4","Df_residuals_4","SumSq_factor_4","SumSq_residuals_4","MeanSq_factor_4","MeanSq_residuals_4","p_value_4","F_value_4")
  for(i in seq(2,ncol(input))){
    anova.summary <- summary(aov(as.numeric(input[,i])~input[,1],data=input))
    anova_results[i-1,1] = anova.summary[[1]][1]$Df[1]
    anova_results[i-1,2] = anova.summary[[1]][1]$Df[2]
    anova_results[i-1,3] = anova.summary[[1]][2]$`Sum Sq`[1]
    anova_results[i-1,4] = anova.summary[[1]][2]$`Sum Sq`[2]
    anova_results[i-1,5] = anova.summary[[1]][3]$`Mean Sq`[1]
    anova_results[i-1,6] = anova.summary[[1]][3]$`Mean Sq`[2]
    anova_results[i-1,7] = anova.summary[[1]][5]$`Pr(>F)`[1]
    anova_results[i-1,8] = anova.summary[[1]][4]$`F value`[1]
  }
  
  anova_results4 = as.data.frame(anova_results)
  
  input = cbind.data.frame(group.phage,phage.5.otu)
  input$group.phage= as.factor(input$group.phage)
  anova_results = matrix(nrow = ncol(input)-1,ncol = 8)
  rownames(anova_results) = colnames(input)[2:ncol(input)]
  colnames(anova_results) = c("Df_factor_5","Df_residuals_5","SumSq_factor_5","SumSq_residuals_5","MeanSq_factor_5","MeanSq_residuals_5","p_value_5","F_value_5")
  for(i in seq(2,ncol(input))){
    anova.summary <- summary(aov(as.numeric(input[,i])~input[,1],data=input))
    anova_results[i-1,1] = anova.summary[[1]][1]$Df[1]
    anova_results[i-1,2] = anova.summary[[1]][1]$Df[2]
    anova_results[i-1,3] = anova.summary[[1]][2]$`Sum Sq`[1]
    anova_results[i-1,4] = anova.summary[[1]][2]$`Sum Sq`[2]
    anova_results[i-1,5] = anova.summary[[1]][3]$`Mean Sq`[1]
    anova_results[i-1,6] = anova.summary[[1]][3]$`Mean Sq`[2]
    anova_results[i-1,7] = anova.summary[[1]][5]$`Pr(>F)`[1]
    anova_results[i-1,8] = anova.summary[[1]][4]$`F value`[1]
  }
  
  anova_results5 = as.data.frame(anova_results)
  
  input = cbind.data.frame(group.phage,phage.6.otu)
  input$group.phage= as.factor(input$group.phage)
  anova_results = matrix(nrow = ncol(input)-1,ncol = 8)
  rownames(anova_results) = colnames(input)[2:ncol(input)]
  colnames(anova_results) = c("Df_factor_6","Df_residuals_6","SumSq_factor_6","SumSq_residuals_6","MeanSq_factor_6","MeanSq_residuals_6","p_value_6","F_value_6")
  for(i in seq(2,ncol(input))){
    anova.summary <- summary(aov(as.numeric(input[,i])~input[,1],data=input))
    anova_results[i-1,1] = anova.summary[[1]][1]$Df[1]
    anova_results[i-1,2] = anova.summary[[1]][1]$Df[2]
    anova_results[i-1,3] = anova.summary[[1]][2]$`Sum Sq`[1]
    anova_results[i-1,4] = anova.summary[[1]][2]$`Sum Sq`[2]
    anova_results[i-1,5] = anova.summary[[1]][3]$`Mean Sq`[1]
    anova_results[i-1,6] = anova.summary[[1]][3]$`Mean Sq`[2]
    anova_results[i-1,7] = anova.summary[[1]][5]$`Pr(>F)`[1]
    anova_results[i-1,8] = anova.summary[[1]][4]$`F value`[1]
  }
  
  anova_results6 = as.data.frame(anova_results)
  
  anova_results_all = cbind(anova_results_all,anova_results0,anova_results3,anova_results4,anova_results5,anova_results6)
  return(anova_results_all)
}



dataframe2matrix  <- function(input){
  dataframe = input
  rownames(dataframe) = dataframe[,1]
  dataframe = dataframe[,-1]
  return(dataframe)
}


rm_anova <- function(data.phage2){
  data.phage2[,2] = as.factor(data.phage2[,2])
  data.phage2[,3] = as.factor(data.phage2[,3])
  data.phage2[,4] = as.factor(data.phage2[,4])
  result = c(NA,NA,NA,NA,NA,NA,NA,NA)
  for ( i in 5:ncol(data.phage2)){
    fit<-aov(data.phage2[,i]~time*disease_outcome+Error(plant/time),data=data.phage2)
    out = tidy(fit)
    name = rep(colnames(data.phage2)[i],nrow(out))
    out2 = cbind.data.frame(name,out)
    result = rbind.data.frame(result,out2)
  }
  result = result[-1,]
}

##########################  Repeated measures ANOVA ############################
data = read.table(""rmanova.data.txt,sep="\t",header = T)
data$time = as.factor(data$time)
data$disease_outcome = as.factor(data$disease_outcome)
data$plant = as.factor(data$plant)

rmanova.result = rm_anova(data)


##########################  one-way ANOVA ############################

data = read.table("anova.data.txt",sep="\t",header = T)
anova=anova_detail(dataframe2matrix(data))

###########################  shannon  ###################
otu.data = read.table("vOTU.data.txt",sep="\t",header = T)

otu.t = t(dataframe2matrix(otu.data))
Shannon=vegan::diversity(otu.t,index='shannon')
alpha = data.frame(sample = rownames(otu.t),Shannon)

###########################  lm  ###################

data = read.table("rmanova.data.txt",sep="\t",head=T,row.names=1,check.names = F,quote = "")

shannon.h = data[5:20,c(2,3,5,6)]
shannon.d = data[25:40,c(2,3,5,6)]
summary(lm(bacteria.Shannon~vOTU.Shannon,shannon.h))
summary(lm(bacteria.Shannon~vOTU.Shannon,shannon.d))

###########################  phage PCA ###################
phage.otu = read.table("vOTU.data.txt",sep="\t",head=T,row.names=1)
phage.0.otu = t(cbind(phage.otu[,1:4],phage.otu[,21:24]))
phage.3.otu = t(cbind(phage.otu[,5:8],phage.otu[,25:28]))
phage.4.otu = t(cbind(phage.otu[,9:12],phage.otu[,29:32]))
phage.5.otu = t(cbind(phage.otu[,13:16],phage.otu[,33:36]))
phage.6.otu = t(cbind(phage.otu[,17:20],phage.otu[,37:40]))
group.phage = as.factor(c("H","H","H","H","D","D","D","D"))
phage.0.pca <- summary(prcomp(phage.0.otu))
phage.3.pca <- summary(prcomp(phage.3.otu))
phage.4.pca <- summary(prcomp(phage.4.otu))
phage.5.pca <- summary(prcomp(phage.5.otu))
phage.6.pca <- summary(prcomp(phage.6.otu))

phage.pca = data.frame(group=rep(c("H","H","H","H","D","D","D","D"),5),rbind.data.frame(phage.0.pca$x,phage.3.pca$x,phage.4.pca$x,phage.5.pca$x,phage.6.pca$x))

group.0.phage = data.frame(group.phage)
rownames(group.0.phage) = rownames(phage.0.otu)
adonis.0.result = adonis2(phage.0.otu~group.phage,data = group.0.phage,permutations = 999,method="bray")
group.3.phage = data.frame(group.phage)
rownames(group.3.phage) = rownames(phage.3.otu)
adonis.3.result = adonis2(phage.3.otu~group.phage,data = group.3.phage,permutations = 999,method="bray")
group.4.phage = data.frame(group.phage)
rownames(group.4.phage) = rownames(phage.4.otu)
adonis.4.result = adonis2(phage.4.otu~group.phage,data = group.4.phage,permutations = 999,method="bray")
group.5.phage = data.frame(group.phage)
rownames(group.5.phage) = rownames(phage.5.otu)
adonis.5.result = adonis2(phage.5.otu~group.phage,data = group.5.phage,permutations = 999,method="bray")
group.6.phage = data.frame(group.phage)
rownames(group.6.phage) = rownames(phage.6.otu)
adonis.6.result = adonis2(phage.6.otu~group.phage,data = group.6.phage,permutations = 999,method="bray")
adonis.6.result

###########################  distance  ###################

phage.otu = read.table("vOTU.data.txt",sep="\t",head=T,row.names=1)
phage.0.otu = t(cbind(phage.otu[,1:4],phage.otu[,21:24]))
phage.3.otu = t(cbind(phage.otu[,5:8],phage.otu[,25:28]))
phage.4.otu = t(cbind(phage.otu[,9:12],phage.otu[,29:32]))
phage.5.otu = t(cbind(phage.otu[,13:16],phage.otu[,33:36]))
phage.6.otu = t(cbind(phage.otu[,17:20],phage.otu[,37:40]))
group.phage = as.factor(c("H","H","H","H","D","D","D","D"))

group.0.phage = data.frame(group.phage)
rownames(group.0.phage) = rownames(phage.0.otu)
distance.0.bray<-vegdist(phage.0.otu,method = 'bray')
group.3.phage = data.frame(group.phage)
rownames(group.3.phage) = rownames(phage.3.otu)
distance.3.bray<-vegdist(phage.3.otu,method = 'bray')
group.4.phage = data.frame(group.phage)
rownames(group.4.phage) = rownames(phage.4.otu)
distance.4.bray<-vegdist(phage.4.otu,method = 'bray')
group.5.phage = data.frame(group.phage)
rownames(group.5.phage) = rownames(phage.5.otu)
distance.5.bray<-vegdist(phage.5.otu,method = 'bray')
group.6.phage = data.frame(group.phage)
rownames(group.6.phage) = rownames(phage.6.otu)
distance.6.bray<-vegdist(phage.6.otu,method = 'bray')

distance.0.bray = as.matrix(distance.0.bray)
group.distance.0 = cbind(as.numeric(distance.0.bray[5:8,1:4]),rep(0,16),c(rep("H-2",4),rep("H-3",4),rep("H-4",4),rep("H-5",4)),rep(c("D-1","D-2","D-3","D-5"),4))
distance.3.bray = as.matrix(distance.3.bray)
group.distance.3 = cbind(as.numeric(distance.3.bray[5:8,1:4]),rep(3,16),c(rep("H-2",4),rep("H-3",4),rep("H-4",4),rep("H-5",4)),rep(c("D-1","D-2","D-3","D-5"),4))
distance.4.bray = as.matrix(distance.4.bray)
group.distance.4 = cbind(as.numeric(distance.4.bray[5:8,1:4]),rep(4,16),c(rep("H-2",4),rep("H-3",4),rep("H-4",4),rep("H-5",4)),rep(c("D-1","D-2","D-3","D-5"),4))
distance.5.bray = as.matrix(distance.5.bray)
group.distance.5 = cbind(as.numeric(distance.5.bray[5:8,1:4]),rep(5,16),c(rep("H-2",4),rep("H-3",4),rep("H-4",4),rep("H-5",4)),rep(c("D-1","D-2","D-3","D-5"),4))
distance.6.bray = as.matrix(distance.6.bray)
group.distance.6 = cbind(as.numeric(distance.6.bray[5:8,1:4]),rep(6,16),c(rep("H-2",4),rep("H-3",4),rep("H-4",4),rep("H-5",4)),rep(c("D-1","D-2","D-3","D-5"),4))
group.distance = as.data.frame(rbind(group.distance.0,group.distance.3,group.distance.4,group.distance.5,group.distance.6))
colnames(group.distance) = c("distance","time","plant1","plant2")

group.distance$time = as.factor(group.distance$time)
group.distance$distance = as.numeric(group.distance$distance)
result1 <- aov(distance~time,data=group.distance)
hsd <-HSD.test(result1,"time",alpha = 0.05)
summary <- summarySE(data=group.distance, measurevar="distance",groupvars="time")
arrange = data.frame(hsd$groups)
arrange2 = cbind(row.names(arrange),arrange)
colnames(arrange2)[1] = "id"
arrange2 = arrange2[order(arrange2[,"id"],decreasing=F),]
summary2=cbind(summary,arrange2)
summary2 = summary2[,-3]

group.distance$time = as.numeric(as.character(group.distance$time))

model = lm(distance~time,subset(group.distance,time != 0))
summary(model)



###########################  SEM  ###################
data.effect.h = read.table("effect.H.txt",sep="\t",header = T)
data.effect.d = read.table("effect.D.txt",sep="\t",header = T)

facilitator_bacteria.h = subset(data.effect.h,group=="facilitator_bacteria")[,-c(2,3,4)]
facilitator_bacteria.d = subset(data.effect.d,group=="facilitator_bacteria")[,-c(2,3,4)]
facilitator_associated_phage.h = subset(data.effect.h,group=="facilitator_associated_phage")[,-c(2,3,4)]
facilitator_associated_phage.d = subset(data.effect.d,group=="facilitator_associated_phage")[,-c(2,3,4)]
inhibitor_bacteria.h = subset(data.effect.h,group=="inhibitor_bacteria")[,-c(2,3,4)]
inhibitor_bacteria.d = subset(data.effect.d,group=="inhibitor_bacteria")[,-c(2,3,4)]
inhibitor_associated_phage.h = subset(data.effect.h,group=="inhibitor_associated_phage")[,-c(2,3,4)]
inhibitor_associated_phage.d = subset(data.effect.d,group=="inhibitor_associated_phage")[,-c(2,3,4)]
Pathogen_phage.h = subset(data.effect.h,group=="Pathogen_phage")[,-c(2,3,4)]
Pathogen_phage.d = subset(data.effect.d,group=="Pathogen_phage")[,-c(2,3,4)]
Pathogen.h = subset(data.effect.h,group=="Pathogen")[,-c(2,3,4)]
Pathogen.d = subset(data.effect.d,group=="Pathogen")[,-c(2,3,4)]

facilitator_bacteria.h.shannon=vegan::diversity(t(dataframe2matrix(facilitator_bacteria.h)),index='shannon')
facilitator_bacteria.d.shannon=vegan::diversity(t(dataframe2matrix(facilitator_bacteria.d)),index='shannon')
facilitator_associated_phage.h.shannon=vegan::diversity(t(dataframe2matrix(facilitator_associated_phage.h)),index='shannon')
facilitator_associated_phage.d.shannon=vegan::diversity(t(dataframe2matrix(facilitator_associated_phage.d)),index='shannon')
inhibitor_bacteria.h.shannon=vegan::diversity(t(dataframe2matrix(inhibitor_bacteria.h)),index='shannon')
inhibitor_bacteria.d.shannon=vegan::diversity(t(dataframe2matrix(inhibitor_bacteria.d)),index='shannon')
inhibitor_associated_phage.h.shannon=vegan::diversity(t(dataframe2matrix(inhibitor_associated_phage.h)),index='shannon')
inhibitor_associated_phage.d.shannon=vegan::diversity(t(dataframe2matrix(inhibitor_associated_phage.d)),index='shannon')
Pathogen_phage.h.shannon=vegan::diversity(t(dataframe2matrix(Pathogen_phage.h)),index='shannon')
Pathogen_phage.d.shannon=vegan::diversity(t(dataframe2matrix(Pathogen_phage.d)),index='shannon')

data.effect.h.shannon = cbind.data.frame(
  data.frame(time=c(rep(0,4),rep(3,4),rep(4,4),rep(5,4),rep(6,4)),
             plant=as.factor(rep(c(2,3,4,5),5))),
  facilitator_associated_phage.h.shannon,
  facilitator_bacteria.h.shannon,
  inhibitor_associated_phage.h.shannon,
  inhibitor_bacteria.h.shannon,
  Pathogen_phage.h.shannon,
  t(dataframe2matrix(Pathogen.h)))
data.effect.d.shannon = cbind.data.frame(
  data.frame(time=c(rep(0,4),rep(3,4),rep(4,4),rep(5,4),rep(6,4)),
             plant=as.factor(rep(c(2,3,4,5),5))),
  facilitator_associated_phage.d.shannon,
  facilitator_bacteria.d.shannon,
  inhibitor_associated_phage.d.shannon,
  inhibitor_bacteria.d.shannon,
  Pathogen_phage.d.shannon,
  t(dataframe2matrix(Pathogen.d)))

colnames(data.effect.h.shannon)[8]="rs"
colnames(data.effect.d.shannon)[8]="rs"

h_psem <- psem(
  
  lm(inhibitor_bacteria.h.shannon ~ inhibitor_associated_phage.h.shannon, data = data.effect.h.shannon),
  lm(rs ~ Pathogen_phage.h.shannon  + inhibitor_bacteria.h.shannon + inhibitor_associated_phage.h.shannon, data = data.effect.h.shannon),
  data = data.effect.h.shannon)
fisherC(h_psem)
AIC(h_psem)
summary(h_psem, .progressBar = T)


d_psem <- psem(
  
  lm(inhibitor_bacteria.d.shannon ~ inhibitor_associated_phage.d.shannon, data = data.effect.d.shannon),
  lm(rs ~ Pathogen_phage.d.shannon    +inhibitor_bacteria.d.shannon +inhibitor_associated_phage.d.shannon , data = data.effect.d.shannon),
  data = data.effect.d.shannon)
fisherC(d_psem)
AIC(d_psem)
summary(d_psem, .progressBar = FALSE)

lm.data.h = lm(rs ~ Pathogen_phage.h.shannon  + inhibitor_bacteria.h.shannon + inhibitor_associated_phage.h.shannon,data = data.effect.h.shannon)
calc.relimp(lm.data.h, type = c("lmg"), rela = F)

lm.data.d = lm(rs ~ Pathogen_phage.d.shannon  + inhibitor_bacteria.d.shannon + inhibitor_associated_phage.d.shannon,data = data.effect.d.shannon)
calc.relimp(lm.data.d, type = c("lmg"), rela = F)

###########################  net-stat  ###################

data.h= read.table("net.H.txt", sep = "\t",head = T,quote = "",as.is = T,check.names = T)
network.h = as.matrix(data.h[,1:2])
graph.h = graph_from_edgelist(network.h,directed = F)

edges.h = length(E(graph.h))
nodes.h = length(V(graph.h))
connectance.h = edge_density(graph.h,loops=FALSE)
edge.connectivity.h = edge_connectivity(graph.h)
centralization.betweenness.h = centralization.betweenness(graph.h)$centralization 
centralization.degree.h = centralization.degree(graph.h)$centralization

data.d= read.table("net.D.txt", sep = "\t",head = T,quote = "",as.is = T,check.names = T)
network.d = as.matrix(data.d[,1:2])
graph.d = graph_from_edgelist(network.d,directed = F)

edges.d = length(E(graph.d))
nodes.d = length(V(graph.d))
connectance.d = edge_density(graph.d,loops=FALSE)
edge.connectivity.d = edge_connectivity(graph.d)
centralization.betweenness.d = centralization.betweenness(graph.d)$centralization 
centralization.degree.d = centralization.degree(graph.d)$centralization

stat.network = data.frame(stat = c("nodes","edges","connectance","connectivity","betweenness centralization","Degree centralization"),
                          h = c(as.numeric(length(E(graph.h))),as.numeric(length(V(graph.h))),as.numeric(edge_density(graph.h,loops=FALSE)),as.numeric(edge_connectivity(graph.h)),as.numeric(centralization.betweenness(graph.h)$centralization) ,as.numeric(centralization.degree(graph.h)$centralization)),
                          d = c(as.numeric(length(E(graph.d))),as.numeric(length(V(graph.d))),as.numeric(edge_density(graph.d,loops=FALSE)),as.numeric(edge_connectivity(graph.d)),as.numeric(centralization.betweenness(graph.d)$centralization) ,as.numeric(centralization.degree(graph.d)$centralization)))





