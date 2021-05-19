# required libraries
library(edgeR)
library(tidyverse)
library(readr)
library(readxl)
library(cowplot) # for plotting both genotypes or density

# load reads mapped to Brassica genome (either v1.5 annotation or v3.0 annotation)
getwd()
# ## for exp1 v1.5annotation
# counts.exp1.v1.5 <- readr::read_csv(file.path("..","v1.5annotation","20170617-samples","input","raw_counts.csv.gz"),col_names=TRUE)
# counts.exp1.v1.5 # make sure this is v1.5 annotation (look target_id column)
# ## for exp1 v3.0annotation
# counts.exp1.v3.0 <- readr::read_csv(file.path("..","v3.0annotation","20170617-samples","input","20170617_V3.0_raw_counts.csv.gz"),col_names=TRUE)
# ### cpm
# cpm.exp1.leaf.v3.0 <- readr::read_csv(file.path("..","v3.0annotation","20170617-samples","output","cpm_wide_20170617_leaf_samples.csv.gz"),col_names=TRUE)
# cpm.exp1.root.v3.0 <- readr::read_csv(file.path("..","v3.0annotation","20170617-samples","output","cpm_wide_20170617_root_samples.csv.gz"),col_names=TRUE)
# 
# counts.exp1.v3.0 # make sure this is v3.0 annotation (look target_id column)
# ## for exp3 v3.0annotation
# counts.exp3.v3.0 <- readr::read_csv(file.path("..","v3.0annotation","20180202-samples","input","20180202_V3.0_raw_counts.csv.gz"),col_names=TRUE)
# counts.exp3.v3.0
# ### cpm
# cpm.exp3.leaf.v3.0<-readr::read_csv(file.path("..","v3.0annotation","20180202-samples","output","cpm_wide_20180202_leaf_samples.csv.gz"),col_names=TRUE)
# cpm.exp3.root.v3.0<-readr::read_csv(file.path("..","v3.0annotation","20180202-samples","output","cpm_wide_20180202_root_samples.csv.gz"),col_names=TRUE)
# 
## for timecourse v3.0annotation
# loading DEG object created in "02_Normalize_DGE.Rmd"
load("../output/timecourseDGE.Rdata")
cpm.timecourse.v3.0 <- cpm(dge) %>% as_tibble() %>% bind_cols(data.frame(transcript_ID=rownames(dge$counts)),.)
cpm.timecourse.v3.0.log <- cpm(dge,log=TRUE) %>% as_tibble() %>% bind_cols(data.frame(transcript_ID=rownames(dge$counts)),.)
# sample files
# # exp1 (20170617-samples)
# sample.description.exp1<-readr::read_csv(file.path("..","v1.5annotation","20170617-samples","output","Br.mbio.e1.sample.description.csv"))
# # exp3 (20180202-samples)
# sample.description.exp3<-readr::read_csv(file.path("..","v3.0annotation","20180202-samples","output","Br.mbio.e3.sample.description.csv"))
# timecourse
sample.description.timecourse<-dge$samples %>% bind_cols(data.frame(sample=rownames(dge$samples)),.)
  
# functions for drawing expression pattern
## for exp1 (v1.5 annotation)
expression.pattern.Br.graph.exp1.v1.5annotation<-function(data=counts.exp1.v1.5,target.genes,sample.description=sample.description.exp1,title="",geno){
  data[is.na(data)] <- 0 #
  # select genes and add sample info
  data.temp<-data %>% filter(target_id %in% target.genes) %>% gather(sample,value,-target_id) %>%
    inner_join(sample.description, by="sample") 
  # 
  if(geno=="both") { # needs to impove this
    # ggplot(data.temp, aes(x=genotype,y=value))  + geom_jitter(alpha = 0.5,aes(colour=trt,shape=tissue) )  + theme_bw() + facet_grid(target_id~tissue,scales="free") + theme(strip.text.y=element_text(angle=0),axis.text.x=element_text(angle=90)) + theme(legend.position="bottom") + labs(title=title)
    p.FPsc<-data.temp %>% filter(genotype=="FPsc") %>% ggplot(aes(x=trt,y=value))  + geom_jitter(alpha = 0.5,aes(colour=trt,shape=tissue) )  + theme_bw() + facet_grid(target_id~tissue,scales="free") + theme(strip.text.y=element_text(angle=0),axis.text.x=element_text(angle=90)) + theme(legend.position="bottom") + labs(title=" ")
    p.R500<-data.temp %>% filter(genotype=="R500") %>% ggplot(aes(x=trt,y=value))  + geom_jitter(alpha = 0.5,aes(colour=trt,shape=tissue) )  + theme_bw() + facet_grid(target_id~tissue,scales="free") + theme(strip.text.y=element_text(angle=0),axis.text.x=element_text(angle=90)) + theme(legend.position="bottom") + labs(title=" ")
    # merge two plots
    plot_grid(p.FPsc,p.R500,labels=c(paste(title,"FPsc"),paste(title,"R500")))
  } else if(geno=="FPsc"|geno=="R500") {
    data.temp %>% filter(genotype==geno) %>% ggplot(aes(x=trt,y=value))  + geom_jitter(alpha = 0.5,aes(colour=trt,shape=tissue) )  + theme_bw() + facet_grid(target_id~tissue,scales="free") + theme(strip.text.y=element_text(angle=0),axis.text.x=element_text(angle=90)) + theme(legend.position="bottom") + labs(title=title)
  }  else {print("Specify genotype.");stop}
}
# for exp3 (v3.0 annotation)
expression.pattern.Br.graph.exp3<-function(data=counts.exp3.v3.0,target.genes,sample.description=sample.description.exp3,title="",dens="both"){
  data[is.na(data)] <- 0 #
  # select genes and add sample info
  data.temp<-data %>% filter(target_id %in% target.genes) %>% gather(sample,value,-target_id) %>%
    inner_join(sample.description, by="sample") 
  # plot (separated by density info)
  if(dens=="both") { # needs to improve using cowplot
    # ggplot(data.temp, aes(x=density,y=value))  + geom_jitter(alpha = 0.5,aes(colour=trt,shape=tissue) )  + theme_bw() + facet_grid(target_id~tissue,scales="free") + theme(strip.text.y=element_text(angle=0),axis.text.x=element_text(angle=90)) + theme(legend.position="bottom") + labs(title=title)
    p.cr<-data.temp %>% filter(density=="cr") %>% ggplot(aes(x=trt,y=value))  + geom_jitter(alpha = 0.5,aes(colour=trt,shape=tissue) )  + theme_bw() + facet_grid(target_id~tissue,scales="free") + theme(strip.text.y=element_text(angle=0),axis.text.x=element_text(angle=90)) + theme(legend.position="bottom") + labs(title=" ")
    p.un<-data.temp %>% filter(density=="un") %>% ggplot(aes(x=trt,y=value))  + geom_jitter(alpha = 0.5,aes(colour=trt,shape=tissue) )  + theme_bw() + facet_grid(target_id~tissue,scales="free") + theme(strip.text.y=element_text(angle=0),axis.text.x=element_text(angle=90)) + theme(legend.position="bottom") + labs(title=" ")
    # merge two plots
    plot_grid(p.cr,p.un,labels=c(paste(title,"cr"),paste(title,"un")))
  } else if(dens=="cr"|dens=="un") {
    data.temp %>% filter(density==dens) %>% ggplot(aes(x=trt,y=value))  + geom_jitter(alpha = 0.5,aes(colour=trt,shape=tissue) )  + theme_bw() + facet_grid(target_id~tissue,scales="free") + theme(strip.text.y=element_text(angle=0),axis.text.x=element_text(angle=90)) + theme(legend.position="bottom") + labs(title=title)
  }  else {print("Specify genotype.");stop}
}
## for exp1 (v3.0 annotation)
expression.pattern.Br.graph.exp1.v3.0annotation<-function(data=counts.exp1.v3.0,target.genes,sample.description=sample.description.exp1,title="",geno,tissue.type="root"){
  data[is.na(data)] <- 0 #
  # select genes and add sample info
  data.temp<-data %>% filter(target_id %in% target.genes) %>% gather(sample,value,-target_id) %>%
    inner_join(sample.description, by="sample") 
  # 
  if(geno=="both") { # needs to impove this
    # ggplot(data.temp, aes(x=genotype,y=value))  + geom_jitter(alpha = 0.5,aes(colour=trt,shape=tissue) )  + theme_bw() + facet_grid(target_id~tissue,scales="free") + theme(strip.text.y=element_text(angle=0),axis.text.x=element_text(angle=90)) + theme(legend.position="bottom") + labs(title=title)
    p.FPsc<-data.temp %>% filter(genotype=="FPsc",tissue==tissue.type) %>% ggplot(aes(x=trt,y=value))  + geom_jitter(alpha = 0.5,aes(colour=trt),width=0.2)  + theme_bw() + facet_grid(target_id~tissue,scales="free") + theme(strip.text.y=element_text(angle=0),axis.text.x=element_text(angle=90)) + theme(legend.position="bottom") + labs(title=" ")
    p.R500<-data.temp %>% filter(genotype=="R500",tissue==tissue.type) %>% ggplot(aes(x=trt,y=value))  + geom_jitter(alpha = 0.5,aes(colour=trt),width=0.2 )  + theme_bw() + facet_grid(target_id~tissue,scales="free") + theme(strip.text.y=element_text(angle=0),axis.text.x=element_text(angle=90)) + theme(legend.position="bottom") + labs(title=" ")
    # merge two plots
    plot_grid(p.FPsc,p.R500,labels=c(paste(title,"FPsc"),paste(title,"R500")))
  } else if(geno=="FPsc"|geno=="R500") {
    data.temp %>% filter(genotype==geno) %>% ggplot(aes(x=trt,y=value))  + geom_jitter(alpha = 0.5,aes(colour=trt,shape=tissue) )  + theme_bw() + facet_grid(target_id~tissue,scales="free") + theme(strip.text.y=element_text(angle=0),axis.text.x=element_text(angle=90)) + theme(legend.position="bottom") + labs(title=title)
  }  else {print("Specify genotype.");stop}
}

# only one tissue with the same y-axis for both genotypes
expression.pattern.Br.graph.exp1.v3.0annotation.2<-function(data=counts.exp1.v3.0,target.genes,sample.description=sample.description.exp1,title="",geno,tissue.type="root"){
  data[is.na(data)] <- 0 #
  # select genes and add sample info
  data.temp<-data %>% filter(target_id %in% target.genes) %>% gather(sample,value,-target_id) %>%
    inner_join(sample.description, by="sample") %>% filter(tissue==tissue.type)
  # 
  if(geno=="both") { # needs to impove this
    # ggplot(data.temp, aes(x=genotype,y=value))  + geom_jitter(alpha = 0.5,aes(colour=trt,shape=tissue) )  + theme_bw() + facet_grid(target_id~tissue,scales="free") + theme(strip.text.y=element_text(angle=0),axis.text.x=element_text(angle=90)) + theme(legend.position="bottom") + labs(title=title)
    p<-data.temp %>% ggplot(aes(x=trt,y=value))  + geom_jitter(alpha = 0.5,aes(colour=trt),width=0.2)  + theme_bw() + facet_grid(target_id~genotype,scales="free") + theme(strip.text.y=element_text(angle=0),axis.text.x=element_text(angle=90)) + theme(legend.position="bottom") + labs(title=" ")
    p
  } else if(geno=="FPsc"|geno=="R500") {
    data.temp %>% filter(genotype==geno) %>% ggplot(aes(x=trt,y=value))  + geom_jitter(alpha = 0.5,aes(colour=trt) )  + theme_bw() + facet_grid(target_id~.,scales="free") + theme(strip.text.y=element_text(angle=0),axis.text.x=element_text(angle=90)) + theme(legend.position="bottom") + labs(title=title)
  }  else {print("Specify genotype.");stop}
}
# only one tissue with the same y-axis for both genotypes with FDR info. (target.genes.FDR is a dataframe with gene name and FDR)
## cf. https://www.datacamp.com/community/tutorials/facets-ggplot-r by using labeller
## cf. 
expression.pattern.Br.graph.exp1.v3.0annotation.FDR.3<-function(data=counts.exp1.v3.0,target.genes.FDR,sample.description=sample.description.exp1,title="",geno,tissue.type="root"){
  data[is.na(data)] <- 0 #
  # select genes and add sample info
  data.temp<-data %>% 
    filter(target_id %in% target.genes.FDR$genes) %>% 
    gather(sample,value,-target_id) %>%
    inner_join(sample.description, by="sample") %>% 
    filter(tissue==tissue.type) %>%
    right_join(target.genes.FDR, by=c("target_id"="genes")) %>%
    mutate(FDR=format(FDR,digits=2,scientific=TRUE)) %>%
    unite(target_id.FDR,target_id,FDR,sep="\n FDR ")
  # for adding FDR,,, but under construction (070419). Making a new data frame for writing FDR on each plot. Not used anynore
   # data.temp2 <- data %>% 
   #  gather(sample,value,-target_id) %>%
   #  inner_join(sample.description, by="sample") %>% 
   #  filter(tissue==tissue.type) %>% 
   #  group_by(target_id,genotype) %>% 
   #  summarise(y.max.plus=max(value)*1.05) %>% 
   #  mutate(x=3) %>% 
   #  right_join(target.genes.FDR, by=c("target_id"="genes")) %>% 
   #  mutate(FDR=format(FDR,digits=2,scientific=TRUE))
  # %>% slice(1) %>% select(target_id,FDR,tissue,genotype) %>% 
####
    if(geno=="both") { # needs to impove this
    # ggplot(data.temp, aes(x=genotype,y=value))  + geom_jitter(alpha = 0.5,aes(colour=trt,shape=tissue) )  + theme_bw() + facet_grid(target_id~tissue,scales="free") + theme(strip.text.y=element_text(angle=0),axis.text.x=element_text(angle=90)) + theme(legend.position="bottom") + labs(title=title)
    #data.temp <- data %>% right_join(target.genes.FDR, by=c("target_id"="genes")) %>% gather(sample,value,-target_id,-FDR) %>%
    #  inner_join(sample.description, by="sample") %>% filter(tissue==tissue.type) 
    #p<-data.temp %>% ggplot(aes(x=trt,y=value))  + geom_jitter(alpha = 0.5,aes(colour=trt),width=0.2)  + theme_bw() + geom_text(aes(label=FDR))+facet_grid(target_id~genotype,scales="free") + theme(strip.text.y=element_text(angle=0),axis.text.x=element_text(angle=90)) + theme(legend.position="bottom") + labs(title=" ")
    p<-data.temp %>% ggplot(aes(x=trt,y=value))  + 
      geom_jitter(alpha = 0.5,aes(colour=trt,shape=as.character(block)),width=0.2)  + 
      theme_bw() +
      facet_grid(target_id.FDR~genotype,scales="free") + 
      theme(strip.text.y=element_text(angle=0),axis.text.x=element_text(angle=90)) +
      theme(legend.position="bottom") + labs(title=" ")
    #p <- p + geom_text(data=data.temp2,aes(x=x,y=y.max.plus,label=FDR))
    p
  } else if(geno=="FPsc"|geno=="R500") {
    data.temp %>% filter(genotype==geno) %>% ggplot(aes(x=trt,y=value))  + geom_jitter(alpha = 0.5,aes(colour=trt) )  + theme_bw() + facet_grid(target_id~.,scales="free") + theme(strip.text.y=element_text(angle=0),axis.text.x=element_text(angle=90)) + theme(legend.position="bottom") + labs(title=title)
  }  else {print("Specify genotype.");stop}
}


# using cpm data (tissue specific)
expression.pattern.Br.graph.exp1.v3.0annotation.3.cpm<-function(data=cpm.exp1.leaf.v3.0,target.genes.FDR,sample.description=sample.description.exp1,title="",geno,tissue.type="leaf"){
  if (tissue.type=="leaf") {
      data<-cpm.exp1.leaf.v3.0
    } else {data<-cpm.exp1.root.v3.0}
  print(paste("data is",data[1:10,]))
  print(paste("tissue.type is",tissue.type))
  
  
  data[is.na(data)] <- 0 #
  # select genes and add sample info
  # data.temp<-data %>% filter(target_id %in% target.genes) %>% gather(sample,value,-target_id) %>%
  #   inner_join(sample.description, by="sample") %>% filter(tissue==tissue.type)
  # 
  # using FDR
  data.temp<-data %>% 
    rename(target_id=transcript_ID) %>%
    filter(target_id %in% target.genes.FDR$genes) %>% 
    gather(sample,value,-target_id) %>%
    inner_join(sample.description, by="sample") %>% 
    filter(tissue==tissue.type) %>%
    right_join(target.genes.FDR, by=c("target_id"="genes")) %>%
    mutate(FDR=format(FDR,digits=2,scientific=TRUE)) %>%
    unite(target_id.FDR,target_id,FDR,sep="\n FDR ")

  # 
  if(geno=="both") { # needs to impove this
    # ggplot(data.temp, aes(x=genotype,y=value))  + geom_jitter(alpha = 0.5,aes(colour=trt,shape=tissue) )  + theme_bw() + facet_grid(target_id~tissue,scales="free") + theme(strip.text.y=element_text(angle=0),axis.text.x=element_text(angle=90)) + theme(legend.position="bottom") + labs(title=title)
    p<-data.temp %>% ggplot(aes(x=trt,y=value))  + 
      geom_jitter(alpha = 0.5,aes(colour=trt,shape=as.character(block)),width=0.2,size=3)  + 
      theme_bw() +
      facet_grid(target_id.FDR~genotype,scales="free") + 
      theme(strip.text.y=element_text(angle=0),axis.text.x=element_text(angle=90)) +
      theme(legend.position="bottom") + labs(title=" ")
    p
  } else if(geno=="FPsc"|geno=="R500") {
    data.temp %>% filter(genotype==geno) %>% ggplot(aes(x=trt,y=value))  + geom_jitter(alpha = 0.5,aes(colour=trt) )  + theme_bw() + facet_grid(target_id~.,scales="free") + theme(strip.text.y=element_text(angle=0),axis.text.x=element_text(angle=90)) + theme(legend.position="bottom") + labs(title=title)
  }  else {print("Specify genotype.");stop}
}

# 
# using cpm data (tissue specific), exp3
expression.pattern.Br.graph.exp3.v3.0annotation.3.cpm<-function(data=cpm.exp3.leaf.v3.0,target.genes.FDR,sample.description=sample.description.exp3,title="",geno,tissue.type="leaf"){
  if (tissue.type=="leaf") {
    data<-cpm.exp3.leaf.v3.0
  } else {data<-cpm.exp3.root.v3.0}
  print(paste("data is",data[1:10,]))
  print(paste("tissue.type is",tissue.type))
  
  data[is.na(data)] <- 0 #
  # select genes and add sample info
  # data.temp<-data %>% filter(target_id %in% target.genes) %>% gather(sample,value,-target_id) %>%
  #   inner_join(sample.description, by="sample") %>% filter(tissue==tissue.type)
  # 
  # using FDR
  data.temp<-data %>% 
    rename(target_id=transcript_ID) %>%
    filter(target_id %in% target.genes.FDR$genes) %>% 
    gather(sample,value,-target_id) %>%
    inner_join(sample.description, by="sample") %>% 
    filter(tissue==tissue.type) %>%
    right_join(target.genes.FDR, by=c("target_id"="genes")) %>%
    mutate(FDR=format(FDR,digits=2,scientific=TRUE)) %>%
    unite(target_id.FDR,target_id,FDR,sep="\n FDR ")
  
  # 
  if(geno=="both") { # needs to impove this
    # ggplot(data.temp, aes(x=genotype,y=value))  + geom_jitter(alpha = 0.5,aes(colour=trt,shape=tissue) )  + theme_bw() + facet_grid(target_id~tissue,scales="free") + theme(strip.text.y=element_text(angle=0),axis.text.x=element_text(angle=90)) + theme(legend.position="bottom") + labs(title=title)
    p<-data.temp %>% ggplot(aes(x=trt,y=value))  + 
      geom_jitter(alpha = 0.5,aes(colour=trt,shape=as.character(block)),width=0.2,size=3)  + 
      theme_bw() +
      facet_grid(target_id.FDR~genotype,scales="free") + 
      theme(strip.text.y=element_text(angle=0),axis.text.x=element_text(angle=90)) +
      theme(legend.position="bottom") + labs(title=" ")
    p
  } else if(geno=="FPsc"|geno=="R500") {
    data.temp %>% filter(genotype==geno) %>% ggplot(aes(x=trt,y=value))  + geom_jitter(alpha = 0.5,aes(colour=trt) )  + theme_bw() + facet_grid(target_id~.,scales="free") + theme(strip.text.y=element_text(angle=0),axis.text.x=element_text(angle=90)) + theme(legend.position="bottom") + labs(title=title)
  }  else {print("Specify genotype.");stop}
}
# timecourse (under construction)
expression.pattern.Br.graph.timecourse.v3.0annotation.cpm<-function(data=cpm.timecourse.v3.0,target.genes.FDR,sample.description=sample.description.timecourse,title=""){
  #print(paste("data is",data[1:10,]))
  #print(paste("tissue.type is root"))
  
  data[is.na(data)] <- 0 #
  # select genes and add sample info
  # data.temp<-data %>% filter(target_id %in% target.genes) %>% gather(sample,value,-target_id) %>%
  #   inner_join(sample.description, by="sample") %>% filter(tissue==tissue.type)
  # 
  # using FDR
  #target.genes.FDR <- target.genes.FDR #%>% slice(1)
  data.temp<-data  %>% 
    dplyr::rename(target_id=transcript_ID) %>%
    filter(target_id %in% target.genes.FDR$genes) %>% 
    gather(sample,value,-target_id) %>%
    inner_join(sample.description, by="sample")  %>%
    right_join(target.genes.FDR[,c("genes","FDR")], by=c("target_id"="genes")) %>%
    mutate(FDR=format(FDR,digits=2,scientific=TRUE)) %>%
    unite(target_id.FDR,target_id,FDR,sep="\n FDR ") 

# needs to impove this
    # ggplot(data.temp, aes(x=genotype,y=value))  + geom_jitter(alpha = 0.5,aes(colour=trt,shape=tissue) )  + theme_bw() + facet_grid(target_id~tissue,scales="free") + theme(strip.text.y=element_text(angle=0),axis.text.x=element_text(angle=90)) + theme(legend.position="bottom") + labs(title=title)
    p<-data.temp %>% ggplot(aes(x=soil_trt,y=value))  + 
      geom_jitter(alpha = 0.5,aes(colour=soil_trt,shape=as.character(block)),width=0.2,size=3)  + 
      theme_bw() +
      facet_grid(sampling_day~sampling_time,scales="free") + 
      theme(strip.text.y=element_text(angle=0),axis.text.x=element_text(angle=90)) +
      theme(legend.position="bottom") + labs(title=data.temp$target_id.FDR[1])
    p
}

# timecourse (under construction)
expression.pattern.Br.graph.timecourse.v3.0annotation.cpm.2<-function(data=cpm.timecourse.v3.0,target.genes,sample.description=sample.description.timecourse,title=""){
  #print(paste("data is",data[1:10,]))
  #print(paste("tissue.type is root"))
  data[is.na(data)] <- 0 #
  data.temp<-data %>% dplyr::filter(transcript_ID %in% target.genes$transcript_ID) %>% 
    gather(sample,value,-transcript_ID) %>%
    inner_join(sample.description, by="sample")
  # needs to impove this
  # ggplot(data.temp, aes(x=genotype,y=value))  + geom_jitter(alpha = 0.5,aes(colour=trt,shape=tissue) )  + theme_bw() + facet_grid(target_id~tissue,scales="free") + theme(strip.text.y=element_text(angle=0),axis.text.x=element_text(angle=90)) + theme(legend.position="bottom") + labs(title=title)
  p<-data.temp %>% ggplot(aes(x=soil_trt,y=value))  + 
    geom_jitter(alpha = 0.5,aes(colour=transcript_ID,shape=as.character(block)),width=0.2,size=3)  + 
    theme_bw() +
    facet_grid(sampling_day~sampling_time,scales="free") + 
    theme(strip.text.y=element_text(angle=0),axis.text.x=element_text(angle=90)) +
    theme(legend.position="bottom") + labs(title=data.temp$target_id.FDR[1])
  p
}

# logFC timecourse
expression.pattern.Br.graph.timecourse.v3.0annotation.logFC<-function(data=cpm.timecourse.v3.0.logFC,target.genes,title="",subset.data="only_two_afternoon"){
  #print(paste("data is",data[1:10,]))
  #print(paste("tissue.type is root"))
  data[is.na(data)] <- 0 #
  data.temp<-data  %>% dplyr::filter(transcript_ID %in% target.genes) 
  
  # if (2-afternoon=TRUE)
  if (subset.data=="only_two_afternoon") {
    p<-data.temp %>% ggplot(aes(x=sampling_day,y=logFC))  + 
      geom_boxplot(alpha = 0.5)  + 
      theme_bw() +
      theme(strip.text.y=element_text(angle=0),axis.text.x=element_text(angle=90)) +
      theme(legend.position="bottom") + labs(title=title)
    p
  } else {print("Define subset.data other than only_two_afternoon.")}
}





