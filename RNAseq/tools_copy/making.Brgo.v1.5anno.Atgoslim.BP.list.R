# common tools for Brassica microbiome projects
## home dir is this file directory
  library(tidyverse)
  ## making Brgo.v1.5anno.Atgoslim.BP.list (run once)  (bug fixed with max_target_seq option problem in blastn)
  ### prerequisit
  # Br.v1.5anno.At.BLAST<-readr::read_csv(file.path("Annotation","input","v1.5annotation","Brapa1.5Davis_vs_At_dc-megablast_out.csv"),col_names=FALSE) # BLASTed by Julin
  # colnames(Br.v1.5anno.At.BLAST) <- c("query","subject","perc_ID","aln_length","mismatch","gap_open","qstart","qend","sstart","send","eval","score")
  # head(Br.v1.5anno.At.BLAST) # no! needs to use Brapa_V1.5Davis_annotated.csv?
  # 
  # Br.v1.5anno.At.BLAST <- Br.v1.5anno.At.BLAST %>% tidyr::separate(subject, c("AGI","splicing"),sep="\\.") #%>% select(AGI)
  # summary(Br.v1.5anno.At.BLAST)
  Br.v1.5anno.At.BLAST<-readr::read_csv(file.path("..","Annotation","output","v1.5annotation","Brapa_V1.5Davis_annotated.csv"),col_names=TRUE) # BLASTed by Julin
  summary(Br.v1.5anno.At.BLAST)
  # there are multiple Br genes found
  Br.v1.5anno.At.BLAST %>% group_by(name) %>% summarise(n=n()) %>% arrange(desc(n))
  # pick only one with highest score within the same Br gene name
  Br.v1.5anno.At.BLAST.highscore <- Br.v1.5anno.At.BLAST %>% group_by(name) %>% arrange(desc(score)) %>% slice(1)
  Br.v1.5anno.At.BLAST.highscore %>% group_by(name) %>% summarise(n=n()) %>% arrange(desc(n))
  
  ### option (find ATHB2("AT4G16780") and phyA ("AT1G09570"))
  #Br.v2.5anno.At.BLAST.s <- Br.v2.5anno.At.BLAST %>% tidyr::separate(subject, c("AGI","splicing"),sep="\\.")
  # update Arabidopsis GO annotation (030419)
  download.file("https://www.arabidopsis.org/download_files/Subscriber_Data_Releases/TAIR_Data_20181231/ATH_GO_GOSLIM.txt.gz",destfile ="ATH_GO_GOSLIM.txt.gz") 
  system("unzip -c ATH_GO_GOSLIM.txt.gz >  ATH_GO_GOSLIM.txt") # this is not gunzip.
  Atgoslim.TAIR<-read_tsv("ATH_GO_GOSLIM.txt",skip=2,col_names = FALSE)
  Atgoslim.TAIR %>% filter(X8=="P",X1=="AT4G38360") %>% dplyr::select(1,5,6,8,9,10)
  Atgoslim.TAIR.BP <-Atgoslim.TAIR%>% filter(X8=="P")
  Atgoslim.TAIR.BP.list<-tapply(Atgoslim.TAIR.BP$X6,Atgoslim.TAIR.BP$X1,c)
  save(Atgoslim.TAIR.BP.list,file=file.path("..","Annotation","input","Atgoslim.TAIR.BP.list.Rdata"))
  system("rm ATH_GO_GOSLIM.txt*") # removce both files
  # read At GO list
  load(file.path("..","Annotation","input","Atgoslim.TAIR.BP.list.Rdata"))
  
  # asign At GO into corresponding Br genes
  Brgo.v1.5anno.Atgoslim.BP.list<-list()
  
  for(i in 1:length(Br.v1.5anno.At.BLAST.highscore$name)) {
    if(is.null(Atgoslim.TAIR.BP.list[[as_vector(Br.v1.5anno.At.BLAST.highscore[i,"AGI"])]])) next else {
      Brgo.v1.5anno.Atgoslim.BP.list[[i]]<-Atgoslim.TAIR.BP.list[[as_vector(Br.v1.5anno.At.BLAST.highscore[i,"AGI"])]]
      names(Brgo.v1.5anno.Atgoslim.BP.list)[[i]]<-as_vector(Br.v1.5anno.At.BLAST.highscore[i,"name"])
    }
  }
  table(sapply(Brgo.v1.5anno.Atgoslim.BP.list,is.null))
  Brgo.v1.5anno.Atgoslim.BP.list<-Brgo.v1.5anno.Atgoslim.BP.list[!sapply(Brgo.v1.5anno.Atgoslim.BP.list,is.null)]
  table(sapply(Brgo.v1.5anno.Atgoslim.BP.list,is.null))
  table(sapply(Brgo.v1.5anno.Atgoslim.BP.list,length))
  
  save(Brgo.v1.5anno.Atgoslim.BP.list,file=file.path("..","Annotation","output","v1.5annotation","Brgo.v1.5anno.Atgoslim.BP.list.Rdata"))

  # comparison with an old version
load(file.path("..","Annotation","output","v1.5annotation","Brgo.v1.5anno.Atgoslim.BP.list_old.Rdata"))
table(sapply(Brgo.v1.5anno.Atgoslim.BP.list,length)) # different from new one. Old one is based on old go slim (Sep 2018?)



