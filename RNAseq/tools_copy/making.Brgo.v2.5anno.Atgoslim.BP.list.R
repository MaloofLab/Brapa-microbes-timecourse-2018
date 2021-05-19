# common tools for Brassica microbiome projects
## home dir is
# "/Volumes/data_work/Data8/NGS_related/Brassica_microbiome/Brassica_microbiome_Kazu"
# if(getwd()=="/Volumes/data_work/Data8/NGS_related/Brassica_microbiome/Brassica_microbiome_Kazu") {
  library(tidyverse)
  ## making Brgo.v2.5anno.Atgoslim.BP.list (run once)
  ### prerequisit
  Br.v2.5anno.At.BLAST<-readr::read_csv(file.path("common","Annotation","input","v2.5annotation","Brapa_V2.5_annotated.csv")) # BLASTed by Julin
  # 
  ### option (find ATHB2("AT4G16780") and phyA ("AT1G09570"))
  #Br.v2.5anno.At.BLAST.s <- Br.v2.5anno.At.BLAST %>% tidyr::separate(subject, c("AGI","splicing"),sep="\\.")
  # read At GO list
  load(file.path("common","Annotation","input","Atgoslim.TAIR.BP.list.Rdata"))
  # asign At GO into corresponding Br genes
  Brgo.v2.5anno.Atgoslim.BP.list<-list()
  for(i in 1:length(Br.v2.5anno.At.BLAST$name)) {
    if(is.null(Atgoslim.TAIR.BP.list[[as_vector(Br.v2.5anno.At.BLAST[i,"AGI"])]])) next else {
      Brgo.v2.5anno.Atgoslim.BP.list[[i]]<-Atgoslim.TAIR.BP.list[[as_vector(Br.v2.5anno.At.BLAST[i,"AGI"])]]
      names(Brgo.v2.5anno.Atgoslim.BP.list)[[i]]<-as_vector(Br.v2.5anno.At.BLAST[i,"name"])
    }
  }
  table(sapply(Brgo.v2.5anno.Atgoslim.BP.list,is.null))
  Brgo.v2.5anno.Atgoslim.BP.list.list<-Brgo.v2.5anno.Atgoslim.BP.list[!sapply(Brgo.v2.5anno.Atgoslim.BP.list,is.null)]
  table(sapply(Brgo.v2.5anno.Atgoslim.BP.list,is.null))
  #save(Brgo.v2.5anno.Atgoslim.BP.list,file="../../Annotation/input/Brgo.v2.5anno.Atgoslim.BP.list.Rdata")
  save(Brgo.v2.5anno.Atgoslim.BP.list,file=file.path("common","Annotation","output","v2.5annotation","Brgo.v2.5anno.Atgoslim.BP.list.Rdata"))
#} else {print("change directory");stop}


