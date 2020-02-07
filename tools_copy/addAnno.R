# functions for adding Arabidopsis thaliana annotation to edgeR topTag output
# Note: use Brassica rapa v3.0 annotation
# Note2: annotated files were written in csv under "../output" directory.
## prep
# annotation file for v3.0annotation
Br.v3.0.At.BLAST <- read_csv(file.path("..","Annotation_copy","output","v3.0annotation","Brapa_v3.0_annotated.csv")) 
# This annotation is redundant with name (Br grene). Eg 
Br.v3.0.At.BLAST %>% filter(name=="BraA01g040570.3C")
# reduce the redundancy (112418)
Br.v3.0anno.At.BLAST.highscore <- Br.v3.0.At.BLAST %>% group_by(name) %>% arrange(desc(score)) %>% dplyr::slice(1)
# function for adding annotation
## get object name https://stackoverflow.com/questions/14577412/how-to-convert-variable-object-name-into-string
# learning how to get object name
# myfunc <- function(v1) {  deparse(substitute(v1))}
# myfunc(foo)
# adding annotation and write_csv adding ".v3.0anno.csv" to the object name.
addAnno<-function(DGE) {
  temp<-left_join(DGE %>%
                    rownames_to_column(var="genes"),Br.v3.0anno.At.BLAST.highscore,by=c(genes="name")) %>%
    dplyr::select(genes,names(DGE),AGI, At_symbol, At_short_description, perc_ID);
  print(deparse(substitute(DGE)));
  write_csv(temp, path=file.path("..","output",paste(deparse(substitute(DGE)),".v3.0anno.csv",sep="")));
  return(temp)
  } 

# an example (found in 06_two_factor_models.Rmd)
# twoafternoon.trtlive.DEGs.all.v3.0anno <- addAnno(twoafternoon.trtlive.DEGs.all)
