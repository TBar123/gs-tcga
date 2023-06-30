##load packages

library(shiny)
library(ggplot2)
library(data.table)
library(fgsea)
library(dplyr)
library(abind)
library(shinyscreenshot)
library(shinyjs)
library("scales")
library(tidyr)
library(survival)
library(survminer)
library(ggplot2)
library(plyr)
##load data
Hallmark <- gmtPathways("h.all.v7.4.symbols.gmt")
C6<-gmtPathways("c6.all.v7.4.symbols.gmt")
C3<-gmtPathways("c3.tft.v7.4.symbols.gmt")
C7<-gmtPathways("c7.all.v7.4.symbols.gmt")

paths<-read.csv("C2_paths.csv")
pathws<-paths$pathway

C2<-gmtPathways("c2.all.v7.4.symbols.gmt")
tcga<-readRDS("tcgaa.rds")

pathway_choices<-c("Hallmark","C2","C3","C7","C6")
pathway_choices<-as.data.frame(pathway_choices)

colheads<-colnames(tcga)
colheads<-gsub("mrna_*","",colheads)
colnames(tcga)<-colheads
inputnames<-colheads[3:20504]



####################step 1 - for each tumour type remove genes in that tumour type
#which have no tumours with 10 or more reads
###then scale. 

a<-split(tcga,f=tcga$cdr_type)


###work out quartile then filter for ones where quartile bigger than 10
i<-1
b<-list()
for(i in 1:32){
  tum<-a[[i]]
  tumour_type<-tum[2,2]
  tum<-tum[,c(-1,-2)]
  tumquar<-apply( tum , 2 , quantile , probs = 0.75 , na.rm = TRUE )
  tumquar<-as.data.frame(tumquar)
  tumquar_filt<-tumquar%>%
    filter(tumquar>=10)
  b[[i]]<-tumquar_filt
  names(b[[i]])<-tumour_type
}

##now want to select genes from big tcga data that are in this list

##extract list of names of genes which have over 10 reads
c<-list()
i<-1
for(i in 1:32){
  tum<-b[[i]]
  genes<-rownames(tum)
  c[[i]]<-genes
}



#select these genes rom the OG TCGA list
list_tcga<-split(tcga,f=tcga$cdr_type)
i<-1
for(i in 1:32){
  tum<-list_tcga[[i]]
  genes<-c[[i]]
  genes2<-c("cdr_bcr_patient_barcode","cdr_type",genes)
  tum<-tum[names(tum) %in% genes2]
  tum<-tum%>%
    ##scale to 0 to 1
    mutate(across(where(is.numeric),~ rescale(.x),na.rm= TRUE))
  list_tcga[[i]]<-tum
}

##save it

saveRDS(list_tcga,"preprocessed_final.RDS")

data<-readRDS("preprocessed_final.rds")


####average score for each tumour for each gene set. 

####average score for each tumour for each gene set. ##########hallmark version

j<-1
k<-1

b<-list()
d<-list()

for (k in 1:32){
  tum<-data[[k]]
  tum = tum[!duplicated(tum$cdr_bcr_patient_barcode),]
  tumour_type<-tum[2]
  barcode<-tum[[1]]
  tum<-tum[,c(-1,-2)]
  len<-nrow(tum)
  for(j in 1:50){
    paths_used<-rep(names(Hallmark[j]),len)
    tum2<-tum%>%
      select(any_of(Hallmark[[j]]))
    mns<-rowMeans(tum2,na.rm=T)
    c<-cbind(barcode,tumour_type,mns,paths_used)
    b[[j]]<-c
    
  }
  
  d[[k]]<-bind_rows(b)
}



hallmark_means<-bind_rows(d)

saveRDS(hallmark_means,file="hallmark_means.rds")
hallmark_means<-readRDS("hallmark_means.rds")

hallmark_means$mns<-as.numeric(hallmark_means$mns)
##make wide
hallmark_means<-hallmark_means%>%
  pivot_wider(names_from = "paths_used",values_from = "mns")



###append cinical data
cdr_clinicaldata<-fread('TCGA-CDR-SupplementalTableS1_sampledata.tsv', verbose=TRUE, data.table=FALSE, header=TRUE)

hallmark_means<-hallmark_means%>%
  left_join(cdr_clinicaldata,by=c("barcode"="bcr_patient_barcode"))

saveRDS(hallmark_means,file="hallmark_means_and_clinical.rds")



############################survival things

hall<-colnames(hallmark_means)
hall<-hall[3:52]

tt<-unique(hallmark_means$cdr_type)


##########work out the cut point
i<-1
p<-list()


for(i in 1:32){
toi<-hallmark_means%>%
  filter(cdr_type==tt[i])
tumtype<-toi$cdr_type
h_cut<-surv_cutpoint(toi,time = "OS.time", event = "OS", variables = hall, minprop=0.1, progressbar = TRUE)
##put the cutpoints on the data
h_surv_cat<-surv_categorize(h_cut)
h_surv_cat<-as.data.frame(h_surv_cat)
together<-cbind(tumtype,h_surv_cat)
p[[i]]<-together
}




##stick all the tumour types into a giant document
h_surv_cat<-rbind(p[[1]],p[[2]],p[[3]],p[[4]],p[[5]],p[[6]],p[[7]],p[[8]],p[[9]],p[[10]],p[[11]],p[[12]],p[[13]]
                  ,p[[14]],p[[15]],p[[16]],p[[17]],p[[18]],p[[19]],p[[20]],p[[21]],p[[22]],p[[23]],p[[24]],p[[25]],p[[26]]
                  ,p[[27]],p[[28]],p[[29]],p[[30]],p[[31]],p[[32]])



saveRDS(h_surv_cat,file="Hallmark_surv_categories.rds")


###i also need to know a list of the cutpoints

hall<-colnames(hallmark_means)
hall<-hall[3:52]

tt<-unique(hallmark_means$cdr_type)

i<-1
cutsya<-list()
h_cutpts<-list()
y<-1


###get the cutpoints
for(i in 1:32){
  toi<-hallmark_means%>%
    filter(cdr_type==tt[i])
  tumtype<-toi$cdr_type
  h_cut<-surv_cutpoint(toi,time = "OS.time", event = "OS", variables = hall, minprop=0.1, progressbar = TRUE)
  ##extract the cut points
  for(y in 1:length(hall)){
    pt<-h_cut[[y]]$estimate
    nm<-hall[y]
    both<-cbind(nm,pt)
    both<-as.data.frame(both)
    h_cutpts[[y]]<-both
    
  }
  
  bound<-bind_rows(h_cutpts)
  bound$tumour_type<-rep(tumtype[1],nrow(bound))
  cutsya[[i]]<-bound
}

hall_cutpoints<-bind_rows(cutsya)

saveRDS(hall_cutpoints,"hallmark_cuts.RDS")

################### do a 50 50 split


hall<-colnames(hallmark_means)
hall<-hall[3:52]

tt<-unique(hallmark_means$cdr_type)

i<-1
p<-1
fiftysplit<-list()
paths_list<-list()
### categorise into high and low around median
for(i in 1:32){
  toi<-hallmark_means%>%
    filter(cdr_type==tt[i])
  tumtype<-toi$cdr_type
  for(p in 1:50){
  path<-hall[p]
  pathselected<-toi%>%select(barcode,path)
  colnames(pathselected)<-c("barcode","path")
  pathselected$path<-as.numeric(pathselected$path)
  
  pathselected<-pathselected%>%
    mutate(category=ifelse(path>=median(path),"high","low"))
  
  
  colnames(pathselected)<-c("barcode",paste0(path,"_value"),paste0(path))
 
  paths_list[[p]]<-pathselected
}
  fiftysplit[[i]]<-join_all(paths_list, by="barcode")
  fiftysplit[[i]]<-cbind(fiftysplit[[i]],tumtype)

  }



h_fiftyfifty<-bind_rows(fiftysplit)
saveRDS(h_fiftyfifty,"h_fifty.rds")
##stick clinical on
hfifty_fifty<-readRDS("h_fifty.rds")
hfifty_fifty<-hfifty_fifty%>%
  left_join(cdr_clinicaldata,by=c("barcode"="bcr_patient_barcode"))
saveRDS(hfifty_fifty,"h_fifty.rds")

a<-readRDS("h_fifty.rds")


####################################################C2


##############################C2
########################################################################################
#######################################################################################
data<-readRDS("preprocessed_final.rds")

j<-1
k<-1

b<-list()
d<-list()
##C2 means
for (k in 1:32){
  tum<-data[[k]]
  tum = tum[!duplicated(tum$cdr_bcr_patient_barcode),]
  tumour_type<-tum[2]
  barcode<-tum[[1]]
  tum<-tum[,c(-1,-2)]
  len<-nrow(tum)
  for(j in 1:6290){
    paths_used<-rep(names(C2[j]),len)
    tum2<-tum%>%
      select(any_of(C2[[j]]))
    mns<-rowMeans(tum2,na.rm=T)
    c<-cbind(barcode,tumour_type,mns,paths_used)
    b[[j]]<-c
    
  }
  
  d[[k]]<-bind_rows(b)
}



c2_means<-bind_rows(d)

saveRDS(c2_means,file="c2_means.rds")
a<-c2_means$mns
a<-as.data.frame(a)
c2_means<-readRDS("c2_means.rds")

c2_means$mns<-as.numeric(c2_means$mns)

c2_means<-c2_means%>%
  pivot_wider(names_from = "paths_used",values_from = "mns")



###append cinical data
cdr_clinicaldata<-fread('TCGA-CDR-SupplementalTableS1_sampledata.tsv', verbose=TRUE, data.table=FALSE, header=TRUE)

c2_means<-c2_means%>%
  left_join(cdr_clinicaldata,by=c("barcode"="bcr_patient_barcode"))

saveRDS(c2_means,file="c2_means_and_clinical.rds")

c2_means<-readRDS("c2_means_and_clinical.rds")



c2_paths<-colnames(c2_means)
c2_paths<-c2_paths[3:6292]

##remove NA paths- all genes 0
c2_paths<-c2_paths[! c2_paths %in% c('REACTOME_ORGANIC_ANION_TRANSPORT',"RUNNE_GENDER_EFFECT_UP",
                                     "SALVADOR_MARTIN_PEDIATRIC_TBD_ANTI_TNF_THERAPY_NONRESPONDER_POST_TREATMENT_DN",
                                     "REACTOME_FORMATION_OF_ATP_BY_CHEMIOSMOTIC_COUPLING",
                                     "REACTOME_ASSEMBLY_OF_THE_ORC_COMPLEX_AT_THE_ORIGIN_OF_REPLICATION"
                                     ,"REACTOME_RELAXIN_RECEPTORS" ,
                                     "REACTOME_MUSCARINIC_ACETYLCHOLINE_RECEPTORS",
                                     "REACTOME_DOPAMINE_RECEPTORS" ,
                                     "WEBER_METHYLATED_ICP_IN_SPERM_UP",
                                     "REACTOME_FREE_FATTY_ACID_RECEPTORS",
                                     "REACTOME_DEFECTIVE_F9_ACTIVATION",
                                     "BRUNEAU_SEPTATION_ATRIAL",
                                     "WEBER_METHYLATED_LCP_IN_FIBROBLAST_UP",
                                     "WEBER_METHYLATED_ICP_IN_SPERM_UP",
                                     "REACTOME_AMINO_ACID_CONJUGATION",
                                     "REACTOME_CONJUGATION_OF_BENZOATE_WITH_GLYCINE" ,
                                     "REACTOME_REGULATION_OF_GENE_EXPRESSION_IN_EARLY_PANCREATIC_PRECURSOR_CELLS"  ,
                                     "REACTOME_OREXIN_AND_NEUROPEPTIDES_FF_AND_QRFP_BIND_TO_THEIR_RESPECTIVE_RECEPTORS",
                                     "REACTOME_SEROTONIN_AND_MELATONIN_BIOSYNTHESIS"   ,
                                     "REACTOME_SEROTONIN_AND_MELATONIN_BIOSYNTHESIS"              ,
                                     "REACTOME_REGULATION_OF_GENE_EXPRESSION_IN_ENDOCRINE_COMMITTED_NEUROG3_PROGENITOR_CELLS",
                                     "REACTOME_FREE_FATTY_ACID_RECEPTORS" ,
                                     "BRUNEAU_SEPTATION_ATRIAL",
                                     "TESAR_ALK_AND_JAK_TARGETS_MOUSE_ES_D4_DN" ,
                                     "WP_GASTRIC_ACID_PRODUCTION"    ,
                                     "REACTOME_CONJUGATION_OF_BENZOATE_WITH_GLYCINE",
                                     "REACTOME_REGULATION_OF_GENE_EXPRESSION_IN_ENDOCRINE_COMMITTED_NEUROG3_PROGENITOR_CELLS",
                                     "REACTOME_MELANIN_BIOSYNTHESIS",
                                     "REACTOME_MINERALOCORTICOID_BIOSYNTHESIS",
                                     "TESAR_ALK_AND_JAK_TARGETS_MOUSE_ES_D4_DN",
                                     "WP_SECRETION_OF_HYDROCHLORIC_ACID_IN_PARIETAL_CELLS",
                                     "REACTOME_TACHYKININ_RECEPTORS_BIND_TACHYKININS" ,
                                     "WEBER_METHYLATED_ICP_IN_SPERM_DN" 
)]




tt<-unique(c2_means$cdr_type)

i<-1
p<-list()
##get the cutpoints and put them in
for(i in 1:32){
  toi<-c2_means%>%
    filter(cdr_type==tt[i])
  tumtype<-toi$cdr_type
  c2_cut<-surv_cutpoint(toi,time = "OS.time", event = "OS", variables = c2_paths, minprop=0.1, progressbar = TRUE)
  
  c2_surv_cat<-surv_categorize(c2_cut)
  c2_surv_cat<-as.data.frame(c2_surv_cat)
  together<-cbind(tumtype,c2_surv_cat)
  p[[i]]<-together
}



##stick all tumour types together

c2_surv_cat<-rbind(p[[1]],p[[2]],p[[3]],p[[4]],p[[5]],p[[6]],p[[7]],p[[8]],p[[9]],p[[10]],p[[11]],p[[12]],p[[13]]
                  ,p[[14]],p[[15]],p[[16]],p[[17]],p[[18]],p[[19]],p[[20]],p[[21]],p[[22]],p[[23]],p[[24]],p[[25]],p[[26]]
                  ,p[[27]],p[[28]],p[[29]],p[[30]],p[[31]],p[[32]])



saveRDS(c2_surv_cat,file="c2_surv_categories.rds")





i<-1
cutsya<-list()
c2_cutpts<-list()
y<-1


##get the cutpoints themselves in a list
for(i in 1:32){
  toi<-c2_means%>%
    filter(cdr_type==tt[i])
  tumtype<-toi$cdr_type
  c2_cut<-surv_cutpoint(toi,time = "OS.time", event = "OS", variables = c2_paths, minprop=0.1, progressbar = TRUE)
  
  for(y in 1:length(c2_paths)){
    pt<-c2_cut[[y]]$estimate
    nm<-c2_paths[y]
    both<-cbind(nm,pt)
    both<-as.data.frame(both)
    c2_cutpts[[y]]<-both
    
  }
  
  bound<-bind_rows(c2_cutpts)
  bound$tumour_type<-rep(tumtype[1],nrow(bound))
  cutsya[[i]]<-bound
}

c2_cutpoints<-bind_rows(cutsya)

saveRDS(c2_cutpoints,"c2_cuts.RDS")


################### do a 50 50 split



i<-1
p<-1
fiftysplit<-list()
paths_list<-list()

for(i in 1:32){
  toi<-c2_means%>%
    filter(cdr_type==tt[i])
  tumtype<-toi$cdr_type
  for(p in 1:6264){
    path<-c2_paths[p]
    pathselected<-toi%>%select(barcode,path)
    colnames(pathselected)<-c("barcode","path")
    pathselected$path<-as.numeric(pathselected$path)

    
    pathselected<-pathselected%>%
      mutate(category=ifelse(path>=median(path),"high","low"))
    
    
    colnames(pathselected)<-c("barcode",paste0(path,"_value"),paste0(path))
    
    paths_list[[p]]<-pathselected
  }
  fiftysplit[[i]]<-join_all(paths_list, by="barcode")
  fiftysplit[[i]]<-cbind(fiftysplit[[i]],tumtype)
  
}

c2_fiftyfifty<-bind_rows(fiftysplit)
saveRDS(c2_fiftyfifty,"c2_fifty.rds")

c2_fiftyfifty<-readRDS("c2_fifty.rds")
c2_fiftyfifty<-c2_fiftyfifty%>%
  left_join(cdr_clinicaldata,by=c("barcode"="bcr_patient_barcode"))
saveRDS(c2_fiftyfifty,"c2_fiftyfifty.rds")

####################################################c3


##############################c3
########################################################################################
#######################################################################################
data<-readRDS("preprocessed_final.rds")

j<-1
k<-1

b<-list()
d<-list()
##get the means
for (k in 1:32){
  tum<-data[[k]]
  tum = tum[!duplicated(tum$cdr_bcr_patient_barcode),]
  tumour_type<-tum[2]
  barcode<-tum[[1]]
  tum<-tum[,c(-1,-2)]
  len<-nrow(tum)
  for(j in 1:1133){
    paths_used<-rep(names(C3[j]),len)
    tum2<-tum%>%
      select(any_of(C3[[j]]))
    mns<-rowMeans(tum2,na.rm=T)
    c<-cbind(barcode,tumour_type,mns,paths_used)
    b[[j]]<-c
    
  }
  
  d[[k]]<-bind_rows(b)
}



c3_means<-bind_rows(d)

saveRDS(c3_means,file="c3_means.rds")
c3_means<-readRDS("c3_means.rds")

c3_means$mns<-as.numeric(c3_means$mns)

c3_means<-c3_means%>%
  pivot_wider(names_from = "paths_used",values_from = "mns")



###append cinical data
cdr_clinicaldata<-fread('TCGA-CDR-SupplementalTableS1_sampledata.tsv', verbose=TRUE, data.table=FALSE, header=TRUE)

c3_means<-c3_means%>%
  left_join(cdr_clinicaldata,by=c("barcode"="bcr_patient_barcode"))

saveRDS(c3_means,file="c3_means_and_clinical.rds")




c3_paths<-colnames(c3_means)
c3_paths<-c3_paths[3:1135]

tt<-unique(c3_means$cdr_type)

i<-1
p<-list()


##remove paths
c3_paths<-c3_paths[! c3_paths %in% c("LAMTOR5_TARGET_GENES"   ,
                                     "RARB_TARGET_GENES", "CASP3_TARGET_GENES",
                                     "F2RL1_TARGET_GENES")]
##categrorise by cutpoint
for(i in 1:32){
  toi<-c3_means%>%
    filter(cdr_type==tt[i])
  tumtype<-toi$cdr_type
  c3_cut<-surv_cutpoint(toi,time = "OS.time", event = "OS", variables = c3_paths, minprop=0.1, progressbar = TRUE)
  
  c3_surv_cat<-surv_categorize(c3_cut)
  c3_surv_cat<-as.data.frame(c3_surv_cat)
  together<-cbind(tumtype,c3_surv_cat)
  p[[i]]<-together
}



##put tumour types together

c3_surv_cat<-rbind(p[[1]],p[[2]],p[[3]],p[[4]],p[[5]],p[[6]],p[[7]],p[[8]],p[[9]],p[[10]],p[[11]],p[[12]],p[[13]]
                   ,p[[14]],p[[15]],p[[16]],p[[17]],p[[18]],p[[19]],p[[20]],p[[21]],p[[22]],p[[23]],p[[24]],p[[25]],p[[26]]
                   ,p[[27]],p[[28]],p[[29]],p[[30]],p[[31]],p[[32]])



saveRDS(c3_surv_cat,file="c3_surv_categories.rds")


###cut points
i<-1
cutsya<-list()
c3_cutpts<-list()
y<-1

##extract cutpoints

for(i in 1:32){
  toi<-c3_means%>%
    filter(cdr_type==tt[i])
  tumtype<-toi$cdr_type
  c3_cut<-surv_cutpoint(toi,time = "OS.time", event = "OS", variables = c3_paths, minprop=0.1, progressbar = TRUE)
  
  for(y in 1:length(c3_paths)){
    pt<-c3_cut[[y]]$estimate
    nm<-c3_paths[y]
    both<-cbind(nm,pt)
    both<-as.data.frame(both)
    c3_cutpts[[y]]<-both
    
  }
  
  bound<-bind_rows(c3_cutpts)
  bound$tumour_type<-rep(tumtype[1],nrow(bound))
  cutsya[[i]]<-bound
}

c3_cutpoints<-bind_rows(cutsya)

saveRDS(c3_cutpoints,"c3_cuts.RDS")


################### do a 50 50 split



i<-1
p<-1
fiftysplit<-list()
paths_list<-list()

c3_paths<-colnames(c3_means)
c3_paths<-c3_paths[3:1135]


c3_paths<-c3_paths[! c3_paths %in% c("LAMTOR5_TARGET_GENES"   ,
                                     "RARB_TARGET_GENES", "CASP3_TARGET_GENES",
                                     "F2RL1_TARGET_GENES")]


for(i in 1:32){
  toi<-c3_means%>%
    filter(cdr_type==tt[i])
  tumtype<-toi$cdr_type
  for(p in 1:1129){
    path<-c3_paths[p]
    pathselected<-toi%>%select(barcode,path)
    colnames(pathselected)<-c("barcode","path")
    pathselected$path<-as.numeric(pathselected$path)
    

    pathselected<-pathselected%>%
      mutate(category=ifelse(path>=median(path),"high","low"))
    
    
    colnames(pathselected)<-c("barcode",paste0(path,"_value"),paste0(path))
    
    paths_list[[p]]<-pathselected
  }
  fiftysplit[[i]]<-join_all(paths_list, by="barcode")
  fiftysplit[[i]]<-cbind(fiftysplit[[i]],tumtype)
  
}

c3_fiftyfifty<-bind_rows(fiftysplit)
saveRDS(c3_fiftyfifty,"c3_fifty.rds")

c3_fiftyfifty<-readRDS("c3_fifty.rds")
c3_fiftyfifty<-c3_fiftyfifty%>%
  left_join(cdr_clinicaldata,by=c("barcode"="bcr_patient_barcode"))
saveRDS(c3_fiftyfifty,"c3_fiftyfifty.rds")

##############################c6
########################################################################################
#######################################################################################
data<-readRDS("preprocessed_final.rds")

j<-1
k<-1

b<-list()
d<-list()
##means
for (k in 1:32){
  tum<-data[[k]]
  tum = tum[!duplicated(tum$cdr_bcr_patient_barcode),]
  tumour_type<-tum[2]
  barcode<-tum[[1]]
  tum<-tum[,c(-1,-2)]
  len<-nrow(tum)
  for(j in 1:189){
    paths_used<-rep(names(C6[j]),len)
    tum2<-tum%>%
      select(any_of(C6[[j]]))
    mns<-rowMeans(tum2,na.rm=T)
    c<-cbind(barcode,tumour_type,mns,paths_used)
    b[[j]]<-c
    
  }
  
  d[[k]]<-bind_rows(b)
}



c6_means<-bind_rows(d)

saveRDS(c6_means,file="c6_means.rds")
c6_means<-readRDS("c6_means.rds")

c6_means$mns<-as.numeric(c6_means$mns)

c6_means<-c6_means%>%
  pivot_wider(names_from = "paths_used",values_from = "mns")



###append cinical data
cdr_clinicaldata<-fread('TCGA-CDR-SupplementalTableS1_sampledata.tsv', verbose=TRUE, data.table=FALSE, header=TRUE)

c6_means<-c6_means%>%
  left_join(cdr_clinicaldata,by=c("barcode"="bcr_patient_barcode"))

saveRDS(c6_means,file="c6_means_and_clinical.rds")



c6_paths<-colnames(c6_means)
c6_paths<-c6_paths[3:191]
#a<-as.data.frame(c6_paths)

tt<-unique(c6_means$cdr_type)

i<-1
p<-list()


##cutpoints and categorise

for(i in 1:32){
  toi<-c6_means%>%
    filter(cdr_type==tt[i])
  tumtype<-toi$cdr_type
  c6_cut<-surv_cutpoint(toi,time = "OS.time", event = "OS", variables = c6_paths, minprop=0.1, progressbar = TRUE)
  
  c6_surv_cat<-surv_categorize(c6_cut)
  c6_surv_cat<-as.data.frame(c6_surv_cat)
  together<-cbind(tumtype,c6_surv_cat)
  p[[i]]<-together
}




##stick tumour types together
c6_surv_cat<-rbind(p[[1]],p[[2]],p[[3]],p[[4]],p[[5]],p[[6]],p[[7]],p[[8]],p[[9]],p[[10]],p[[11]],p[[12]],p[[13]]
                   ,p[[14]],p[[15]],p[[16]],p[[17]],p[[18]],p[[19]],p[[20]],p[[21]],p[[22]],p[[23]],p[[24]],p[[25]],p[[26]]
                   ,p[[27]],p[[28]],p[[29]],p[[30]],p[[31]],p[[32]])



saveRDS(c6_surv_cat,file="c6_surv_categories.rds")

##extract cutpoints
i<-1
cutsya<-list()
c6_cutpts<-list()
y<-1



for(i in 1:32){
  toi<-c6_means%>%
    filter(cdr_type==tt[i])
  tumtype<-toi$cdr_type
  c6_cut<-surv_cutpoint(toi,time = "OS.time", event = "OS", variables = c6_paths, minprop=0.1, progressbar = TRUE)
  
  for(y in 1:length(c6_paths)){
    pt<-c6_cut[[y]]$estimate
    nm<-c6_paths[y]
    both<-cbind(nm,pt)
    both<-as.data.frame(both)
    c6_cutpts[[y]]<-both
    
  }
  
  bound<-bind_rows(c6_cutpts)
  bound$tumour_type<-rep(tumtype[1],nrow(bound))
  cutsya[[i]]<-bound
}

c6_cutpoints<-bind_rows(cutsya)

saveRDS(c6_cutpoints,"c6_cuts.RDS")




################### do a 50 50 split



i<-1
p<-1
fiftysplit<-list()
paths_list<-list()

c6_paths<-colnames(c6_means)
c6_paths<-c6_paths[3:191]



for(i in 1:32){
  toi<-c6_means%>%
    filter(cdr_type==tt[i])
  tumtype<-toi$cdr_type
  for(p in 1:189){
    path<-c6_paths[p]
    pathselected<-toi%>%select(barcode,path)
    colnames(pathselected)<-c("barcode","path")
    pathselected$path<-as.numeric(pathselected$path)
    
   
    pathselected<-pathselected%>%
      mutate(category=ifelse(path>=median(path),"high","low"))
    
    
    colnames(pathselected)<-c("barcode",paste0(path,"_value"),paste0(path))
    
    paths_list[[p]]<-pathselected
  }
  fiftysplit[[i]]<-join_all(paths_list, by="barcode")
  fiftysplit[[i]]<-cbind(fiftysplit[[i]],tumtype)
  
}

c6_fiftyfifty<-bind_rows(fiftysplit)
saveRDS(c6_fiftyfifty,"c6_fifty.rds")

c6_fiftyfifty<-readRDS("c6_fifty.rds")
c6_fiftyfifty<-c6_fiftyfifty%>%
  left_join(cdr_clinicaldata,by=c("barcode"="bcr_patient_barcode"))
saveRDS(c6_fiftyfifty,"c6_fiftyfifty.rds")

##############################c7
########################################################################################
#######################################################################################
data<-readRDS("preprocessed_final.rds")

j<-1
k<-1

b<-list()
d<-list()
##means
for (k in 1:32){
  tum<-data[[k]]
  tum = tum[!duplicated(tum$cdr_bcr_patient_barcode),]
  tumour_type<-tum[2]
  barcode<-tum[[1]]
  tum<-tum[,c(-1,-2)]
  len<-nrow(tum)
  for(j in 1:5219){
    paths_used<-rep(names(C7[j]),len)
    tum2<-tum%>%
      select(any_of(C7[[j]]))
    mns<-rowMeans(tum2,na.rm=T)
    c<-cbind(barcode,tumour_type,mns,paths_used)
    b[[j]]<-c
    
  }
  
  d[[k]]<-bind_rows(b)
}



c7_means<-bind_rows(d)

saveRDS(c7_means,file="c7_means.rds")
c7_means<-readRDS("c7_means.rds")

c7_means$mns<-as.numeric(c7_means$mns)

c7_means<-c7_means%>%
  pivot_wider(names_from = "paths_used",values_from = "mns")



###append cinical data
cdr_clinicaldata<-fread('TCGA-CDR-SupplementalTableS1_sampledata.tsv', verbose=TRUE, data.table=FALSE, header=TRUE)

c7_means<-c7_means%>%
  left_join(cdr_clinicaldata,by=c("barcode"="bcr_patient_barcode"))

saveRDS(c7_means,file="c7_means_and_clinical.rds")

c7_means<-readRDS("c7_means_and_clinical.rds")



c7_paths<-colnames(c7_means)
c7_paths<-c7_paths[3:5221]
a<-as.data.frame(c7_paths)

tt<-unique(c7_means$cdr_type)

i<-1
p<-list()

c7_paths<-c7_paths[! c7_paths %in% c("FISCHER_BLOOD_PLASMA_RVSV_EBOV_AGE_18_55YO_HIGH_DOSE_3DY_DN",
                                     "QI_CD4_POSITIVE_ALPHA_BETA_MEMORY_T_CELL_ZOSTAVAX_AGE_52_75YO_CD4_T_CELL_VS_NAIVE_CD4_T_CELL_7_TO_9DY_DN"
                                     ,"QI_NAIVE_T_CELL_ZOSTAVAX_AGE_52_75YO_CD4_T_CELL_VS_NAIVE_CD4_T_CELL_7_TO_9DY_DN" ,
                                     "QI_NAIVE_T_CELL_ZOSTAVAX_AGE_52_75YO_CD4_T_CELL_VS_NAIVE_CD4_T_CELL_7_TO_9DY_UP")]



##cutpoints and categorise

for(i in 1:32){
  toi<-c7_means%>%
    filter(cdr_type==tt[i])
  tumtype<-toi$cdr_type
  c7_cut<-surv_cutpoint(toi,time = "OS.time", event = "OS", variables = c7_paths, minprop=0.1, progressbar = TRUE)
  
  c7_surv_cat<-surv_categorize(c7_cut)
  c7_surv_cat<-as.data.frame(c7_surv_cat)
  together<-cbind(tumtype,c7_surv_cat)
  p[[i]]<-together
}



##put tumour types together

c7_surv_cat<-rbind(p[[1]],p[[2]],p[[3]],p[[4]],p[[5]],p[[6]],p[[7]],p[[8]],p[[9]],p[[10]],p[[11]],p[[12]],p[[13]]
                   ,p[[14]],p[[15]],p[[16]],p[[17]],p[[18]],p[[19]],p[[20]],p[[21]],p[[22]],p[[23]],p[[24]],p[[25]],p[[26]]
                   ,p[[27]],p[[28]],p[[29]],p[[30]],p[[31]],p[[32]])



saveRDS(c7_surv_cat,file="c7_surv_categories.rds")

i<-1
cutsya<-list()
c7_cutpts<-list()
y<-1

##extract cutpoints

for(i in 1:32){
  toi<-c7_means%>%
    filter(cdr_type==tt[i])
  tumtype<-toi$cdr_type
  c7_cut<-surv_cutpoint(toi,time = "OS.time", event = "OS", variables = c7_paths, minprop=0.1, progressbar = TRUE)
  
  for(y in 1:length(c7_paths)){
    pt<-c7_cut[[y]]$estimate
    nm<-c7_paths[y]
    both<-cbind(nm,pt)
    both<-as.data.frame(both)
    c7_cutpts[[y]]<-both
    
  }
  
  bound<-bind_rows(c7_cutpts)
  bound$tumour_type<-rep(tumtype[1],nrow(bound))
  cutsya[[i]]<-bound
}

c7_cutpoints<-bind_rows(cutsya)

saveRDS(c7_cutpoints,"c7_cuts.RDS")


################### do a 50 50 split



i<-1
p<-1
fiftysplit<-list()
paths_list<-list()





for(i in 1:32){
  toi<-c7_means%>%
    filter(cdr_type==tt[i])
  tumtype<-toi$cdr_type
  for(p in 1:5215){
    path<-c7_paths[p]
    pathselected<-toi%>%select(barcode,path)
    colnames(pathselected)<-c("barcode","path")
    pathselected$path<-as.numeric(pathselected$path)
    
   
    pathselected<-pathselected%>%
      mutate(category=ifelse(path>=median(path),"high","low"))
    
    
    colnames(pathselected)<-c("barcode",paste0(path,"_value"),paste0(path))
    
    
    paths_list[[p]]<-pathselected
  }
  fiftysplit[[i]]<-join_all(paths_list, by="barcode")
  fiftysplit[[i]]<-cbind(fiftysplit[[i]],tumtype)
  
}

c7_fiftyfifty<-bind_rows(fiftysplit)
saveRDS(c7_fiftyfifty,"c7_fifty.rds")

c7_fiftyfifty<-readRDS("c7_fifty.rds")
c7_fiftyfifty<-c7_fiftyfifty%>%
  left_join(cdr_clinicaldata,by=c("barcode"="bcr_patient_barcode"))
saveRDS(c7_fiftyfifty,"c7_fiftyfifty.rds")



##########putting things together

h_surv<-readRDS("Hallmark_surv_categories.rds")
write_fst(h_surv,"h_opt.fst")
c2_surv<-readRDS("C2_surv_categories.rds")
write_fst(c2_surv,"c2_opt.fst")
c3_surv<-readRDS("c3_surv_categories.rds")
write_fst(c3_surv,"c3_opt.fst")
c6_surv<-readRDS("c6_surv_categories.rds")
write_fst(c6_surv,"c6_opt.fst")
c7_surv<-readRDS("c7_surv_categories.rds")
write_fst(c7_surv,"c7_opt.fst")

h_50<-readRDS("h_fifty.rds")
write_fst(h_50,"h.fst")
c2_50<-readRDS("c2_fiftyfifty.rds")
write_fst(c2_50,"c2.fst")
c3_50<-readRDS("c3_fiftyfifty.rds")
write_fst(c3_50,"c3.fst")
c6_50<-readRDS("c6_fiftyfifty.rds")
write_fst(c6_50,"c6.fst")
c7_50<-readRDS("c7_fiftyfifty.rds")
write_fst(c7_50,"c7.fst")

h_cuts<-readRDS("hallmark_cuts.rds")
c2_cuts<-readRDS("c2_cuts.rds")
c3_cuts<-readRDS("c3_cuts.rds")
c6_cuts<-readRDS("c6_cuts.rds")
c7_cuts<-readRDS("c7_cuts.RDS")

cuts<-list(h_cuts,c2_cuts,c3_cuts,c6_cuts,c7_cuts)
names(cuts)<-c("Hallmark","C2","C3","C6","C7")

#dataa<-list(h_surv,c2_surv,c3_surv,c6_surv,c7_surv)
#dataa50<-list(h_50,c2_50,c3_50,c6_50,c7_50)
#names(dataa)<-c("Hallmark","C2","C3","C6","C7")
#names(dataa50)<-c("Hallmark","C2","C3","C6","C7")

#saveRDS(dataa,"opt_surv.RDS")
#saveRDS(dataa50,"fifty_surv.RDS")
saveRDS(cuts,"cutpoints.RDS")

##load in the pathway data- paths and the genes in it
Hallmark <- gmtPathways("h.all.v7.4.symbols.gmt")
C6<-gmtPathways("c6.all.v7.4.symbols.gmt")
C3<-gmtPathways("c3.tft.v7.4.symbols.gmt")
C7<-gmtPathways("c7.all.v7.4.symbols.gmt")
C2<-gmtPathways("c2.all.v7.4.symbols.gmt")

all_paths<-list(Hallmark,C2,C3,C6,C7)
names(all_paths)<-c("Hallmark","C2","C3","C6","C7")
saveRDS(all_paths,"all_paths.RDS")