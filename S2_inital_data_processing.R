
require(data.table)
require(stringr)
library(reshape)
library(survival)
library(ggplot2)
library(fst)
library(dplyr)

###make the clinical data file contain only what we need
cdr_clinicaldata<-fread('TCGA-CDR-SupplementalTableS1_sampledata.csv', verbose=TRUE, data.table=FALSE, header=TRUE)

cdr_clinicaldata<-cdr_clinicaldata%>%
  select(bcr_patient_barcode,type,OS.time,OS)
#write.table(cdr_clinicaldata,'clinical_data.tsv',sep="\t")


###load in rna data 
rnaseq<-fread('EBPlusPlusAdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp_transposed.tsv', verbose=TRUE, data.table=FALSE, header=TRUE) #Load transposed tcga rnaseq data
colnames(rnaseq)<-sub("\\|.*", "", colnames(rnaseq))
colnames(rnaseq)<-paste0("mrna_",colnames(rnaseq)) #Appends mrna_ before all colnames in table


rnaseq_genes<-colnames(rnaseq[,2:length(rnaseq)]) #Makes a list of all genes (for later column selection)
rnaseq$mrna_gene_id<-substring(rnaseq$mrna_gene_id,1,15) #Makes the TCGA barcode compatible with the methylcs table tcga barcode (cuts off everything to the right of the 15th character which is redundant information)


####make sure its just primary tumour
rnaseq_primary_tumour<-rnaseq[substring(rnaseq$mrna_gene_id,14,15)=='01',]
rnaseq_primary_tumour$patient_id<-substring(rnaseq_primary_tumour$mrna_gene_id,1,12) # Create a Patient ID column for subsequent merge with clinical data (ie cut off everything from the TCGA barcode to the right of the 12th character - ie we don't need sample type any more)


tcga<-merge(cdr_clinicaldata, rnaseq_primary_tumour, by.x="bcr_patient_barcode", by.y="patient_id", all.x=FALSE, all.y=FALSE)

tcga<-rename(tcga,cdr_bcr_patient_barcode=bcr_patient_barcode)
tcga<-rename(tcga,cdr_type=type)
names<-colnames(tcga)

names<-gsub("mrna_*","",names)
colnames(tcga)<-names

tcga_w_barcode<-tcga[,c(1:2,35:20536)]


####ensure gene symbols match between tcga data and msigDB
tcga_genes<-colnames(tcga_w_barcode)

##use chip file
chip<-read.table("Human_Gene_Symbol_with_Remapping_MSigDB.v2023.2.Hs.chip",sep="\t", header=TRUE,quote = "")

tcga_genes<-as.data.frame(tcga_genes)

##merge chip with current gene names
genes_chip<-merge(tcga_genes,chip,by.x="tcga_genes",by.y="Probe.Set.ID")
genes_chip<-t(genes_chip)
colnames(genes_chip)<-genes_chip[1, ]
genes_chip<-as.data.frame(genes_chip)
genes_chip$cdr_bcr_patient_barcode=NA
genes_chip$cdr_type<-NA


tcga_in_chip<-tcga_w_barcode%>%
  dplyr::select(all_of(colnames(genes_chip)))

tcga_and_chip<-rbind(genes_chip,tcga_in_chip)

tcga_and_chip<-tcga_and_chip[c(4:9683),]



tcga_and_chip<-tcga_and_chip[,c(20271,20272,1:20270)]

tcga_and_chip[3:20272]<-sapply(tcga_and_chip[3:20272],as.numeric)

write_fst(tcga_and_chip,"tcga_w_barcode2.fst")

p<-read_fst("tcga_w_barcode2.fst")


