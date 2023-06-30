
require(data.table)
require(stringr)
library(reshape)
library(survival)
library(ggplot2)

rnaseq<-fread('EBPlusPlusAdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp_transposed.tsv', verbose=TRUE, data.table=FALSE, header=TRUE) #Load transposed tcga rnaseq data
colnames(rnaseq)<-sub("\\|.*", "", colnames(rnaseq))
colnames(rnaseq)<-paste0("mrna_",colnames(rnaseq)) #Appends mrna_ before all colnames in table


rnaseq_genes<-colnames(rnaseq[,2:length(rnaseq)]) #Makes a list of all genes (for later column selection)
rnaseq$mrna_gene_id<-substring(rnaseq$mrna_gene_id,1,15) #Makes the TCGA barcode compatible with the methylcs table tcga barcode (cuts off everything to the right of the 15th character which is redundant information)

cdr_clinicaldata<-fread('TCGA-CDR-SupplementalTableS1_sampledata.csv', verbose=TRUE, data.table=FALSE, header=TRUE)
colnames(cdr_clinicaldata)<-paste0("cdr_",colnames(cdr_clinicaldata))  #Appends cdr_ before all colnames in table 

clinicaldata_cols<-colnames(cdr_clinicaldata) # list of clinical data columns for future subsetting/headscratching


rnaseq_primary_tumour<-rnaseq[substring(rnaseq$mrna_gene_id,14,15)=='01',]
rnaseq_primary_tumour$patient_id<-substring(rnaseq_primary_tumour$mrna_gene_id,1,12) # Create a Patient ID column for subsequent merge with clinical data (ie cut off everything from the TCGA barcode to the right of the 12th character - ie we don't need sample type any more)


tcga<-merge(cdr_clinicaldata, rnaseq_primary_tumour, by.x="cdr_bcr_patient_barcode", by.y="patient_id", all.x=FALSE, all.y=FALSE)


names<-colnames(tcga)

names<-gsub("mrna_*","",names)
colnames(tcga)<-names

tcgaa<-tcga[,c(3,65:20566)]
#saveRDS(tcgaa, file = "tcgaa.rds")



write_fst(tcgaa,"tcga.fst")


tcga_w_barcode<-tcga[,c(1,3,65:20566)]

write_fst(tcga_w_barcode,"tcga_w_barcode.fst")