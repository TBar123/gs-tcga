# Load packages ----------------------------------------------------------------
library(BiocManager)
options(repos = BiocManager::repositories())


library(shiny)
library(ggplot2)
library(data.table)
library(fgsea)
library(dplyr)
library(abind)
library(shinyscreenshot)
library(shinyjs)
library(survival)
library(survminer)
library(bslib)
library(shinydashboard)
library(shinyWidgets)
library(BiocParallel)
library(fst)
library("scales")
library(readr)
library(tidyr)
library(mltools)
library(tidyverse)


options(rsconnect.max.bundle.size=3145728000)
# Load data --------------------------------------------------------------------




##options of groups of pathways for later
pathway_choices<-c("Hallmark","C2","C3","C7","C6")
pathway_choices<-as.data.frame(pathway_choices)

##this is a list of all pathways and what group they go in- provides the options for later
paths<-read.csv("paths.csv")
pathws<-paths$pathway

##this is all the pathways and the genes in them
all_paths<-readRDS("all_paths.RDS")
names(all_paths)<-c("Hallmark","C2","C3","C6","C7")


##gene expression tumour type patient barcode
tcga<-read_fst("tcga_w_barcode2.fst")

##clinical data
cdr_clinicaldata<-read.delim("clinical_data.tsv",sep="\t")

##list of tumour names
tumnames<-tcga%>%
  select(cdr_bcr_patient_barcode,cdr_type)
tum_choices<-as.data.frame(unique(tumnames$cdr_type))
colnames(tum_choices)<-"tum_choices"



#############UI
ui <- fluidPage(
  useShinydashboard(),
  #shinyjs and inlineCSS allow loading page
  useShinyjs(),

  #Custom CSS
  
  ###this is to put the bit under the histogram saying high and low
  tags$footer(tags$style("
                     
                     #footer {
                    display: flex;
                    justify-content: space-between;
                    color:black;
                    height:30px;
                    align-items:center;
                    
               }
                
                     
                     
                     ")),
  
  
  #GS-Surv GS-Surv (Custom) page box
  tags$head(tags$style("
                     
                     #mybox5{height:1300px !important;}
                 
                
                     ")),
  
  #GS-Surv loading box
  tags$head(tags$style("
                     
                     #mybox8{height:1300px !important;}
                 
                
                     ")),
  
  #CC-GSEA CC-GSEA page box
  tags$head(tags$style("
                     
                     #mybox6{height:800px !important;}
                 
                                 ")),
  
  ##GS-Corr GS-Corr page box
  tags$head(tags$style("
                     
                     #mybox7{height:775px !important;}
                 
                
                     ")),
  
  #############################GS-Surv page####################
  navbarPage("GS-TCGA",selected = "GS-Surv",
             theme = bs_theme(bootswatch = "sandstone"),
             ##link to other pages
             tabPanel(HTML("</a></li><li><a href=\"https://gs-tcga.shinyapps.io/home\">Home")),
             
             ##this page
             tabPanel("GS-Surv",
                      sidebarLayout(
                       ##side bar 
                        # Inputs: Select variables to plot
                        sidebarPanel(
                          
                          
                          
                          # Select variable for GSEA dataset
                          selectInput(inputId = "surv_datset", 
                                      label = "List of Gene Sets",
                                      choices = c("Hallmark Gene Sets"="Hallmark","C2: Curated Gene Sets"="C2","C3: Regulatory Target Gene Sets"="C3","C6: Oncogenic Signature Gene Sets"="C6","C7: Immunologic Signature Gene Sets"="C7"),
                                      selected = "Hallmark"),
                          
                          
                          #Select variable for tumour type
                          selectInput(inputId = "ttumour", 
                                      label = "Tumour Type",
                                      choices = c("ACC","GBM","OV","LUAD","PRAD","UCEC","LUSC","BLCA","TGCT","ESCA","PAAD","KIRP","LIHC","CESC","SARC","BRCA","THYM","MESO","COAD","STAD","SKCM","CHOL","KIRC","THCA","HNSC","READ","LGG","DLBC","KICH","UCS","PCPG","UVM"),
                                      selected = "ACC"),
                          
                          # Select variable for gene
                          selectInput(inputId = "gene", 
                                      label = "Gene Set",
                                      choices= NULL),
                          
                          ### select tumours to exclude
                          selectizeInput(inputId = "tumour", 
                                         label = "Select Patient Cases to Exlcude (Optional)",
                                         choices= NULL,
                                         multiple=TRUE),
                          
                          ##go button
                          actionButton("run", label = "Run GS-Surv"),
                          width = "3",
                          
                         ##free text for info 
                          HTML("<p>For more information about tumour type abbreviations please click <a href='https://gdc.cancer.gov/resources-tcga-users/tcga-code-tables/tcga-study-abbreviations'>here</a></p> 
         <p>For information about specifc gene sets please click <a href='https://www.gsea-msigdb.org/gsea/msigdb/genesets.jsp'>here</a></p>
                               <p>For annotations of individual cases please click <a href='https://portal.gdc.cancer.gov/annotations'>here</a></p>")
                        ),
                        
                        # main panels
                        mainPanel(
                          tabsetPanel(
                            
                          
                            
                            tabPanel("Split by thirds",
                                     
                                     br(),
                                     HTML("<p>This tab performs survival analysis of patients split into three equally sized groups based on average gene set expression.</p>"),
                                   
                                     
                                     ##kaplan
                                     plotOutput("third_plot1"),
                                     
                                     ##download survival
                                     downloadButton('downloadPlotsurvthird', 'Download Survival Plot'),
                                     
                                     br(),
                                     br(),
                                     
                                     ##coxPH Ps and HR
                                     fluidRow(
                                       column(
                                         dataTableOutput("coxph_ps"),width = 4
                                       )
                                     ),
                                   
                                     br(),
                                     
                                     ##download 
                                     downloadButton('download_coxph_ps','Download Categorical CoxPH Table'),
                                     br(),
                                     br(),
                                     
                                     fluidRow(
                                       column(
                                         dataTableOutput("cont_coxph_ps"),width = 4
                                       )
                                     ),
                                     br(),
                                     
                                     ##download
                                     downloadButton('download_cont_coxph_ps','Download Continuous CoxPH Table'),
                                     
                                     br(),
                                     br(),
                                     
                                     ##histogram
                                     plotOutput("third_hist"),
                                     ##little thing on the bottom of the histogram saying up and down
                                     HTML('<div class="row">
                                           <div id="footer">
                                           <div >Gene set expression- low</div>
                                           <div ></div>
                                           <div >Gene set expression- high</div>
                                           </div>
                                           </div>
                                           '),
                                     
                                     br(),
                                     
                                     HTML("<p>The histogram shows the median gene expression of each tumour after scaling. The red line shows the cut off points used to separate high, medium and low tumours for survival analysis.</p>"),
                                     
                                     

                                     

                                     ##download hist and tumour info
                                     downloadButton('downloadPlotthirdhist','Download Average Expression Plot'),
                                     downloadButton('Download_tumour_IDs2','Download Tumour Case Data'),

                                     
                                     
                                     ##genes in pathway
                                     fluidRow(
                                       column(
                                         dataTableOutput("Genes3"),width = 4
                                       )
                                     ),
                                     ##download genes in pathway
                                     downloadButton('downloadgenethir','Download Genes in Pathway'),
                                     
                                     
                                     #############
                                     box(
                                       # title = h4("GS-Surv"),
                                       width = 12,
                                       "The survival analysis tool allows the user to see the relationship of a user inputted gene set with survival of patients in the TCGA in a Kaplan-Meier plot. The median expression of genes in the selected gene set and tumour type is calculated, 
               this is then used to generate a Kaplan-Meier plot either splitting data at the median or at the optimal point for discrimination between high and low gene set expression. A histogram of distributions of gene expression and the cut-off is also shown.",
                                       br(),
                                       id = "mybox8",
                                       imageOutput("gs_surv_img4")
                                     )
                                     ############
                                     
                                     
                                     
                            ),
                            #median page
                            tabPanel("Median split",
                                     
                                     br(),
                                     HTML("<p>This tab performs survival analysis of patients split into two groups using a cut-off at the median average gene set expression.</p>"),
                              
                                     
                                     ##kaplan plot
                                     plotOutput("fifty_plot1"),
                                     
                                     
                                     ##download survival
                                     downloadButton('downloadPlotsurvfif', 'Download Survival Plot'),
                                     br(),
                                     br(),
                                     
                                     fluidRow(
                                       column(
                                         dataTableOutput("median_ps"),width = 4
                                       )
                                     ),
                                     br(),
                                     
                                     ##download
                                     downloadButton('download_med_ps','Download Statistics Table'),
                                     
                                     br(),
                                     br(),
                                     
                                     ##histogram median
                                     plotOutput("ave_hist"),
                                     ##a little bit saying high or low beneath
                                     HTML('<div class="row">
                                           <div id="footer">
                                           <div >Gene set expression- low</div>
                                           <div ></div>
                                           <div >Gene set expression- high</div>
                                           </div>
                                           </div>
                                           '),
                                     
                                     br(),
                                     
                                     HTML("<p>The histogram shows the median gene expression of each tumour after scaling. The red line shows the median cut off point used to separate high and low tumours for survival analysis.</p>"),
                                     
                                     
                                     
                                     
                                     ##download hist and tumour info
                                     downloadButton('downloadPlotfifhist','Download Average Expression Plot'),
                                     downloadButton('Download_tumour_IDs','Download Tumour Case Data'),
                                     
                                     
                                     ##genes in pathway table 
                                     fluidRow(
                                       column(
                                         dataTableOutput("Genes2"),width = 4
                                       )
                                     ),
                                     ##download genes in pathway
                                     downloadButton('downloadgenefif'),
                                     
                            ),
                            
                            
                           
                          ),
                          
                        )
                      )
             ),
             ##the other tabs
             #####################CC-GSEA landing page       
             tabPanel(title = "CC-GSEA",
                      
                      h4(strong("CC-GSEA")),
                      p(style="text-align: justify; font-size = 50px",strong("Co-Correlative Gene Set Enrichment Analysis: "),
                        "CC-GSEA allows generation of novel hypotheses of gene function through performing GSEA on co-correlated genes."),
                      p(style="text-align: justify; font-size = 50px",
                        strong("Click the link below to open CC-GSEA. Loading may take a minute or two...")),
                      
                      
                      
                      shiny::actionButton(inputId='ab1', label="Open CC-GSEA", 
                                          icon = icon("th"), 
                                          onclick ="window.location=('https://gs-tcga.shinyapps.io/CC-GSEA/')"),
                      
                      br(),
                      br(),
                      
                      ##CC-GSEA diagram and info                 
                      box(
                        
                        width = 12,
                        "Genes involved in similar cellular processes/pathways may correlate with one another. Co-correlative Gene Set Enrichment Analysis aims to suggest novel gene functions through performing GSEA on a ranked list of genes which correlate with the gene of interest.
               The user can select a tumour type and gene of interest, then the expression of the gene of interest in the TCGA is correlated with all other genes, forming the ranked list used as an input to GSEA. 
               GSEA  analysis outputs a table of gene sets along with their enrichment scores and other statistics. Gene sets with a high normalised enrichment score are enriched in genes that correlate positively with the gene of interest. The user can then select a gene set from this list to explore in enrichment plots and export summary statistics.",
                        id = "mybox6",
                        br(),
                        
                        imageOutput("cc_gsea_img2")
                      )
                      
                      
             ),   
             
             ####################################GS-corr landing page
             
             tabPanel(title = "GS-Corr",
                      
                      h4(strong("GS-Corr")),
                      p(style="text-align: justify; font-size = 50px",strong("Gene Set Correlative Analysis: "),
                        "GS-Corr calculates the average expression of a gene set and correlates this with individual genes."),
                      p(style="text-align: justify; font-size = 50px",
                        strong("Click the link below to open GS-Corr. Loading may take a minute or two...")),
                      
                      
                      
                      shiny::actionButton(inputId='ab1', label="Open GS-Corr", 
                                          icon = icon("th"), 
                                          onclick ="window.location=('https://gs-tcga.shinyapps.io/gene_set_correlates/')"),
                      
                      br(),
                      br(),
                      
                      ##GS-Corr diagram and info      
                      box(
                        
                        width = 12,
                        "The Gene Set Correlates tool aims to suggest novel gene functions through correlating the expression of genes in a pre defined gene set with expression of all other genes. 
               The user can select a tumour type and gene set of interest, then the average expression of genes in the gene set of interest is calculated and correlated with expression of all other genes in the TCGA.",
                        id = "mybox7",
                        br(),
                        
                        imageOutput("gs_corr_img2")
                      )
                      
                      
                      
             ),   
             
             
             
             
             
             
             ###############################GS-Surv custom page 
             
             
             tabPanel(title = "GS-Surv (Custom)",
                      
                      h4(strong("GS-Surv (Custom)")),
                      p(style="text-align: justify; font-size = 50px",strong("GS-Surv (Custom): "),
                        "This function allows you to upload your own gene set for GS-Surv survival analysis. "),
                      p(style="text-align: justify; font-size = 50px",
                        strong("Click the link below to open GS-Surv (Custom). Loading may take a minute or two...")),
                      
                      
                      
                      shiny::actionButton(inputId='ab1', label="Open GS-Surv (Custom)", 
                                          icon = icon("th"), 
                                          onclick ="window.location=('https://gs-tcga.shinyapps.io/upload/')"),
                      
                      br(),
                      br(),
                      
                      ##GS-Surv diagram and info     
                      box(
                        # title = h4("GS-Surv"),
                        width = 12,
                        "The survival analysis tool allows the user to see the relationship of a user inputted gene set with survival of patients in the TCGA in a Kaplan-Meier plot. The median expression of genes in the selected gene set and tumour type is calculated, 
               this is then used to generate a Kaplan-Meier plot either splitting data at the median or at the optimal point for discrimination between high and low gene set expression. A histogram of distributions of gene expression and the cut-off is also shown.",
                        br(),
                        id = "mybox5",
                        imageOutput("gs_surv_img3")
                      )
                      
                      
             ),
             
                      )
             )
  

server <- function(input, output, session) {
  
  
  
  
  # reactive UI-- onyl show gene sets in right group
  path_options_surv <- reactive({
    ## pick out the set of gene sets chosent
    pathwaysel<-pathway_choices%>%
      filter(pathway_choices==input$surv_datset)
    pathwaysel<-as.character(pathwaysel)
    ##filter the list of pathway options to pick out only the set of gene sets chosen
    pathwayz<-paths%>%
      filter(type==pathwaysel)
  })
  ##update the options available for gene set input
  observeEvent(path_options_surv(), {
    choices <- unique(path_options_surv()$pathway)
    updateSelectizeInput(inputId = "gene", choices = choices) 
  })
  
  
  

  ####REV 
  ##get just the selected set of tumour type tumours to show up in the box 
  #reactive UI
  tum_options <- reactive({
    ## pick out the tumour type chosen
    tumsel<-tum_choices%>%
      dplyr::filter(tum_choices==input$ttumour)
    tumsel<-as.character(tumsel)
    ##dplyr::filter the list of tumour options to pick out only the tumour barcodes in the selected tumour type
    tumz<-tumnames%>%
      dplyr::filter(cdr_type==tumsel)
  })
  ##update the options available for input$tumour
  observeEvent(tum_options(), {
    choices_tum <- unique(tum_options()$cdr_bcr_patient_barcode)
    updateSelectInput(inputId = "tumour", choices = choices_tum) 
  })
  
  ##################################################preprocessing
  preprocess<- eventReactive(
    eventExpr = input$run, {
      ##select pathway and genes in pathway
      pathw<-as.data.frame(pathws)
      pathway<-pathw%>%
        filter(pathws==input$gene)
      pathway<-as.character(pathway)
      
      
      to_use<-all_paths[[input$surv_datset]]
      genes_to_use<-to_use[[pathway]]
      genes_to_use<-c("cdr_bcr_patient_barcode",genes_to_use)
      
      
      
      tumours_excluded<-as.character(input$tumour)
      ##pick out the tumour type and remove excluded tumours
      tumour_type<-tcga%>%
        filter(cdr_type==input$ttumour)%>%
        select(any_of(genes_to_use))%>%
        dplyr::filter(!cdr_bcr_patient_barcode %in% tumours_excluded)
      
      tumour_type <- tumour_type[!duplicated(tumour_type$cdr_bcr_patient_barcode), ] 
      ##remove genes with upper quartile less than 10
      tumquar<-apply( tumour_type[,-1] , 2 , quantile , probs = 0.75 , na.rm = TRUE )
      tumquar<-as.data.frame(tumquar)
      tumquar_filt<-tumquar%>%
        filter(tumquar>=10)
      over10<-rownames(tumquar_filt)
      over10<-c("cdr_bcr_patient_barcode",over10)
      ##pick the over 10 genes then scale them between 0 and 1
      tum_ov10<-tumour_type%>%
        select(any_of(over10))%>%
        mutate(across(where(is.numeric),~ rescale(.x)))
      rownames(tum_ov10)<-tum_ov10$cdr_bcr_patient_barcode
      tum_ov10<-tum_ov10[,-1]
      
      tum_ov10<-as.matrix(tum_ov10)
      ##calculate median expression
      mns<-matrixStats::rowMedians(tum_ov10,na.rm=T)
      
      mns<-as.data.frame(mns)
      
      mns$barcode<-rownames(mns)
      
      ##add the clinical on
      means_and_clin<-mns%>%
        left_join(cdr_clinicaldata,by=c("barcode"="bcr_patient_barcode"))
      return(means_and_clin)
  
    })
  
  

  
  
  ###########################################################################################################
  ###########################################################################################################
  

  ##make genes in pathway table
  genetab<-eventReactive(
    eventExpr = input$run, {
    ##select pathway
    pathw<-as.data.frame(pathws)
    pathway<-pathw%>%
      filter(pathws==input$gene)
    pathway<-as.character(pathway)
    
    
    to_use<-all_paths[[input$surv_datset]]
    ##make table
    path_for_table<-to_use[[pathway]]  
    path_for_table<-as.data.frame(path_for_table)
    colnames(path_for_table)<-"Genes in Pathway"
    return(path_for_table)
  })

  
  ##output it
  output$Genes3<-renderDataTable({
    genetab()
  })
  
  ##output it
  output$Genes2<-renderDataTable({
    genetab()
  })
  
  ################################################################################
  ###############################################################################
  ##########5050 plot
  
  ##median cut off categorising
  fifty<-reactive({
    
    
    means_and_clin=preprocess()
    ##put high and low above and below median
    fifty_cut<-means_and_clin%>%
      mutate(Category=ifelse(mns>=median(mns),"high","low"))%>%
      select(OS.time,OS,Category)
    
    colnames(fifty_cut)<-c("OS.time","OS","Category")
 
    
    
    return(fifty_cut)
    
    
  })
  
  
  
  ##kaplan plot it
  output$fifty_plot1<-renderPlot({
    fit=survfit(Surv(OS.time,OS)~Category,data=fifty())
    to_use<-fifty()
    nhigh<-sum(to_use$Category=="high")
    nhigh<-as.numeric(nhigh)
    nlow<-sum(to_use$Category=="low")
    nlow<-as.numeric(nlow)
    legend<-c(paste0("High, n=",nhigh),paste0("Low, n=",nlow))
    ggsurvplot(fit,data=fifty(),xlab="Time (Days)",legend="right",legend.labs=legend,legend.title="Gene Set Expression")
    
    
    
    
    
  })
  ##histogram
  hist_50<-reactive({
    
  
  means_and_clin=preprocess()
  
##hist
  toplot<-ggplot(means_and_clin,aes(x=mns))+
    geom_histogram(bins=50)+
    xlab("Average Pathway Expression")+
    theme_bw()+
    geom_vline(xintercept=median(means_and_clin$mns),linetype="dashed",color="red")
  print(toplot)
  
  
  })
  
  ##output hist
  output$ave_hist<-renderPlot({hist_50()})
  

  
  #################################################################################################
  #########################thirds##################################################################
  third<-reactive({
    
    
    means_and_clin=preprocess()
    
    third_cut<-means_and_clin%>%
      arrange(mns)
    ##put into groups and label
    third_cut<-split(third_cut, cut_number(third_cut$mns, n=3, labels=1:3))
    lowest<-as.data.frame(third_cut[[1]])
    lowest$Category<-"Low"
    middle<-as.data.frame(third_cut[[2]])
    middle$Category<-"Middle"
    high<-as.data.frame(third_cut[[3]])
    high$Category<-"High"
    
    third_cut<-rbind(lowest,middle,high)  
    
    third_cut<-third_cut%>%
      select(OS.time,OS,Category)
    
    colnames(third_cut)<-c("OS.time","OS","Category")

    
    return(third_cut)
    
    
  })
  
  ##kaplan plot it
  output$third_plot1<-renderPlot({
    fit=survfit(Surv(OS.time,OS)~Category,data=third())  
    to_use<-third()
    nhigh<-sum(to_use$Category=="High")
    nhigh<-as.numeric(nhigh)
    nlow<-sum(to_use$Category=="Low")
    nlow<-as.numeric(nlow)
    nmiddle<-sum(to_use$Category=="Middle")
    legend<-c(paste0("High, n=",nhigh),paste0("Low, n=",nlow),paste0("Middle, n=",nmiddle))
    ggsurvplot(fit,data=third(),xlab="Time (Days)",legend="right",legend.labs=legend,legend.title="Gene Set Expression")
    
    
    
    
  })
  
thirds_pvals<-reactive({ 
  
  third_cuts<-third()
  third_cuts$Category<-as.factor(third_cuts$Category)
  third_cuts$Category<-relevel(third_cuts$Category,"Low")
  
  coxfit <- coxph(
    Surv(OS.time,OS) ~ Category,
    data = third_cuts,
    ties = 'exact')
  
 # summary(coxfit)
  
 fit<-summary(coxfit)
 a<-fit$coefficients
 a<-as.data.frame(a)
 a<-a[,c(2,5)]
 
 third_cut2<-third()
 third_cut2$Category<-as.factor(third_cut2$Category)
 third_cut2<-third_cut2%>%
   filter(Category=="Middle"|Category=="High")
 
 third_cut2$Category<-relevel(third_cut2$Category,"Middle")
 
 coxfit2 <- coxph(
   Surv(OS.time,OS) ~ Category,
   data = third_cut2,
   ties = 'exact')
 
 # summary(coxfit)
 
 fit2<-summary(coxfit2)
 a2<-fit2$coefficients
 a2<-as.data.frame(a2)
 a2<-a2[1,c(2,5)]
 
 
 colnames(a)<-c("CoxPH HR","CoxPH p value")
 colnames(a2)<-c("CoxPH HR","CoxPH p value")
 rownames(a)<-c("High vs Low", "Middle vs Low")
 rownames(a2)<-c("High vs Middle")
 combined<-rbind(a,a2)
 
 combined$Comparison<-rownames(combined)
 combined<-combined[,c(3,1,2)]
 return(combined)
  
})  

##output it
output$coxph_ps<-renderDataTable({
  thirds_pvals()
})



  
  
  ##histogram
  hist_30<-reactive({
    
    
    means_and_clin=preprocess()
    
    third_cut<-means_and_clin%>%
      arrange(mns)
    
    third_cut<-split(third_cut, cut_number(third_cut$mns, n=3, labels=1:3))
    lowest<-as.data.frame(third_cut[[1]])
    middle<-as.data.frame(third_cut[[2]])
    
    
    
    ##hist
    toplot<-ggplot(means_and_clin,aes(x=mns))+
      geom_histogram(bins=50)+
      xlab("Average Pathway Expression")+
      theme_bw()+
      geom_vline(xintercept=max(lowest$mns),linetype="dashed",color="red")+
      geom_vline(xintercept=max(middle$mns),linetype="dashed",color="red")
    print(toplot)
    
    
  })
  
  ##output hist
  output$third_hist<-renderPlot({hist_30()})
  
  
  ##################################################################################################
  ##################################################################################################
  #####################################COXPH############################################################################
 ##get coxph Hazard ratio
   cox<-reactive({
    means_and_clin=preprocess()
    
   
    
#####categorical cox
    
  
    coxfit <- coxph(
      Surv(OS.time,OS) ~ Category,
      data = fifty(),
      ties = 'exact')
    
    # summary(coxfit)
    
    fit<-summary(coxfit)
    a<-fit$coefficients
    a<-as.data.frame(a)
    a<-a[,c(2,5)]
    colnames(a)<-c("CoxPH HR - categorical (median split)","CoxPH p value - categorical (median split)")
    a<-a[,c(1,2)]
    
    
    ###categorical logrank
    fifty_cut<-fifty()
    fit=survfit(
    Surv(OS.time,OS)~Category,
    data=fifty_cut) 
    
    fit=survfit(Surv(OS.time,OS)~Category,data=fifty())
    p<-surv_pvalue(fit,fifty_cut)
    p<-p$pval
    p<-as.data.frame(p)
    colnames(p)<-"Logrank p value (median split)"
    
    pvals<-cbind(a,p)
    
    return(pvals)
    

  })
  

  ##output it
  output$median_ps<-renderDataTable({
    cox()
  })
   
  
  
  #######################################################################
 
    cont_cox<-reactive({
    
  means_and_clin=preprocess()
  fit<-coxph(Surv(OS.time, OS)~mns, data=means_and_clin)
  fit<-summary(fit)
  f<-fit$coefficients
  f<-as.data.frame(f)
  f<-f[,c(2,5)]
  colnames(f)<-c("CoxPH HR- median gene set expression as continuous variable","CoxPH p value - median gene set expression as continuous variable")
  f<-f[,c(1,2)]
  
  return(f)
  
})
  
  output$cont_coxph_ps<-renderDataTable({
    cont_cox()
  })
  
  ############################################################
  
  
  
  
  
  
  
  #####################################################################################################################
  
  #######################download tumour info#####################################################################
  tums<-reactive({
    means_and_clin<-preprocess()
    means_and_clin<-means_and_clin%>%
      select(barcode,mns)
    colnames(means_and_clin)<-c("barcode","Median Gene Set Expression")
    return(means_and_clin)
  })
  
  output$Download_tumour_IDs<-downloadHandler(
    filename = "Tumour_data.csv",
    content = function(file) {
      write.csv(tums(), file,row.names=FALSE)
    }
  )
  output$Download_tumour_IDs2<-downloadHandler(
    filename = "Tumour_data.csv",
    content = function(file) {
      write.csv(tums(), file,row.names=FALSE)
    }
  )

  ##############################################################################################################
  

  
  ##download plots
  
  
 
  
  output$download_coxph_ps<-downloadHandler(
    filename = function(){
      paste0(input$ttumour,"_",input$gene,"_CoxPH.csv")
    },
    content = function(file) {
      write.csv(thirds_pvals(), file,row.names=FALSE)
    }
  )
  
  output$download_cont_coxph_ps<-downloadHandler(
    filename = function(){
      paste0(input$ttumour,"_",input$gene,"_continuous_CoxPH.csv")
    },
    content = function(file) {
      write.csv(cont_cox(), file,row.names=FALSE)
    }
  )
  
  output$download_med_ps<-downloadHandler(
    filename = function(){
      paste0(input$ttumour,"_",input$gene,"_statistics.csv")
    },
    content = function(file) {
      write.csv(cox(), file,row.names=FALSE)
    }
  )
  
  output$downloadPlotsurvthird <- downloadHandler(
    filename = function(){
      paste0(input$ttumour,"_",input$gene,"_third_cutoff_survival_plot.svg")
    },
    content = function(file) {
      svg(file)
      print(ggsurvplot((survfit(Surv(OS.time,OS)~Category,data=third())  ),data=third()))
      dev.off()
    }
  )
  
  output$downloadPlotthirdhist <- downloadHandler(
    filename = function(){
      paste0(input$ttumour,"_",input$gene,"_third_cutoff_histogram.svg")
    },
    content = function(file) {
      svg(file)
      print(hist_30())
      dev.off()
    }
  )
  
  

  output$downloadPlotsurvfif <- downloadHandler(
    filename = function(){
      paste0(input$ttumour,"_",input$gene,"_median_cutoff_survival_plot.svg")
    },
    content = function(file) {
      svg(file)
      print(ggsurvplot((survfit(Surv(OS.time,OS)~Category,data=fifty())  ),data=fifty(),pval=TRUE,title =  paste0(input$ttumour," ",input$gene)))
      dev.off()
    }
  )
  
  
  output$downloadPlotfifhist <- downloadHandler(
    filename = function(){
      paste0(input$ttumour,"_",input$gene,"_median_cutoff_histogram.svg")
    },
    content = function(file) {
      svg(file)
      print(hist_50())
      dev.off()
    }
  )
  

  
 

 
  
  
  output$downloadgenefif <- downloadHandler(
    filename = function(){
      paste0(input$gene,"_genes.csv")
    },
    content = function(file) {
      vroom::vroom_write(genetab(), file)
    }
  )
  
  output$downloadgenethir <- downloadHandler(
    filename = function(){
      paste0(input$gene,"_genes.csv")
    },
    content = function(file) {
      vroom::vroom_write(genetab(), file)
    }
  )
  ###############GS-Surv custom page, GS-Surv image
  output$gs_surv_img3 <- renderImage({
    
    list(src = "www/gs_surv_diagram.png",
         height=1200,
         width=800)
    
  }, deleteFile = F)
  ##GS-Surv image for loading
  output$gs_surv_img4 <- renderImage({
    
    list(src = "www/gs_surv_diagram.png",
         height=1200,
         width=800)
    
  }, deleteFile = F)
  
  #########CC-GSEA page CC-GSEA diagram
  
  output$cc_gsea_img2 <- renderImage({
    
    list(src = "www/ccgsea_diagram.png",
         height=700,
         width=750)
    
  }, deleteFile = F)
  
  #############GS-Corr page GS-Corr diagram
  output$gs_corr_img2 <- renderImage({
    
    list(src = "www/GS-Corr_diagram.png",
         height=700,
         width=750)
    
  }, deleteFile = F)
  
  

}

# Create a Shiny app object ----------------------------------------------------

shinyApp(ui = ui, server = server)

