# Load packages ----------------------------------------------------------------
######options(rsconnect.max.bundle.size=...)

library(BiocManager)
options(repos = BiocManager::repositories())

library(shiny)
library(ggplot2)
library(data.table)
library(dplyr)
library(abind)
library(shinyscreenshot)
library(shinyjs)
library(bslib)
library(shinydashboard)
library(shinyWidgets)
library(fst)
library(ggrepel)
library(tidyr)
library(fgsea)
library(DESeq2)

options(rsconnect.max.bundle.size=3145728000)

# Load data --------------------------------------------------------------------





##the pathways and what group they are in
paths<-read.csv("paths.csv")
pathws<-paths$pathway

##choices for later
pathway_choices<-c("Hallmark","C2","C3","C7","C6")
pathway_choices<-as.data.frame(pathway_choices)

##gene expression data with barcode
tcga<-read_fst("tcga_w_barcode2.fst")

#get just a list of genes
colheads<-colnames(tcga)

inputnames<-colheads[3:length(colheads)]

###all the paths and the genes in them 

all_paths<-readRDS("all_paths.RDS")
names(all_paths)<-c("Hallmark","C2","C3","C6","C7")


## list of tumour barcodes
tumnames<-tcga%>%
  select(cdr_bcr_patient_barcode,cdr_type)
tum_choices<-as.data.frame(unique(tumnames$cdr_type))
colnames(tum_choices)<-"tum_choices"


#######UI

ui <- fluidPage(
  useShinydashboard(),
  #shinyjs and inlineCSS allow loading page
  useShinyjs(),
  
  
  #GS-Surv GS-Surv page box
  tags$head(tags$style("
                     
                     #mybox4{height:1300px !important;}
                 
                
                     ")),
  
  #GS-Surv GS-Surv (Custom) page box
  tags$head(tags$style("
                     
                     #mybox5{height:1300px !important;}
                 
                
                     ")),
  #CC-GSEA CC-GSEA page box
  tags$head(tags$style("
                     
                     #mybox6{height:800px !important;}
                 
                                 ")),
  
  
  ##GS-Corr GS-Corr page box
  tags$head(tags$style("
                     
                     #mybox7{height:775px !important;}
                 
                
                     ")),

  
   navbarPage("GS-TCGA",selected ="GS-Corr",
             theme = bs_theme(bootswatch = "sandstone"),
             
             
     ##link to other pages       
             
             tabPanel(HTML("</a></li><li><a href=\"https://gs-tcga.shinyapps.io/home\">Home")),
             
     ########################GS-Surv landing page
     tabPanel(title = "GS-Surv",
              
              h4(strong("GS-Surv")),
              p(style="text-align: justify; font-size = 50px",strong("Gene Set Survival Analysis: "),
                "GS-Surv allows the user to investigate how the average expression of genes in a specified gene set relates to overall survival in patient data. "),
              p(style="text-align: justify; font-size = 50px",
                strong("Click the link below to open GS-Surv. Loading may take a minute or two...")),
              
              
              
              shiny::actionButton(inputId='ab1', label="Open GS-Surv", 
                                  icon = icon("th"), 
                                  onclick ="window.location=('https://gs-tcga.shinyapps.io/GS-Surv/')"),
              
              br(),
              br(),
              
              ##GS-Surv diagram and info     
              box(
                # title = h4("GS-Surv"),
                width = 12,
                "The survival analysis tool allows the user to see the relationship of their selected gene set with survival of patients in the TCGA in a Kaplan-Meier plot. The mean expression of genes in the selected gene set and tumour type is calculated, 
               this is then used to generate a Kaplan-Meier plot either splitting data at the median or at the optimal point for discrimination between high and low gene set expression. A histogram of distributions of gene expression and the cut-off is also shown. Users can upload custom gene sets for analysis on the GS-Surv (Custom) page.",
                br(),
                id = "mybox4",
                imageOutput("gs_surv_img2")
              )
              
              
     ),
     
     
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
     
             
       ##the GS-Corr page      
             tabPanel("GS-Corr",
                      sidebarLayout(
                        
                        # Inputs: Select variables to plot
                        sidebarPanel(
                          ###GS-corr side panel
                          conditionalPanel(condition="input.conditionedPanels == 'GS-Corr'",  
                        
                          #Select variable for GSEA dataset
                          selectInput(inputId = "c", 
                                      label = "List of Gene Sets",
                                      choices = c("Hallmark Gene Sets"="Hallmark","C2: Curated Gene Sets"="C2","C3: Regulatory Target Gene Sets"="C3","C6: Oncogenic Signature Gene Sets"="C6","C7: Immunologic Signature Gene Sets"="C7"),
                                      selected = "Hallmark"),
                          
                          
                          # Select variable for gene set
                          selectInput(inputId = "y", 
                                      label = "Gene Set",
                                      choices= NULL),
                          
                          #Select variable for tumour type
                          selectInput(inputId = "z", 
                                      label = "Tumour Type",
                                      choices = c("UVM","ACC","GBM","OV","LUAD","PRAD","UCEC","LUSC","BLCA","TGCT","ESCA","PAAD","KIRP","LIHC","CESC","SARC","BRCA","THYM","MESO","COAD","STAD","SKCM","CHOL","KIRC","THCA","HNSC","READ","LGG","DLBC","KICH","UCS","PCPG"),
                                      selected = "UVM"),
                          
                          ## 
                          ###select correlation method
                          selectInput(inputId = "corr", 
                                      label = "Correlation method",
                                      choices = c("Pearson","Spearman"),
                                      selected = "Pearson"),
                          
                          
                          ###- select tumours to exclude
                          selectizeInput(inputId = "tumour", 
                                         label = "Select Patient Cases to Exlcude (Optional)",
                                         choices= NULL,
                                         multiple=TRUE),
                          
                          ##make a go button
                          actionButton("run", label = "Run GS-Corr"),
                          width = "3",
                          
                          
                          
                          
                          
                          
                          HTML("<p>For more information about tumour type abbreviations please click <a href='https://gdc.cancer.gov/resources-tcga-users/tcga-code-tables/tcga-study-abbreviations'>here</a></p> 
         <p>For information about specifc gene sets please click <a href='https://www.gsea-msigdb.org/gsea/msigdb/genesets.jsp'>here</a></p>
                               <p>For annotations of individual cases please click <a href='https://portal.gdc.cancer.gov/annotations'>here</a></p>")
                        ),
                        
                      
                        
                        ###GS-Custom side panel
                        conditionalPanel(condition="input.conditionedPanels == 'GS-Corr (Custom)'",
                                         
                                         ##select if you want a file or a copy paste box
                                         selectInput("choose", "Choose Input", choices = c("File", "Text"), selected = NULL),
                                         
                                         HTML("<p>Please enter human gene symbols.</p>"),
                                         
                                         
                                         ##file upload
                                         conditionalPanel(
                                           "input.choose=='File'",
                                           fileInput("file1", "Choose CSV file", multiple = TRUE, accept = c(
                                             "text/csv",
                                             "text/comma-separated-values,text/plain",
                                             ".csv")),
                                           tags$hr(),
                                           checkboxInput("header", "Header", TRUE),
                                           
                                           
                                         ),
                                         ##text upload
                                         conditionalPanel(
                                           "input.choose=='Text'",
                                           textAreaInput("caption", "List of Genes", "Please enter each gene on a separate row", width = "1000px")
                                         ),
                                         
                                         
                                         #Select variable for tumour type
                                         selectInput(inputId = "z2", 
                                                     label = "Tumour Type",
                                                     choices = c("UVM","ACC","GBM","OV","LUAD","PRAD","UCEC","LUSC","BLCA","TGCT","ESCA","PAAD","KIRP","LIHC","CESC","SARC","BRCA","THYM","MESO","COAD","STAD","SKCM","CHOL","KIRC","THCA","HNSC","READ","LGG","DLBC","KICH","UCS","PCPG"),
                                                     selected = "UVM"),
                                         
                                         ## 
                                         ###select correlation method
                                         selectInput(inputId = "corr2", 
                                                     label = "Correlation method",
                                                     choices = c("Pearson","Spearman"),
                                                     selected = "Pearson"),
                                         
                                         
                                         ### - select tumours to exclude
                                         selectizeInput(inputId = "tumour2", 
                                                        label = "Select Patient Cases to Exlcude (Optional)",
                                                        choices= NULL,
                                                        multiple=TRUE),
                                         
                                         ##make a go button
                                         actionButton("run_custom", label = "Run GS-Corr (Custom)"),
                                         width = "3",
                                         
                                         
                                         
                                         
                                         
                                         
                                         HTML("<p>For more information about tumour type abbreviations please click <a href='https://gdc.cancer.gov/resources-tcga-users/tcga-code-tables/tcga-study-abbreviations'>here</a></p> 
         <p>For information about specifc gene sets please click <a href='https://www.gsea-msigdb.org/gsea/msigdb/genesets.jsp'>here</a></p>
                               <p>For annotations of individual cases please click <a href='https://portal.gdc.cancer.gov/annotations'>here</a></p>")
                        
              
                        ),
                        ###Gene set membership side panel
                        conditionalPanel(condition="input.conditionedPanels == 'Gene Set Membership'",
                                         # Select gene
                                         selectizeInput(inputId = "x", 
                                                        label = "Gene",
                                                        choices=inputnames),  
                                         
                                         #Select variable for GSEA dataset
                                         selectInput(inputId = "c2", 
                                                     label = "List of Gene Sets",
                                                     choices = c("Hallmark Gene Sets"="Hallmark","C2: Curated Gene Sets"="C2","C3: Regulatory Target Gene Sets"="C3","C6: Oncogenic Signature Gene Sets"="C6","C7: Immunologic Signature Gene Sets"="C7"),
                                                     selected = "Hallmark"),
                                         
                                         ##make a go button
                                         actionButton("run_membership", label = "Run Gene Set Membership"),
                                         width = "3",
                        ),
                        ),
                        # main panel
                        mainPanel(
                          
                          tabsetPanel(
                            
                            
                            
                            #correlation table
                            tabPanel("GS-Corr",
                                     br(),
                                     
                                     HTML("<p>This is the correlation coefficient of the gene set of interest correlated with all genes in the TCGA for the chosen tumour type. The genes in the selected gene set have been removed from the list.</p> 
         ")
                                     ,
                                     ##table
                                     fluidRow(
                                       column(
                                         dataTableOutput("Cor_Table"),width = 4
                                       )
                                     ),
                                     
                                     #download button- gs corr table and tumour IDs used
                                     br(),
                                     downloadButton('Download1','Download GS Corr Table'),
                                     downloadButton('Download_tumour_IDs','Download Tumour IDs'),
                                     br(),
                                     br(),
                           
                     ##genes in pathway       
                        fluidRow(
                          column(
                            dataTableOutput("Genes"),width = 4
                          )
                        ),
                     br(),
                     ##download genes in gene set
                        downloadButton('downloadgeneopt','Download genes in gene set'),
                             
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
                     
                     ##############custom GS-Corr panel
                     tabPanel("GS-Corr (Custom)",
                              
                              HTML("<p>This is the correlation coefficient of the user uploaded gene set correlated with all genes in the TCGA for the chosen tumour type. The genes in the user inputted gene set have been removed from the list.</p> 
         ")
                              ,
                              ##table
                              fluidRow(
                                column(
                                  dataTableOutput("Cor_Table2"),width = 4
                                )
                              ),
                              br(),
                              #download button- table and tumours analysed
                              downloadButton('Download12','Download GS Corr Table'),
                              downloadButton('Download_tumour_IDs2','Download Tumour IDs'),
                              
                              
                              
                     ),
                     ############make a panel of 'gene set membership' - for any given gene what gene sets is it in?
                     tabPanel("Gene Set Membership",
                              HTML("<p>This tab shows which gene sets contain the selected gene. If no data is available in the table then this gene is not in any gene sets.</p> 
         ")
                              ,
                              
                              ###table of genes sets containing gene of interest
                              fluidRow(
                                column(
                                  dataTableOutput("Gene_set_membership"),width = 4
                                )
                              ),
                              ##download
                              downloadButton('Download_membership','Download Table'),
                              
                     ),
                     
                     id = "conditionedPanels"  
                          ),  
                          ),       
                            
                            
                            
                          
                        
                      )
             ),
     
     ##more pages
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
                "The survival analysis tool allows the user to see the relationship of a user inputted gene set with survival of patients in the TCGA in a Kaplan-Meier plot. The mean expression of genes in the selected gene set and tumour type is calculated, 
               this is then used to generate a Kaplan-Meier plot either splitting data at the median or at the optimal point for discrimination between high and low gene set expression. A histogram of distributions of gene expression and the cut-off is also shown.",
                br(),
                id = "mybox5",
                imageOutput("gs_surv_img3")
              )
              
              
     ),
     
             
  )
)
  
  
  


server <- function(input, output, session) {
  
  

  #reactive UI- make only the set of pathways chosen show up in the gene set box
  path_options <- reactive({
    pathwaysel<-pathway_choices%>%
      filter(pathway_choices==input$c)
    pathwaysel<-as.character(pathwaysel)
    pathwayz<-paths%>%
      filter(type==pathwaysel)
  })
  observeEvent(path_options(), {
    choices <- unique(path_options()$pathway)
    updateSelectInput(inputId = "y", choices = choices) 
  })
  

  ##get just the selected set of tumour type tumours to show up in the box 
  #reactive UI
  tum_options <- reactive({
    ## pick out the selected tumour type
    tumsel<-tum_choices%>%
      dplyr::filter(tum_choices==input$z)
    tumsel<-as.character(tumsel)
    ##dplyr::filter the list of tumours to have only selected tumour type
    tumz<-tumnames%>%
      dplyr::filter(cdr_type==tumsel)
  })
  ##update the options available for input$tumour
  observeEvent(tum_options(), {
    choices_tum <- unique(tum_options()$cdr_bcr_patient_barcode)
    updateSelectInput(inputId = "tumour", choices = choices_tum) 
  })
  
  ##########################################################################################################################################
 
  #pre processing and GS-Corr table----------------------------------

  
  j<- eventReactive(
    eventExpr = input$run, {
    data=tcga
    
    
    # select tumour type for tcga
    tcga_dat<-tcga%>%
      filter(cdr_type==input$z)%>%
      distinct(cdr_bcr_patient_barcode,.keep_all=TRUE)
    
    ## get rid of excluded tumours
    tumours_excluded<-as.character(input$tumour)
    tcga_dat<-tcga_dat%>%
      dplyr::filter(!cdr_bcr_patient_barcode %in% tumours_excluded)
    
    
    ##select gene set
    pathw<-as.data.frame(pathws)
    pathway<-pathw%>%
      filter(pathws==input$y)
    pathway<-as.character(pathway)
    
    
    to_use<-all_paths[[input$c]]
     gsoi<-to_use[[pathway]]
    
    
    
    ###select genes with UQ>10
    to_analyse<-tcga_dat%>%
       select(cdr_bcr_patient_barcode,any_of(gsoi))
    tum<-to_analyse[,c(-1)]
    tumquar<-apply(tum , 2 , quantile , probs = 0.75 , na.rm = TRUE )
    tumquar<-as.data.frame(tumquar)
    tumquar_filt<-tumquar%>%
      filter(tumquar>=10)
    genes<-rownames(tumquar_filt)
    genes<-c("cdr_bcr_patient_barcode",genes)
    genes_selected<-to_analyse[names(to_analyse) %in% genes]
    
    rm(tumquar)
    rm(tumquar_filt)
    rm(tum)
      
    path_selected_genes<-genes_selected%>%
      ##scale to 0 to 1- expression of each gene
      mutate(across(where(is.numeric),~scales::rescale(.x)))
    
    rownames(path_selected_genes)<-path_selected_genes$cdr_bcr_patient_barcode
    path_selected_genes<-path_selected_genes%>%
      select(-cdr_bcr_patient_barcode)

    ###calc median gene set exp
    path_selected_genes<-as.matrix(path_selected_genes)
    meds<-matrixStats::rowMedians(path_selected_genes,na.rm=T)
    meds<-as.data.frame(meds)
    meds$cdr_bcr_patient_barcode<-rownames(meds)
    
    rm(path_selected_genes)
    
    all_to_corr<-merge(tcga_dat,meds,by="cdr_bcr_patient_barcode")
    
    rm(tcga_dat)
    rm(meds)
    
    ##correlate
    cor_coef<-cor(all_to_corr[3:20272], all_to_corr$meds, method=tolower(input$corr), use="everything") 
    
    ##make sensible table
    cor_coef<-as.data.frame(cor_coef)
    cor_coef$gene<-rownames(cor_coef)
    cor_coef <- cor_coef[, c(2,1)]
    cor_coef<-na.omit(cor_coef)
    colnames(cor_coef)<-c("Gene","Correlation Coefficient")
    cor_coef<-cor_coef %>%
      as_tibble() %>%
      arrange(desc(`Correlation Coefficient`))
    
   
    
    ##take the genes in the pathway out of the list
    cor_coef<-cor_coef%>%
      filter(!(Gene %in% gsoi))
    
    
    
    return(cor_coef)
   
  })
  
  
  
  
  
  
  
  
  
  ##output GS-Corr table
  output$Cor_Table <-renderDataTable({
    j()
  })
  

  ##download
  output$Download1 <- downloadHandler(
    filename = function(){
      paste0(input$y,"_",input$z,"_correlation_table.csv")
    },
    content = function(file) {
      write.csv(j(),file,row.names=FALSE)
    }
  )
 
  
  ##genes in pathway
  genetab<- eventReactive(
    eventExpr = input$run, {
    
    
    #reactive({
 ##select pathway
       pathw<-as.data.frame(pathws)
    pathway<-pathw%>%
      filter(pathws==input$y)
    pathway<-as.character(pathway)
    
    
    to_use<-all_paths[[input$c]]
    ##make table
    path_for_table<-to_use[[pathway]]  
    path_for_table<-as.data.frame(path_for_table)
    colnames(path_for_table)<-"Genes in Pathway"
    return(path_for_table)
  })
##output  genes in pathway
  output$Genes<-renderDataTable({
    genetab()
  })
 ##download genes in pathway
  output$downloadgeneopt <- downloadHandler(
    filename = function(){
      paste0(input$y,"_genes.csv")
    },
    content = function(file) {
      vroom::vroom_write(genetab(), file)
    }
  )
  
  
  
  
  ############################################################################################################
  ###########################################################################################################
  ####################### gene set membership 
  
  gsmem<-eventReactive(
    eventExpr = input$run_membership, {
      gene<-as.character(input$x)
      paths<-all_paths[[input$c2]]
      ##put name in list if its in pathway, otherwise say not present
      i<-1
      paths_to_fill<-list()
      for(i in 1:length(paths)){
        pathw<-paths[[i]]
        if(gene %in% pathw){paths_to_fill[[i]]<-names(paths)[[i]]} else{paths_to_fill[[i]]<-"not present"}
      }
      
      paths_to_fill<-as.data.frame(unlist(paths_to_fill))
      colnames(paths_to_fill)<-"Gene sets containing gene"
      paths_to_fill<-paths_to_fill%>%
        filter(!`Gene sets containing gene`=="not present")
      
      return(paths_to_fill)
      
    })

  
  ##print table
  output$Gene_set_membership <-renderDataTable({
    gsmem()
  })
####download gene set membership
  output$Download_membership<-downloadHandler(
    filename = function(){
      paste0(input$x,"_",input$c2,"_gene_set_membership.csv")
    },
    content = function(file) {
      write.csv(gsmem(), file,row.names=FALSE)
    }
  )
  
  
  ####download tumours IDs
  tummys<-reactive({
    tumours_excluded<-as.character(input$tumour)
    tumour_type<-tcga%>%
      dplyr::filter(cdr_type==input$z)%>%
      dplyr::filter(!cdr_bcr_patient_barcode %in% tumours_excluded)%>%
      select(cdr_bcr_patient_barcode)
    return(tumour_type)
  })
  
  output$Download_tumour_IDs<-downloadHandler(
    filename = "Tumours_IDs_in_analysis.csv",
    content = function(file) {
      write.csv(tummys(), file,row.names=FALSE)
    }
  )
  
  
  
  
  
  ######################################################################
  ############################ custom
  #####Get the uploaded gene set
  data <- reactiveVal()
  ##make the action button work and do the right things for text and files
  observeEvent(
    eventExpr = input$run_custom,
    handlerExpr = {
      switch(input$choose,
             File = read.csv(input$file1$datapath),
             Text = data.frame(Genes = c(input$caption))%>%
               separate_longer_delim(Genes,delim="\n")#delim=" "
      ) %>% data()
    }
  ) 
  

  
  ####REV 
  ##get just the selected set of tumour type tumours to show up in the box 
  tum_options2 <- reactive({
    ## pick out the set of tumour barcodes selected
    tumsel2<-tum_choices%>%
      dplyr::filter(tum_choices==input$z2)
    tumsel2<-as.character(tumsel2)
    ##dplyr::filter to show only selected tumour type in box
    tumz2<-tumnames%>%
      dplyr::filter(cdr_type==tumsel2)
  })
  ##update the options available for input$tumour2
  observeEvent(tum_options2(), {
    choices_tum2 <- unique(tum_options2()$cdr_bcr_patient_barcode)
    updateSelectInput(inputId = "tumour2", choices = choices_tum2) 
  })
  
  
  #####pre process and make GS-Corr custom table
  j2<- eventReactive(
    eventExpr = input$run_custom, {
      data=tcga
      
      
      # select tumour type for tcga
      tcga_dat<-tcga%>%
        filter(cdr_type==input$z2)%>%
        distinct(cdr_bcr_patient_barcode,.keep_all=TRUE)
      
      ## get rid of excluded tumours
      tumours_excluded<-as.character(input$tumour2)
      tcga_dat<-tcga_dat%>%
        dplyr::filter(!cdr_bcr_patient_barcode %in% tumours_excluded)
      
      
      ##select gene set
      data<-data()
      colnames(data)<-"Genes" ##
      data$Genes<-toupper(data$Genes)
      gsoi<-data$Genes
      gsoi<-as.character(gsoi)
      

      ###select genes with UQ>10
      to_analyse<-tcga_dat%>%
        select(cdr_bcr_patient_barcode,any_of(gsoi))
      tum<-to_analyse[,c(-1)]
      tumquar<-apply(tum , 2 , quantile , probs = 0.75 , na.rm = TRUE )
      tumquar<-as.data.frame(tumquar)
      tumquar_filt<-tumquar%>%
        filter(tumquar>=10)
      genes<-rownames(tumquar_filt)
      genes<-c("cdr_bcr_patient_barcode",genes)
      genes_selected<-to_analyse[names(to_analyse) %in% genes]
      
      rm(tumquar)
      rm(tumquar_filt)
      rm(tum)
      
      
      path_selected_genes<-genes_selected%>%
        ##scale to 0 to 1- expression of each gene
        mutate(across(where(is.numeric),~scales::rescale(.x)))
      
      rownames(path_selected_genes)<-path_selected_genes$cdr_bcr_patient_barcode
      path_selected_genes<-path_selected_genes%>%
        select(-cdr_bcr_patient_barcode)
      
      ###do median gene set exp
      path_selected_genes<-as.matrix(path_selected_genes)
      meds<-matrixStats::rowMedians(path_selected_genes,na.rm=T)
      meds<-as.data.frame(meds)
      meds$cdr_bcr_patient_barcode<-rownames(meds)
      
      rm(path_selected_genes)
      
      all_to_corr<-merge(tcga_dat,meds,by="cdr_bcr_patient_barcode")
      
      rm(tcga_dat)
      rm(meds)
      
      ##correlate
      cor_coef<-cor(all_to_corr[3:20272], all_to_corr$meds, method=tolower(input$corr), use="everything") 
      
      ##make sensible table
      cor_coef<-as.data.frame(cor_coef)
      cor_coef$gene<-rownames(cor_coef)
      cor_coef <- cor_coef[, c(2,1)]
      cor_coef<-na.omit(cor_coef)
      colnames(cor_coef)<-c("Gene","Correlation Coefficient")
      cor_coef<-cor_coef %>%
        as_tibble() %>%
        arrange(desc(`Correlation Coefficient`))
      
      
      
      ##take the genes in the pathway out of the list
      cor_coef<-cor_coef%>%
        filter(!(Gene %in% gsoi))
      
      
      
      return(cor_coef)
      
    })
  
  
  
  
  
  ##output custom table
  output$Cor_Table2 <-renderDataTable({
    j2()
  })
  
  
  ##download custom
  output$Download12 <- downloadHandler(
    filename = function(){
      paste0(input$z2,"_custom_correlation_table.csv")
    },
    content = function(file) {
      write.csv(j2(),file,row.names=FALSE)
    }
  )
  
  
  
  ####download tumours IDs
  tummys2<-reactive({
    tumours_excluded<-as.character(input$tumour2)
    tumour_type<-tcga%>%
      dplyr::filter(cdr_type==input$z2)%>%
      dplyr::filter(!cdr_bcr_patient_barcode %in% tumours_excluded)%>%
      select(cdr_bcr_patient_barcode)
    return(tumour_type)
  })
  
  output$Download_tumour_IDs2<-downloadHandler(
    filename = "Tumours_IDs_in_analysis.csv",
    content = function(file) {
      write.csv(tummys2(), file,row.names=FALSE)
    }
  ) 
  
  
  
  
  ###############GS-Surv page, GS-Surv image
  output$gs_surv_img2 <- renderImage({
    
    list(src = "www/gs_surv_diagram.png",
         height=1200,
         width=800)
    
  }, deleteFile = F)
  
  ###############GS-Surv custom page, GS-Surv image
  output$gs_surv_img3 <- renderImage({
    
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
  
  
  ###GS-Corr image
  
  #############GS-Corr page GS-Corr diagram
  output$gs_corr_img2 <- renderImage({
    
    list(src = "www/GS-Corr_diagram.png",
         height=700,
         width=750)
    
  }, deleteFile = F)
  
}

# Create a Shiny app object ----------------------------------------------------

shinyApp(ui = ui, server = server)

