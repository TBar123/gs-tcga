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

# Load data --------------------------------------------------------------------





##the pathways and what group they are in
paths<-read.csv("paths.csv")
pathws<-paths$pathway

##the mean expression of every pathway in every patient
c2<-read_fst("C2_means.fst")
c3<-readRDS("c3_means.rds")
c6<-readRDS("C6_means.rds")
c7<-readRDS("c7_means.rds")
hall<-readRDS("hallmark_means.rds")
##put it together
path_means<-list(hall,c2,c3,c6,c7)
names(path_means)<-c("Hallmark","C2","C3","C6","C7")
##get rid of individual
rm(c2)
rm(c3)
rm(c6)
rm(c7)
rm(hall)

##choices for later
pathway_choices<-c("Hallmark","C2","C3","C7","C6")
pathway_choices<-as.data.frame(pathway_choices)

##gene expression data with barcode
tcga<-read_fst("tcga_w_barcode.fst")

###all the paths and the genes in them 
all_paths<-readRDS("all_paths.RDS")
names(all_paths)<-c("Hallmark","C2","C3","C6","C7")



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
                          
                          
       ##side bar                   
                          #Select variable for GSEA dataset
                          selectInput(inputId = "c", 
                                      label = "GSEA Dataset",
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
                          
                          
                          
                          
                          
                          
                          
                          HTML("<p>For more information about tumour type abbreviations please click <a href='https://gdc.cancer.gov/resources-tcga-users/tcga-code-tables/tcga-study-abbreviations'>here</a></p> 
         <p>For information about specifc gene sets please click <a href='https://www.gsea-msigdb.org/gsea/msigdb/genesets.jsp'>here</a></p>")
                        ),
                        
                        
                        # main panel
                        mainPanel(
                          
                                     HTML("<p>This is the correlation coefficient of the gene set of interest correlated with all genes in the TCGA for the chosen tumour type. The genes in the selected gene set have been removed from the list</p> 
         ")
                                     ,
                                     ##tale
                                     fluidRow(
                                       column(
                                         dataTableOutput("Cor_Table"),width = 4
                                       )
                                     ),
                                     
                                     #download button
                                     downloadButton('Download1'),
                           
                     ##genes in pathway       
                        fluidRow(
                          column(
                            dataTableOutput("Genes"),width = 4
                          )
                        ),
                        downloadButton('downloadgeneopt'),
                                     
                                     
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
  
  
  ###########make the correlation table - j()
  j<-reactive({
    data=tcga
    
    
   # select tumour type for tcga
    tcga_dat<-tcga%>%
      filter(cdr_type==input$z)%>%
      distinct(cdr_bcr_patient_barcode,.keep_all=TRUE)
    
    #select tumour type and gene set of interest from pre processed data
    to_use<-path_means[[input$c]]
    
    gsoi<-to_use%>%
      filter(paths_used==input$y)%>%
      filter(cdr_type==input$z)%>%
      left_join(tcga_dat,by=c("barcode" ="cdr_bcr_patient_barcode"))

    ##correlate
    cor_coef<-cor(gsoi[6:20507], gsoi$mns, method="pearson", use="everything") 
    
   
    ##make sensible table
    cor_coef<-as.data.frame(cor_coef)
    cor_coef$gene<-rownames(cor_coef)
    cor_coef <- cor_coef[, c(2,1)]
    cor_coef<-na.omit(cor_coef)
    colnames(cor_coef)<-c("Gene","Correlation Coefficient")
    cor_coef<-cor_coef %>%
      as_tibble() %>%
      arrange(desc(`Correlation Coefficient`))
    
    ##get the genes in the pathway
    pathw<-as.data.frame(pathws)
    pathway<-pathw%>%
      filter(pathws==input$y)
    pathway<-as.character(pathway)
    
    
    to_usepaths<-all_paths[[input$c]]
    
    path_for_table<-to_usepaths[[pathway]]  

    ##take the genes in the pathway out of the list
    cor_coef<-cor_coef%>%
      filter(!(Gene %in% path_for_table))
    
    
    
    return(cor_coef)
  })
  ##output
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
  genetab<-reactive({
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
##output  
  output$Genes<-renderDataTable({
    genetab()
  })
 ##download 
  output$downloadgeneopt <- downloadHandler(
    filename = function(){
      paste0(input$y,"_genes.csv")
    },
    content = function(file) {
      vroom::vroom_write(genetab(), file)
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
  
}

# Create a Shiny app object ----------------------------------------------------

shinyApp(ui = ui, server = server)

