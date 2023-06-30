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

########### UI prep
###make a loading page
appCSS <- "
#loading-content {
  position: fixed;
  top: 50%;
  background: #FFFFFF;
  opacity: 0.9;
  z-index: 100;
  left: 0;
  right: 0;
  height: 100%;
  text-align: center;
  color: #000000;
}
"

#######UI

ui <- fluidPage(
  useShinydashboard(),
  #shinyjs and inlineCSS allow loading page
  useShinyjs(),
  inlineCSS(appCSS),
  
  # Loading message
  div(
    id = "loading-content",
    h2("Loading... This may take a minute or two...")
  ),
  

  
   navbarPage("GS-TCGA",selected ="GS-Corr",
             theme = bs_theme(bootswatch = "sandstone"),
             
             
     ##link to other pages       
             
             tabPanel(HTML("</a></li><li><a href=\"https://gs-tcga.shinyapps.io/home\">Home")),
             
             tabPanel(HTML("</a></li><li><a href=\"https://gs-tcga.shinyapps.io/GS-Surv\">GS-Surv")),
             
             tabPanel(HTML("</a></li><li><a href=\"https://gs-tcga.shinyapps.io/CC-GSEA\">CC-GSEA")),
             
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
             tabPanel(HTML("</a></li><li><a href=\"https://gs-tcga.shinyapps.io/upload\">GS-Surv (Custom)")),
             
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
  
  
  # Hide the loading message when the rest of the server function has executed
  hide(id = "loading-content", anim = TRUE, animType = "fade")    
  show("app-content")
  
}

# Create a Shiny app object ----------------------------------------------------

shinyApp(ui = ui, server = server)

