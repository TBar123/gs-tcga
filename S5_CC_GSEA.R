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
library(ggrepel)
library(readr)
library(tidyr)
library("scales")


options(rsconnect.max.bundle.size=3145728000)

# Load data --------------------------------------------------------------------



##load in the pathway data- paths and the genes in it
all_paths<-readRDS("all_paths.RDS")
names(all_paths)<-c("Hallmark","C2","C3","C6","C7")


##this is a file of all gene sets and which set of gene sets they belong to
paths<-read.csv("paths.csv")
pathws<-paths$pathway



##this is for selection later
pathway_choices<-c("Hallmark","C2","C3","C7","C6")
pathway_choices<-as.data.frame(pathway_choices)

##this is all gene expression plus the tumour type
tcga<-read_fst("tcga_w_barcode2.fst")
#get just a list of genes
colheads<-colnames(tcga)
inputnames<-colheads[3:length(colheads)]

##list of tumour names
tumnames<-tcga%>%
  select(cdr_bcr_patient_barcode,cdr_type)
tum_choices<-as.data.frame(unique(tumnames$cdr_type))
colnames(tum_choices)<-"tum_choices"





#######UI

ui <- fluidPage(
  useShinydashboard(),
  #shinyjs and inlineCSS allow loading page
  useShinyjs(),

  #GS-Surv GS-Surv (Custom) page box
  tags$head(tags$style("
                     
                     #mybox5{height:1300px !important;}
                 
                
                     ")),
 
  
  ##GS-Corr GS-Corr page box
  tags$head(tags$style("
                     
                     #mybox7{height:775px !important;}
                 
                
                     ")),
  
  
  ##CCGSEA home page box
  tags$head(tags$style("
                     
                     #mybox{height:900px !important;}
                 
                                 ")),
 
  

  ##set up the website with the current page selected to start
  navbarPage("GS-TCGA",selected ="CC-GSEA",
             ##nice colours
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
             
##CCSEA page   


             tabPanel("CC-GSEA",
                      
                       
                      
                      ##make the side bar
                      sidebarLayout(
                        
                        # Inputs: Select variables to plot
                        ##side bar of correlation table page
                        sidebarPanel(
                          conditionalPanel(condition="input.conditionedPanels == 'Correlation Table'",       
                                          
                        
                          
                          
                          
                          #Select variable for GSEA dataset
                          selectInput(inputId = "c", 
                                      label = "List of Gene Sets",
                                      choices = c("Hallmark Gene Sets"="Hallmark","C2: Curated Gene Sets"="C2","C3: Regulatory Target Gene Sets"="C3","C6: Oncogenic Signature Gene Sets"="C6","C7: Immunologic Signature Gene Sets"="C7"),
                                      selected = "Hallmark"),
                          
                          
                          # Select gene
                          selectizeInput(inputId = "x", 
                                         label = "Gene",
                                         choices=inputnames),
                          
                          #Select variable for tumour type
                          selectInput(inputId = "z", 
                                      label = "Tumour Type",
                                      choices = c("UVM","ACC","GBM","OV","LUAD","PRAD","UCEC","LUSC","BLCA","TGCT","ESCA","PAAD","KIRP","LIHC","CESC","SARC","BRCA","THYM","MESO","COAD","STAD","SKCM","CHOL","KIRC","THCA","HNSC","READ","LGG","DLBC","KICH","UCS","PCPG"),
                                      selected = "UVM"),
                          
                          
                          ###select correlation method
                          selectInput(inputId = "corr", 
                                      label = "Correlation method",
                                      choices = c("Pearson","Spearman"),
                                      selected = "Pearson"),
                          
                          
                          ### select tumours to exclude
                          selectizeInput(inputId = "tumour", 
                                         label = "Select Patient Cases to Exlcude (Optional)",
                                         choices= NULL,
                                         multiple=TRUE),
                          
                          
                          ##make a go button
                          actionButton("run", label = "Run CC-GSEA"),
                          width = "3",
                          
                          
                     ###free text to add extra info     
                          HTML("<p>For more information about tumour type abbreviations please click <a href='https://gdc.cancer.gov/resources-tcga-users/tcga-code-tables/tcga-study-abbreviations'>here</a></p> 
         <p>For information about specifc gene sets please click <a href='https://www.gsea-msigdb.org/gsea/msigdb/genesets.jsp'>here</a></p>
              <p>For annotations of individual cases please click <a href='https://portal.gdc.cancer.gov/annotations'>here</a></p>")
                        ),
                       
                     ###side bar for custom page
                     conditionalPanel(condition="input.conditionedPanels == 'CC-GSEA (Custom)'",       
                                      
                                      
                                      
      
                                      # Select gene
                                      selectizeInput(inputId = "x2", 
                                                     label = "Gene",
                                                     choices=inputnames),
                                      
                                      #Select variable for tumour type
                                      selectInput(inputId = "tum_type", 
                                                  label = "Tumour Type",
                                                  choices = c("UVM","ACC","GBM","OV","LUAD","PRAD","UCEC","LUSC","BLCA","TGCT","ESCA","PAAD","KIRP","LIHC","CESC","SARC","BRCA","THYM","MESO","COAD","STAD","SKCM","CHOL","KIRC","THCA","HNSC","READ","LGG","DLBC","KICH","UCS","PCPG"),
                                                  selected = "UVM"),
                                      
                                      ## REV
                                      ###select correlation method
                                      selectInput(inputId = "corr2", 
                                                  label = "Correlation method",
                                                  choices = c("Pearson","Spearman"),
                                                  selected = "Pearson"),
                                      
                                      
                                      ###REV - select tumours to exclude
                                      selectizeInput(inputId = "tumour_excl", 
                                                     label = "Select Patient Cases to Exlcude (Optional)",
                                                     choices= NULL,
                                                     multiple=TRUE),
                                      
                                      
                                      
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
                                      
                                      
                                      
                                      ##make a go button
                                      actionButton("run_custom", label = "Run CC-GSEA (Custom)"),
                                      width = "3",
                                      
                                      ###free text to add extra info     
                                      HTML("<p>For more information about tumour type abbreviations please click <a href='https://gdc.cancer.gov/resources-tcga-users/tcga-code-tables/tcga-study-abbreviations'>here</a></p> 
         <p>For information about specifc gene sets please click <a href='https://www.gsea-msigdb.org/gsea/msigdb/genesets.jsp'>here</a></p>
              <p>For annotations of individual cases please click <a href='https://portal.gdc.cancer.gov/annotations'>here</a></p>")
                     ),
                     
                      ),    
                       
                     
                     ##make the results area 
                        mainPanel(
                          
                          ##results page 1- correlation table
                          tabsetPanel(
                            
                            
                            
                            #correlation table
                            tabPanel("Correlation Table",
                                     br(),
                                     HTML("<p>This is the correlation coefficient of the gene of interest correlated with all other genes in the TCGA for the chosen tumour type.
                                     This is used as the input for GSEA analysis in the CC-GSEA tab.</p> 
         "),
                                     HTML("<p> Click the <b>CC-GSEA tab</b> to perform CC-GSEA analysis of your chosen gene. This analysis may take a while, please be patient... </p> 
         "),
                                     
                                     fluidRow(
                                       column(
                                         dataTableOutput("Cor_Table"),width = 4
                                       )
                                     ),
                                     
                                     #download button- corr table and tumour IDs used
                                     downloadButton('Download1'),
                                     downloadButton('Download_tumour_IDs','Download Tumour IDs'),
                                     
                                     br(),
                                     box(
                                       
                                       width = 12,
                                       "Genes involved in similar cellular processes/pathways may correlate with one another. Co-correlative Gene Set Enrichment Analysis aims to suggest novel gene functions through performing GSEA on a ranked list of genes which correlate with the gene of interest.
               The user can select a tumour type and gene of interest, then the expression of the gene of interest in the TCGA is correlated with all other genes, forming the ranked list used as an input to GSEA. 
               GSEA  analysis outputs a table of gene sets along with their enrichment scores and other statistics. Gene sets with a high normalised enrichment score are enriched in genes that correlate positively with the gene of interest. The user can then select a gene set from this list to explore in enrichment plots and export summary statistics.",
                                       id = "mybox",
                                       br(),
                                       
                                       imageOutput("cc_gsea_img2")
                                     )
                            ),
                            
                            
                            
                            ##results page 2
                            #CC GSEA table
                            tabPanel("CC-GSEA",
                                     br(),
                                     HTML("<p> A rank ordered list of the correlation coefficients of the gene of interest correlated with all other genes in the chosen tumour type was used as the input for GSEA analyses. This provides insights into gene function in cancer through looking at the known functions of correlated genes in the TCGA.</p> 
         "),
                                     HTML("<p> This analysis is computationally intensive and may take a while, <b> please wait </b></p> 
         "),
                                    ####a little reminder of the inputs selected
                                     fluidRow(
                                       column(2,
                                              DT::dataTableOutput("Inputs")
                                       )
                                     ),
                                     ###table of CCGSEA results
                                      fluidRow(
                                       column(8,
                                              DT::dataTableOutput("Hallmark_GSEA_Table")
                                       )
                                     ),
                                     
                                     
                                     
                                     
                                     #download button - gsea table
                                     downloadButton('Download2','Download GSEA Table'),
                                     
                                     
                                     ###volcano
                                     plotOutput("volcano_plot"), 
                                     
                                     #download button volcano
                                     downloadButton('volc_download','Download Volcano Plot'),
                                    
                                      # Select variable for bottom graphs
                                           selectInput(inputId = "y", 
                                                      label = "Gene Set",
                                                     choices= NULL),
                                     
                                     ##plot the graphs
                                     ##enrichment plot and then a plot of rank ordered correlation coefs
                                     plotOutput("H_plot1"),
                                     plotOutput("H_plot2"),
                                     ##summary stats table
                                     fluidRow(
                                       column(4,
                                              dataTableOutput("Hallmark_table1"))
                                     ),
                                     
                                     ##download buttons for enrichment plot, rank plot, stats table and leading edge and tumour IDS
                                     downloadButton('downloadPlot1', 'Download Enrichment Plot'),
                                     downloadButton('Download_rankplot','Download Correlation Plot'),
                                     downloadButton('Download4','Download Statistics Table'),
                                     downloadButton('Download_LE','Download Leading Edge'),
                                     downloadButton('Download_tumour_IDs2','Download Tumour IDs')
                            ),
                            
                            
                            
                            
                          
                        
                      
                        
                      
         



########################################################################################
tabPanel("CC-GSEA (Custom)",
         
        
             ##CCGSEA custom table
                        br(),
                        HTML("<p> This runs CC-GSEA analysis using the gene, tumour type, correlation method and excluded patients selected in the sidebar with the entered custom gene set. </p> 
         "),
         HTML("<p> This analysis is computationally intensive and may take a while, <b> please wait </b></p> 
         "),
                        fluidRow(
                          column(8,
                                 DT::dataTableOutput("Hallmark_GSEA_Table2")
                          )
                        ),
                        
                        
                        
                        
                        #download button- gsea table
                        downloadButton('Download_GSEA_custom','Download GSEA Table'),
                        
                        
                        
                        
                        ##plot the graphs- enrichment and rank ordered list
                        
                        plotOutput("H_plot12"),
                        plotOutput("H_plot22"),
                        ##summary stats table
                        fluidRow(
                          column(4,
                                 dataTableOutput("Hallmark_table12"))
                        ),
                        
                        ##download buttons for enrichment plot, rank plot, stats table and leading edge and tumour IDs
                        downloadButton('download_enrich', 'Download Enrichment Plot'),
                        downloadButton('Download_rankplot2','Download Correlation Plot'),
                        downloadButton('Download_stats','Download Statistics Table'),
                        downloadButton('Download_LE2','Download Leading Edge'),
                        downloadButton('Download_tumour_IDs3','Download Tumour IDs')
         
                ),
id = "conditionedPanels"  
                      ),
)
               )
),
               
               
               
               
           



##################################################################################
######################################################################################
######################################################################

##the rest of the pages
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
           "The survival analysis tool allows the user to see the relationship of a user inputted gene set with survival of patients in the TCGA in a Kaplan-Meier plot. The mean expression of genes in the selected gene set and tumour type is calculated, 
               this is then used to generate a Kaplan-Meier plot either splitting data at the median or at the optimal point for discrimination between high and low gene set expression. A histogram of distributions of gene expression and the cut-off is also shown.",
           br(),
           id = "mybox5",
           imageOutput("gs_surv_img3")
         )
         
         
),


             
  )
)
  
    
#######################################################################
#######################################################################

server <- function(input, output, session) {


  ##get just the selected set of pathways to show up in the box on page 2
  #reactive UI
  path_options <- reactive({
    ## pick out the set of gene sets chosen
    pathwaysel<-pathway_choices%>%
      dplyr::filter(pathway_choices==input$c)
    pathwaysel<-as.character(pathwaysel)
    ##dplyr::filter the list of pathway options to pick out only the set of gene sets chosen
    pathwayz<-paths%>%
      dplyr::filter(type==pathwaysel)
  })
  ##update the options available for input$y
  observeEvent(path_options(), {
    choices <- unique(path_options()$pathway)
    updateSelectInput(inputId = "y", choices = choices) 
  })
  
  
 ####REV 
  ##get just the selected set of tumour type tumours to show up in the box 
  #reactive UI
  tum_options <- reactive({
    ## pick out the chosen tumour type
    tumsel<-tum_choices%>%
      dplyr::filter(tum_choices==input$z)
    tumsel<-as.character(tumsel)
    ##dplyr::filter to pick only barcodes of tumour type chosen
    tumz<-tumnames%>%
      dplyr::filter(cdr_type==tumsel)
  })
  ##update the options available for input$tumour
  observeEvent(tum_options(), {
    choices_tum <- unique(tum_options()$cdr_bcr_patient_barcode)
    updateSelectInput(inputId = "tumour", choices = choices_tum) 
  })
  
  
  ###########make the very initial correlation table - j()
##make it happen on the go button
  
 j<- eventReactive(
    eventExpr = input$run, {

      
      ##select the right things- tumour type
        data=tcga
      tumour_type<-tcga%>%
        dplyr::filter(cdr_type==input$z)
      
      ##exclude tumours 
      tumours_excluded<-as.character(input$tumour)
      tumour_type<-tumour_type%>%
        dplyr::filter(!cdr_bcr_patient_barcode %in% tumours_excluded)
      
      ##select gene
      goi<-tumour_type%>%
        select(input$x)
      gene<-colnames(goi)
      ### do the correlation
      cor_coef<-cor(tumour_type[3:20272], subset(tumour_type,select=gene), method=tolower(input$corr), use="pairwise.complete.obs") 
      
      ##make it a nice table
      cor_coef<-as.data.frame(cor_coef)
      cor_coef$gene<-rownames(cor_coef)
      cor_coef <- cor_coef[, c(2,1)]
      cor_coef<-na.omit(cor_coef)
      colnames(cor_coef)<-c("Gene","Correlation Coefficient")
      cor_coef<-cor_coef %>%
        as_tibble() %>%
        arrange(desc(`Correlation Coefficient`))
       return(cor_coef)
    }
  )

##output the table
  
  output$Cor_Table <-renderDataTable({
    j()
  })
  
  ##make a little table of the selected inputs
  inps<-reactive({
 
    
    tumour_type<-as.character(input$z)
    tumour_type<-as.data.frame(tumour_type)
    gene<-as.character(input$x)
    gene<-as.data.frame(gene)
    
    
    inputs<-cbind(gene, tumour_type)
    colnames(inputs)<-c("Gene","Tumour Type")
    return(inputs)
          })
  
  output$Inputs <-DT::renderDataTable({
    inps()
  })
  
  #######use j to make a hallmark GSEA table of results 
  b<-reactive({
    
    
    ##make j into a ranked list
    cor_coef<-j()
    colnames(cor_coef)<-c("gene","Correlation_Coefficient")
    cor_coef<-as.data.frame(cor_coef)
    res_for_GSEA<-cor_coef$gene
    res_for_GSEA<-as.data.frame(res_for_GSEA)
    
    
    res_for_GSEA<-cbind(res_for_GSEA,cor_coef$Correlation_Coefficient)
    colnames(res_for_GSEA)<-c("gene_name","cor")
    res_for_GSEA<-na.omit(res_for_GSEA)
    
    ranks <- res_for_GSEA$cor
    names(ranks) <- res_for_GSEA$gene_name
    
    a<-sort(ranks, decreasing = T)
    
    ##select the pathways to use
    pathwaysfgsea<-all_paths[[input$c]]
    
    ##gsea it and make a nice table
    hallmark_fgseaRes<-fgseaMultilevel(pathways=pathwaysfgsea, stats=a, BPPARAM=SnowParam(4))
    hallmark_fgseaRes<-as.data.frame(hallmark_fgseaRes)
    colnames(hallmark_fgseaRes)<-c("Pathway","pval","padj","log2err","ES","NES","Size","Leading Edge")
    
    hallmark_fgseaRes<-hallmark_fgseaRes %>%
      as_tibble() %>%
      arrange(desc(NES))
    
    
    return(hallmark_fgseaRes)
    
  })
  ##output the table
  output$Hallmark_GSEA_Table <-DT::renderDataTable({
    b()
    
  },
  options = list(scrollX = TRUE))##make it scrollable across
  
  
  
  
  
  
  
  #################volcano plot fom b
  volcano<-reactive({
    volc<-b()
    
    ## pick out the diff expressed ones to be a different colour
    volc$diffexpressed<-"NO"
   volc$diffexpressed[volc$padj<0.05] <- "UP"
    
   
   ##pick out ones to label
    volc<-volc[order(volc$padj),]
    ppaths<-volc[1:3,]
    ppaths<-ppaths$Pathway
    volc<-volc[order(volc$NES),]
    npathsdn<-volc[1:3,]
    npathsdn<-npathsdn$Pathway
    volc<-volc[order(-volc$NES),]
    npathsup<-volc[1:3,]
    npathsup<-npathsup$Pathway
    paths<-c(ppaths,npathsdn,npathsup)
    paths<-unique(paths)
    
    
    volc$delabel <- NA
    volc$delabel[volc$Pathway %in% paths] <- volc$Pathway[volc$Pathway %in% paths]
 
    
    ##make volcano
    volcplot<-ggplot(volc,aes(x=NES,y=-log10(padj),col=diffexpressed,label=delabel))+geom_point()+
       theme_classic()+
      geom_hline(yintercept=-log10(0.05), col="red",linetype="dashed")+
      scale_color_manual(values=c("black", "red"))+
      geom_text_repel() +
      theme(legend.position = "none")
    
    print(volcplot)
    
  })
  
##output the volcano  
  output$volcano_plot<-renderPlot({
    volcano()
  })
  
  
  
  
  
  
  
  
  
  ############ use j to plot an enrichment pathway
  d<-reactive({
   
    ##get the selected pathway
     pathw<-as.data.frame(pathws)
    pathway<-pathw%>%
      dplyr::filter(pathws==input$y)
    pathway<-as.character(pathway)
    
    ##get the correlation coefficients -j
    cor_coef<-j()
    colnames(cor_coef)<-c("gene","Correlation_Coefficient")
    
    cor_coef<-as.data.frame(cor_coef)
    res_for_GSEA<-cor_coef$gene
    res_for_GSEA<-as.data.frame(res_for_GSEA)
    
    
    res_for_GSEA<-cbind(res_for_GSEA,cor_coef$Correlation_Coefficient)
    colnames(res_for_GSEA)<-c("gene_name","cor")
    res_for_GSEA<-na.omit(res_for_GSEA)
    
    
    ##make ranked list
    ranks <- res_for_GSEA$cor
    names(ranks) <- res_for_GSEA$gene_name
    
    a<-sort(ranks, decreasing = T)
    
    #isolate the corect list of pathways to use in enrichment
    
    pathwaysfgsea<-all_paths[[input$c]]
##plot enrichment
    p<-plotEnrichment(pathway = pathwaysfgsea[[pathway]], a)+ labs(title=input$y)
    print(p)
  })
  
  ##plot 
  output$H_plot1<-renderPlot({
    d()
  })
  #############use j to plot the bar graph of correlation
  e<-reactive({
    
    cor_coef<-j()
    colnames(cor_coef)<-c("gene","Correlation_Coefficient")
    ranks <- cor_coef$Correlation_Coefficient
    names(ranks) <- cor_coef$gene
    ##rank
    a<-sort(ranks, decreasing = T)
    
    
    
    
    return(a)
  })

  ##plot the bar graph
  output$H_plot2<-renderPlot({
     data=e()
    names(data)<-NULL
    data<-as.data.frame(data)
    data$rank<-rev(row_number(data))
    ggplot(data,aes(x=rank,y=data))+geom_bar(stat="identity",fill="black",color="black")+
      theme_classic()+
      xlab("Rank")+
      ylab("Correlation Coefficient")+
      ggtitle("Genes Ranked by Correlation Coefficient")

    
  })

  
  
  
  ########make the table of extra info bits
  t<-reactive({
    #pathway
    pathw<-as.data.frame(pathws)
    pathway<-pathw%>%
      dplyr::filter(pathws==input$y)
    pathway<-as.character(pathway)
    #tumour type
    tumour_type<-tcga%>%
      dplyr::filter(cdr_type==input$z)
    ##collection of gene sets
    pathwaysfgsea<-all_paths[[input$c]]
    ##GSEA table
    hallmark_fgseaRes<-b()
    hallmark_fgseaRes<-as.data.frame(hallmark_fgseaRes)
    ##get the relevant info out of the gsea table
    POI<-hallmark_fgseaRes%>%
      dplyr::filter(Pathway==input$y)
    
    ##number of tumours analysed
    tumours_excluded<-as.character(input$tumour)
    tumour_type<-tumour_type%>%
      dplyr::filter(!cdr_bcr_patient_barcode %in% tumours_excluded)
    r<-c(nrow(tumour_type))
    #pathway
    s<-c(POI$Pathway)
    #padj
    u<-c(POI$padj)
    #NES
    v<-c(POI$NES)
    #put it together to make a table
    w<-rbind(r,s,u,v)
    w<-as.data.frame(w)
    w<-transpose(w)
    colnames(w)<-c("Number of Tumours","Pathway","padj","NES")
    return(w)
    
    
  })
  ##output the table
  output$Hallmark_table1<-renderDataTable({
    t()
  })
  
  ###download the leading edge- get it first
  lead_edge<-reactive({
   #GSEA table
     hallmark_fgseaRes<-b()
    hallmark_fgseaRes<-as.data.frame(hallmark_fgseaRes)
    ##select the pathway
    POI<-hallmark_fgseaRes%>%
      dplyr::filter(Pathway==input$y)
    ##get the leading edge
    POI<-hallmark_fgseaRes%>%
      select(`Leading Edge`)
    ##get it suitable for a csv
    POI<-as.data.frame(unlist(POI))
    colnames(POI)<-paste0(input$y,"_",input$x,"_leading_edge")
    
    
    return(POI)
    
  })
  #########################################downloads
  
##download leading edge  
output$Download_LE<- downloadHandler(
  filename = function(){
    paste0(input$c,"_",input$x,"_",input$z,"_",input$y,"_leading_edge.csv",sep="")
  },
  content = function(file) {
    write.csv(lead_edge(),file,row.names=FALSE)
  }
)
  
  ##donload corr table
  output$Download1 <- downloadHandler(
    filename = function(){
      paste0(input$x,"_correlation_table.csv")
    },
    content = function(file) {
      write.csv(j(),file,row.names=FALSE)
    }
  )
  #doenload gsea table
  output$Download2 <- downloadHandler(
    filename = function(){
      paste0(input$c,"_",input$x,"_",input$z,"_GSEA_table.csv",sep="")
    },
    content = function(file) {
      write.csv(b()[,1:7], file,row.names=FALSE)
    }
  )
  
  
  ##download summary stats
  output$Download4 <- downloadHandler(
    filename = function(){
      paste0(input$c,"_",input$x,"_",input$z,"_",input$y,"_statistics.csv",sep="")
    },
    content = function(file) {
      write.csv(t(), file,row.names=FALSE)
    }
  )
  
  
  ##download volcano
  
  output$volc_download<- downloadHandler(
    filename = function(){
      paste0(input$c,"_",input$x,"_",input$z,"_volcano_plot.svg")
    },
    content = function(file) {
      svg(file)
      print(volcano())
      dev.off()
    }
  )
  
  
  ##download enrichment plot
  
  output$downloadPlot1 <- downloadHandler(
    filename = function(){
      paste0(input$c,"_",input$x,"_",input$z,"_",input$y,"_enrichment_plot.svg")
    },
    content = function(file) {
      svg(file)
      print(d())
      dev.off()
    }
  )
  
  
  ##download the rank plot
  output$Download_rankplot<-downloadHandler(
    filename = function(){
      paste0(input$x,"_",input$z,"_correlation_plot.svg")
    },
    content = function(file) {
      svg(file)
      data=e()
      names(data)<-NULL
      data<-as.data.frame(data)
      data$rank<-rev(row_number(data))
      plot<-ggplot(data,aes(x=rank,y=data))+geom_bar(stat="identity",fill="black",color="black")+
        theme_classic()+
        xlab("Rank")+
        ylab("Correlation Coefficient")+
        ggtitle("Genes Ranked by Correlation Coefficient")
      print(plot)
      dev.off()
    }
  )
    ### download the tumour IDs
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
  
   output$Download_tumour_IDs2<-downloadHandler(
     filename = "Tumours_IDs_in_analysis.csv",
     content = function(file) {
       write.csv(tummys(), file,row.names=FALSE)
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
  
  
  
  #############GS-Corr page GS-Corr diagram
  output$gs_corr_img2 <- renderImage({
    
    list(src = "www/GS-Corr_diagram.png",
         height=700,
         width=750)
    
  }, deleteFile = F)
  
  
  
  
#############################################################################  
  #############################################################################
  #############################################################################
##CC-GSEA custom
  
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
  
  ##get just the selected set of pathways to show up in the box on page 2
  #reactive UI
 
  
  
  
    ###
  ##get just the selected set of tumour type tumours to show up in the box 
  tum_options2 <- reactive({
    ## pick out the tumour type chosen
    tumsel2<-tum_choices%>%
      dplyr::filter(tum_choices==input$tum_type)
    tumsel2<-as.character(tumsel2)
    ##dplyr::filter pick only tumours in selected tumour type
    tumz2<-tumnames%>%
      dplyr::filter(cdr_type==tumsel2)
  })
  ##update the options available for input$tumour_excl
  observeEvent(tum_options2(), {
    choices_tum2 <- unique(tum_options2()$cdr_bcr_patient_barcode)
    updateSelectInput(inputId = "tumour_excl", choices = choices_tum2) 
  })
  
  
  ###########make the very initial correlation table - j2()
  ##make it happen on the go button
  
  j2<- eventReactive(
    eventExpr = input$run_custom, {
      
      
      ##select the right things
      data=tcga
      tumour_type<-tcga%>%
        dplyr::filter(cdr_type==input$tum_type)
      
      ##remove excluded tumours
      tumours_excluded<-as.character(input$tumour_excl)
      tumour_type<-tumour_type%>%
        dplyr::filter(!cdr_bcr_patient_barcode %in% tumours_excluded)
      
      
      goi<-tumour_type%>%
        select(input$x2)
      gene<-colnames(goi)
      ### do the correlation
      cor_coef<-cor(tumour_type[3:20272], subset(tumour_type,select=gene), method=tolower(input$corr2), use="pairwise.complete.obs") 
      
      ##make it a nice table
      cor_coef<-as.data.frame(cor_coef)
      cor_coef$gene<-rownames(cor_coef)
      cor_coef <- cor_coef[, c(2,1)]
      cor_coef<-na.omit(cor_coef)
      colnames(cor_coef)<-c("Gene","Correlation Coefficient")
      cor_coef<-cor_coef %>%
        as_tibble() %>%
        arrange(desc(`Correlation Coefficient`))
      return(cor_coef)
    }
  )
  
  ##output the table
  
  output$Cor_Table2 <-renderDataTable({
    j2()
  })
  
  
 
  

  #######use j to make a hallmark GSEA table of results 
  b2<-reactive({
   
    
    ##make j into a ranked list
    cor_coef<-j2()
    colnames(cor_coef)<-c("gene","Correlation_Coefficient")
    cor_coef<-as.data.frame(cor_coef)
    res_for_GSEA<-cor_coef$gene
    res_for_GSEA<-as.data.frame(res_for_GSEA)
    
    
    res_for_GSEA<-cbind(res_for_GSEA,cor_coef$Correlation_Coefficient)
    colnames(res_for_GSEA)<-c("gene_name","cor")
    res_for_GSEA<-na.omit(res_for_GSEA)
    
    ranks <- res_for_GSEA$cor
    names(ranks) <- res_for_GSEA$gene_name
    
    a<-sort(ranks, decreasing = T)
    
    ##select the pathways to use
    data<-data()
    colnames(data)<-"Genes"
    data$Genes<-toupper(data$Genes)
    genes_to_use<-data$Genes
    genes_to_use<-as.character(genes_to_use)
    user_inputted_list<-list(genes_to_use)
    names(user_inputted_list)<-"User_Inputted_List"

    
    ##gsea it and make a nice table
    hallmark_fgseaRes2<-fgsea(pathways=user_inputted_list, stats=a)
    hallmark_fgseaRes2<-as.data.frame(hallmark_fgseaRes2)
    colnames(hallmark_fgseaRes2)<-c("Pathway","pval","padj","log2err","ES","NES","Size","Leading Edge")
    
    
    return(hallmark_fgseaRes2)
    
  })
  ##output the table
  output$Hallmark_GSEA_Table2 <-DT::renderDataTable({
    b2()
    
  },
  options = list(scrollX = TRUE))##make it scrollable across
  
  
  
  
  ############ use j to plot an enrichment pathway
  d2<-reactive({
    
    ##get the correlation coefficients -j
    cor_coef<-j2()
    colnames(cor_coef)<-c("gene","Correlation_Coefficient")
    
    cor_coef<-as.data.frame(cor_coef)
    res_for_GSEA<-cor_coef$gene
    res_for_GSEA<-as.data.frame(res_for_GSEA)
    
    
    res_for_GSEA<-cbind(res_for_GSEA,cor_coef$Correlation_Coefficient)
    colnames(res_for_GSEA)<-c("gene_name","cor")
    res_for_GSEA<-na.omit(res_for_GSEA)
    
    
    ##make ranked list
    ranks <- res_for_GSEA$cor
    names(ranks) <- res_for_GSEA$gene_name
    
    a<-sort(ranks, decreasing = T)
    
    #isolate the corect list of pathways to use in enrichment
    data<-data()
    colnames(data)<-"Genes"
    data$Genes<-toupper(data$Genes)
    genes_to_use<-data$Genes
    genes_to_use<-as.character(genes_to_use)
    user_inputted_list<-list(genes_to_use)
    names(user_inputted_list)<-"User_Inputted_List"
    ##plot enrichment
    p<-plotEnrichment(pathway = user_inputted_list[[1]], a)+ labs(title="Enrichment Plot")
    print(p)
  })
  
  ##plot 
  output$H_plot12<-renderPlot({
    d2()
  })
  #############use j to plot the bar graph of correlation
  e2<-reactive({
    
    cor_coef<-j2()
    colnames(cor_coef)<-c("gene","Correlation_Coefficient")
    ranks <- cor_coef$Correlation_Coefficient
    names(ranks) <- cor_coef$gene
    ##rank
    a<-sort(ranks, decreasing = T)
    
    
    
    
    return(a)
  })
  
  ##plot the bar graph
  output$H_plot22<-renderPlot({
    data=e2()
    names(data)<-NULL
    data<-as.data.frame(data)
    data$rank<-rev(row_number(data))
    ggplot(data,aes(x=rank,y=data))+geom_bar(stat="identity",fill="black",color="black")+
      theme_classic()+
      xlab("Rank")+
      ylab("Correlation Coefficient")+
      ggtitle("Genes Ranked by Correlation Coefficient")
    
    
  })
  
  
  
  
  ########make the table of extra info bits
  t2<-reactive({
    #tumour type
    tumour_type<-tcga%>%
      dplyr::filter(cdr_type==input$tum_type)
    ##GSEA table
    hallmark_fgseaRes<-b2()
    hallmark_fgseaRes<-as.data.frame(hallmark_fgseaRes)
    ##get the relevant info out of the gsea table
    POI<-hallmark_fgseaRes%>%
      dplyr::filter(Pathway=="User_Inputted_List")
    
    ##number of tumours analysed
    tumours_excluded<-as.character(input$tumour_excl)
    tumour_type<-tumour_type%>%
      dplyr::filter(!cdr_bcr_patient_barcode %in% tumours_excluded)
    r<-c(nrow(tumour_type))
    #padj
    u<-c(POI$padj)
    #NES
    v<-c(POI$NES)
    #put it together to make a table
    w<-rbind(r,u,v)
    w<-as.data.frame(w)
    w<-transpose(w)
    colnames(w)<-c("Number of Tumours","padj","NES")
    return(w)
    
    
  })
  ##output the table
  output$Hallmark_table12<-renderDataTable({
    t2()
  })
  
  ###download the leading edge- get it first
  lead_edge2<-reactive({
    #GSEA table
    hallmark_fgseaRes<-b2()
    hallmark_fgseaRes<-as.data.frame(hallmark_fgseaRes)
    ##select the pathway
    POI<-hallmark_fgseaRes%>%
      dplyr::filter(Pathway=="User_Inputted_List")
    ##get the leading edge
    POI<-hallmark_fgseaRes%>%
      select(`Leading Edge`)
    ##get it suitable for a csv
    POI<-as.data.frame(unlist(POI))
    colnames(POI)<-"leading_edge"
    
    
    return(POI)
    
  })
  #########################################downloads
  
  ##download leading edge  
  output$Download_LE2<- downloadHandler(
    filename = function(){
      paste0(input$tum_type,"_custom_leading_edge.csv",sep="")
    },
    content = function(file) {
      write.csv(lead_edge2(),file,row.names=FALSE)
    }
  )
  
  ##donload corr table
  output$Download_corr_table <- downloadHandler(
    filename = function(){
      paste0(input$x2,"_",input$tum_type,"_correlation_table.csv")
    },
    content = function(file) {
      write.csv(j2(),file,row.names=FALSE)
    }
  )
  #doenload gsea table
  output$Download_GSEA_custom <- downloadHandler(
    filename = function(){
      paste0(input$x2,"_",input$tum_type,"_GSEA_table.csv",sep="")
    },
    content = function(file) {
      write.csv(b2()[,1:7], file,row.names=FALSE)
    }
  )
  
  
  ##download summary stats
  output$Download_stats <- downloadHandler(
    filename = function(){
      paste0(input$x2,"_",input$tum_type,"_statistics.csv",sep="")
    },
    content = function(file) {
      write.csv(t2(), file,row.names=FALSE)
    }
  )
  
  

  
  
  ##download enrichment plot
  
  output$download_enrich <- downloadHandler(
    filename = function(){
      paste0(input$x2,"_",input$tum_type,"_enrichment_plot.svg")
    },
    content = function(file) {
      svg(file)
      print(d2())
      dev.off()
    }
  )
  
  
  ##download the rank plot
  output$Download_rankplot2<-downloadHandler(
    filename = function(){
      paste0(input$x2,"_",input$tum_type,"_correlation_plot.svg")
    },
    content = function(file) {
      svg(file)
      data=e2()
      names(data)<-NULL
      data<-as.data.frame(data)
      data$rank<-rev(row_number(data))
      plot<-ggplot(data,aes(x=rank,y=data))+geom_bar(stat="identity",fill="black",color="black")+
        theme_classic()+
        xlab("Rank")+
        ylab("Correlation Coefficient")+
        ggtitle("Genes Ranked by Correlation Coefficient")
      print(plot)
      dev.off()
    }
  )
  ### download the tumour IDS

  tummys2<-reactive({
    tumours_excluded<-as.character(input$tumour_excl)
    tumour_type<-tcga%>%
      dplyr::filter(cdr_type==input$tum_type)%>%
      dplyr::filter(!cdr_bcr_patient_barcode %in% tumours_excluded)%>%
      select(cdr_bcr_patient_barcode)
    return(tumour_type)
  })
  
  output$Download_tumour_IDs3<-downloadHandler(
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
  
  
  
  #############GS-Corr page GS-Corr diagram
  output$gs_corr_img2 <- renderImage({
    
    list(src = "www/GS-Corr_diagram.png",
         height=700,
         width=750)
    
  }, deleteFile = F)
  
  
  #########CC-GSEA page CC-GSEA diagram
  
  output$cc_gsea_img2 <- renderImage({
    
    list(src = "www/ccgsea_diagram.png",
         height=700,
         width=750)
    
  }, deleteFile = F)
  
  
  
  
  
  #############################################################################
  #############################################################################
  #############################################################################
  
  
  
  
  
}

# Create a Shiny app object ----------------------------------------------------

shinyApp(ui = ui, server = server)

