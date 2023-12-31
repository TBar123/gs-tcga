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

# Load data --------------------------------------------------------------------

##gene expression tumour type patient barcode
tcga<-read_fst("tcga_w_barcode.fst")
##extract gene names
colheads<-colnames(tcga)
inputnames<-colheads[3:20503]
##clinical data
cdr_clinicaldata<-fread('TCGA-CDR-SupplementalTableS1_sampledata.tsv', verbose=TRUE, data.table=FALSE, header=TRUE)




#######UI

ui <- fluidPage(
  useShinydashboard(),
  #shinyjs and inlineCSS allow loading page
  useShinyjs(),
 
  
  #GS-Surv GS-Surv page box
  tags$head(tags$style("
                     
                     #mybox4{height:1300px !important;}
                 
                
                     ")),
  
  
  #CC-GSEA CC-GSEA page box
  tags$head(tags$style("
                     
                     #mybox6{height:800px !important;}
                 
                                 ")),
  
  ##GS-Corr GS-Corr page box
  tags$head(tags$style("
                     
                     #mybox7{height:775px !important;}
                 
                
                     ")),
  
  

  
  
  ###make a thing to go under my histogram saying high or low
  tags$footer(tags$style("
                     
                     #footer {
                    display: flex;
                    justify-content: space-between;
                    color:black;
                    height:30px;
                    align-items:center;
                    
               }
                
                     
                     
                     ")),
  
  navbarPage("GS-TCGA",selected="GS-Surv (Custom)",
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
     
     
     
     
     
     
       ##current page      
             tabPanel("GS-Surv (Custom)",
                      
                      ##side bar
                      sidebarLayout(
                        # Inputs: Select variables to plot
                        sidebarPanel(
                          
                        
                          #Select variable for tumour type
                          selectInput(inputId = "upld_tt", 
                                      label = "Tumour Type",
                                      choices = c("UVM","ACC","GBM","OV","LUAD","PRAD","UCEC","LUSC","BLCA","TGCT","ESCA","PAAD","KIRP","LIHC","CESC","SARC","BRCA","THYM","MESO","COAD","STAD","SKCM","CHOL","KIRC","THCA","HNSC","READ","LGG","DLBC","KICH","UCS","PCPG"),
                                      selected = "UVM"),
                          
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
                       
                          ##go button
                             actionButton("run", label = "Run GS-Surv"),
                          width = "3",
                     
                          
                          
                         
                          
                          HTML("<p>For more information about tumour type abbreviations please click <a href='https://gdc.cancer.gov/resources-tcga-users/tcga-code-tables/tcga-study-abbreviations'>here</a></p>")
                        ),
                          
            

                          
                        
                      
                      ########################################################
                      #######################################################main bit
                        mainPanel(
                          
                          ## a page showing the input
                          tabsetPanel(
                            tabPanel("Uploaded Genes",
                          tableOutput("contents")
                        ),
                        ###median split panel
                        tabPanel("Median split",
                                 
                                 ##kaplan
                                 plotOutput("fifty_plot1"),
                                 
                                 downloadButton('downloadPlotsurvfif', 'Download Survival Plot'),
                                 
                                 ##histogram
                                 plotOutput("ave_hist"),
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
                                 
                                 downloadButton('downloadPlotfifhist','Download Average Expression Plot'),
                                 
                                 HTML("<p>The histogram shows the mean gene expression of each tumour after scaling. The red line shows the median cut off point used to separate high and low tumours for survival analysis.</p>"),
                    
                                 
                        ),
                        
                        #optimal split
                        tabPanel("Optimal Split",
                                 
                             
                                 ##kaplan
                                 plotOutput("opt_plot1"),
                                 
                                 downloadButton('downloadPlotopt_surv', 'Download Survival Plot'),
                                 
                                 
                                 ##hist
                                 plotOutput("ave_hist_cut"),
                                 ##thing for bottom of hist
                                 HTML('<div class="row">
                                           <div id="footer">
                                           <div >Gene set expression- low</div>
                                           <div ></div>
                                           <div >Gene set expression- high</div>
                                           </div>
                                           </div>
                                           '),
                                 
                                 br(),
                                 
                                 downloadButton('downloadPlotopt_hist','Download Average Expression Plot'),
                                 
                                 HTML("<p>The histogram shows the mean gene expression of each tumour after scaling. The red line shows the cut off point which provides the optimal p value used to separate high and low tumours for survival analysis.</p>"),
                     
                                 
                                 
                        )
                        
                          )
                        
                          
                        )
                      )
                      
)
)
)




  
  
  


server <- function(input, output, session) {
  
  data <- reactiveVal()
  ##make the action button work and do the right things for text and files
  observeEvent(
    eventExpr = input$run,
    handlerExpr = {
      switch(input$choose,
             File = read.csv(input$file1$datapath),
             Text = data.frame(Genes = c(input$caption))%>%
               separate_longer_delim(Genes,delim="\n")#delim=" "
      ) %>% data()
    }
  )

##the uploaded genes output
  output$contents <- renderTable({
    data()
  })
  
##pre processing
  preprocess<-reactive({
  
    ##pick out the genes and barcode
    data<-data()
    colnames(data)<-"Genes"
    data$Genes<-toupper(data$Genes)
    genes_to_use<-data$Genes
    genes_to_use<-c("cdr_bcr_patient_barcode",genes_to_use)
    
    
    ##pick out the tumour type
    tumour_type<-tcga%>%
    filter(cdr_type==input$upld_tt)%>%
      select(any_of(genes_to_use))
  
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
  
  
  ##calculate mean expression
  mns<-rowMeans(tum_ov10,na.rm=T)

  mns<-as.data.frame(mns)
  
  mns$barcode<-rownames(mns)
  
  ##add the clinical on
  means_and_clin<-mns%>%
    left_join(cdr_clinicaldata,by=c("barcode"="bcr_patient_barcode"))
  means_and_clin
  
  })
  
 
  
  
  
  
 ############################### optimal calculating cutpoint and getting it ready
  opt<-reactive({
    
    
    means_and_clin=preprocess()
    ##calculate cut point
    opt_cut<-surv_cutpoint(means_and_clin,time = "OS.time", event = "OS", variables = "mns", minprop=0.1, progressbar = TRUE)
    ##put into high and low
    opt_surv_cat<-surv_categorize(opt_cut)
    opt_surv_cat<-as.data.frame(opt_surv_cat)
  

    colnames(opt_surv_cat)<-c("OS.time","OS","Category")
    
    return(opt_surv_cat)
    
    
  })
  
  
  ##kaplan meier plot
  output$opt_plot1<-renderPlot({
    fit=survfit(Surv(OS.time,OS)~Category,data=opt())  
    ggsurvplot(fit,data=opt(),pval=TRUE,xlab="Time (Days)")
    
  })
  
  
  
  ##optimal histogram
  
  hist<-reactive({
    
    to_use<-preprocess()
   ##get cutpoint 
    opt_cut<-surv_cutpoint(to_use,time = "OS.time", event = "OS", variables = "mns", minprop=0.1, progressbar = TRUE)


  ##extract the cut point  
    cutpt<-opt_cut[[1]]$estimate
 
    
   ##plot histogram 
    toplot<-ggplot(to_use,aes(x=mns))+
      geom_histogram(bins=50)+
      xlab("Average Pathway Expression")+
      theme_bw()+
      geom_vline(xintercept=cutpt,linetype="dashed",color="red")
    print(toplot)
    
  })
  
  
##output histogram  
  output$ave_hist_cut<-renderPlot({hist()})
  
  ##median cut off
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
    ggsurvplot(fit,data=fifty(),pval=TRUE,xlab="Time (Days)")
    
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
  
  ###downlaods
  output$downloadPlotopt_surv <- downloadHandler(
    filename = "opt_cutoff_survival_plot.svg"
    ,
    content = function(file) {
      svg(file)
      print(ggsurvplot((survfit(Surv(OS.time,OS)~Category,data=opt())  ),data=opt(),pval=TRUE))
      dev.off()
    }
  )
  
  
  
  
  output$downloadPlotsurvfif <- downloadHandler(
    filename = "median_cutoff_survival_plot.svg"
    ,
    content = function(file) {
      svg(file)
      print(ggsurvplot((survfit(Surv(OS.time,OS)~Category,data=fifty())  ),data=fifty(),pval=TRUE))
      dev.off()
    }
  )
  
  
  output$downloadPlotfifhist <- downloadHandler(
    filename = "median_cutoff_histogram.svg"
    ,
    content = function(file) {
      svg(file)
      print(hist_50())
      dev.off()
    }
  )
  
  
  
  output$downloadPlotopt_hist <- downloadHandler(
    filename = "opt_cutoff_histogram.svg"
    ,
    content = function(file) {
      svg(file)
      print(hist())
      dev.off()
    }
  )
  
  
  ###############GS-Surv page, GS-Surv image
  output$gs_surv_img2 <- renderImage({
    
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

