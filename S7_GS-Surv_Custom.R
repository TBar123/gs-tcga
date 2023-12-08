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

options(rsconnect.max.bundle.size=3145728000)
# Load data --------------------------------------------------------------------

##gene expression tumour type patient barcode
tcga<-read_fst("tcga_w_barcode2.fst")

##clinical data
cdr_clinicaldata<-read.delim("clinical_data.tsv",sep="\t")

##make list of patient barcodes
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
     
     
     
     
     
  #########################GS-Surv Custom UI-------------------------   
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
                          
                          ###REV - select tumours to exclude
                          selectizeInput(inputId = "tumour", 
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
                       
                          ##go button
                             actionButton("run", label = "Run GS-Surv"),
                          width = "3",
                     
                          
                          
                         
                          
                          HTML("<p>For more information about tumour type abbreviations please click <a href='https://gdc.cancer.gov/resources-tcga-users/tcga-code-tables/tcga-study-abbreviations'>here</a></p>
                               <p>For annotations of individual cases please click <a href='https://portal.gdc.cancer.gov/annotations'>here</a></p>")
                        ),
                          
            

                          
                        
                      
                      ########################################################
                      #######################################################main bit
                        mainPanel(
                          
                          ## a page showing the input genes
                          tabsetPanel(
                            tabPanel("Uploaded Genes",
                          tableOutput("contents")
                        ),
                        ###median split panel
                        tabPanel("Median split",
                                 
                                 ##kaplan
                                 plotOutput("fifty_plot1"),
                                
                                 ##download survival plot
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
                                 ##download histogram and tumour info
                                 downloadButton('downloadPlotfifhist','Download Average Expression Plot'),
                                 downloadButton('Download_tumour_IDs','Download Tumour Case Data'),
                                 
                                 HTML("<p>The histogram shows the median gene expression of each tumour after scaling. The red line shows the median cut off point used to separate high and low tumours for survival analysis.</p>"),
                    
                                 
                        ),
                        
                        ###thirds split panel
                        tabPanel("Split by thirds",
                                 
                                 ##kaplan
                                 plotOutput("third_plot1"),
                                 
                                 ##download survival plot
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
                                 downloadButton('download_coxph_ps','Download CoxPH Table'),
                                 
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
                                 ##download histogram and tumour info
                                 downloadButton('downloadPlotthirdhist','Download Average Expression Plot'),
                                 downloadButton('Download_tumour_IDs2','Download Tumour Case Data'),
                                 
                                 HTML("<p>The histogram shows the median gene expression of each tumour after scaling. The red line shows the cut off points used to separate high, medium and low tumours for survival analysis.</p>"),
                                 
                                 
                        ),
                        
                        
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

 
  ##get just the selected set of tumour type tumours to show up in the box 
  #reactive UI
  tum_options <- reactive({
    ## pick out the tumour type chosen
    tumsel<-tum_choices%>%
      dplyr::filter(tum_choices==input$upld_tt)
    tumsel<-as.character(tumsel)
    ##dplyr::filter the list of barcode options to pick out only ones of selected tumour type
    tumz<-tumnames%>%
      dplyr::filter(cdr_type==tumsel)
  })
  ##update the options available for input$tumour
  observeEvent(tum_options(), {
    choices_tum <- unique(tum_options()$cdr_bcr_patient_barcode)
    updateSelectInput(inputId = "tumour", choices = choices_tum) 
  })
  
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
    
    tumours_excluded<-as.character(input$tumour)
    ##pick out the tumour type and exclude any tumours
    tumour_type<-tcga%>%
    filter(cdr_type==input$upld_tt)%>%
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
  
 

  ###########################third cut##################################################################
  ##median cut off
  third<-reactive({
    
    
    means_and_clin=preprocess()
   
    third_cut<-means_and_clin%>%
      arrange(mns)
    ###cut it and label each third
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
    third_cut$Category <- factor(third_cut$Category, levels = c('Low','Middle','High'))
    
    return(third_cut)
    
    
  })
  
  ##kaplan plot it- thirds
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
  
  ##histogram- thirds
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
  
  ##output hist - thirds
  output$third_hist<-renderPlot({hist_30()})
  
  ##coxph pvals
  thirds_pvals<-reactive({ 
    coxfit <- coxph(
      Surv(OS.time,OS) ~ Category,
      data = third(),
      ties = 'exact')
    
    # summary(coxfit)
    
    fit<-summary(coxfit)
    a<-fit$coefficients
    a<-as.data.frame(a)
    a<-a[,c(2,5)]
    colnames(a)<-c("CoxPH HR","CoxPH p value")
    rownames(a)<-c("Middle vs Low", "High vs Low")
    a$Comparison<-rownames(a)
    a<-a[,c(3,1,2)]
    return(a)
    
  })  
  
  ##output it
  output$coxph_ps<-renderDataTable({
    thirds_pvals()
  })
  
  
  
  ######################################################################################################
  
  
  
  
  
  ##median cut off
  fifty<-reactive({
    
    
    means_and_clin=preprocess()
    ##put high and low above and below median
    fifty_cut<-means_and_clin%>%
      mutate(Category=ifelse(mns>=median(mns),"high","low"))%>%
      select(OS.time,OS,Category)
    
    colnames(fifty_cut)<-c("OS.time","OS","Category")
    fifty_cut$Category <- factor(fifty_cut$Category, levels = c('low','high'))
    
    return(fifty_cut)
    
    
  })
  
  
  
  
  
  ##kaplan plot it- median
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
  ##histogram- median
  hist_50<-reactive({
    
  
  means_and_clin=preprocess()
  
##hist - median
  toplot<-ggplot(means_and_clin,aes(x=mns))+
    geom_histogram(bins=50)+
    xlab("Average Pathway Expression")+
    theme_bw()+
    geom_vline(xintercept=median(means_and_clin$mns),linetype="dashed",color="red")
  print(toplot)
  
  
  })
  
  ##output hist- median
  output$ave_hist<-renderPlot({hist_50()})
  
  #####################################COXPH############################################################################
  ##get coxph Hazard ratio
  ##### median as continuous var
  cox<-reactive({
    means_and_clin=preprocess()
    
    fit<-coxph(Surv(OS.time, OS)~mns, data=means_and_clin)
    fit<-summary(fit)
    f<-fit$coefficients
    f<-as.data.frame(f)
    f<-f[,c(2,5)]
    colnames(f)<-c("CoxPH HR- median gene set expression as continuous variable","CoxPH p value - median gene set expression as continuous variable")
    f<-f[,c(1,2)]
    
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
    colnames(a)<-c("CoxPH HR - categorical","CoxPH p value - categorical")
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
    colnames(p)<-"Logrank p value"
    
    pvals<-cbind(f,a,p)
    rownames(pvals)<-NULL
    
    return(pvals)
    
    
  })
  
  
  ##output it
  output$median_ps<-renderDataTable({
    cox()
  })
  
  
  ###downlaods
 
  
  output$downloadPlotsurvfif <- downloadHandler(
    filename = "median_cutoff_survival_plot.svg"
    ,
    content = function(file) {
      svg(file)
      print(ggsurvplot((survfit(Surv(OS.time,OS)~Category,data=fifty())  ),data=fifty(),pval=TRUE))
      dev.off()
    }
  )
  
  
  output$downloadPlotsurvthird <- downloadHandler(
    filename = "third_cutoff_survival_plot.svg"
    ,
    content = function(file) {
      svg(file)
      print(ggsurvplot((survfit(Surv(OS.time,OS)~Category,data=third())  ),data=third()))
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
  
  output$downloadPlotthirdhist <- downloadHandler(
    filename = "third_cutoff_histogram.svg"
    ,
    content = function(file) {
      svg(file)
      print(hist_30())
      dev.off()
    }
  )
  
  
  output$download_coxph_ps<-downloadHandler(
    filename = function(){
      paste0(input$upld_tt,"_CoxPH.csv")
    },
    content = function(file) {
      write.csv(thirds_pvals(), file,row.names=FALSE)
    }
  )
  
  output$download_med_ps<-downloadHandler(
    filename = function(){
      paste0(input$upld_tt,"_statistics.csv")
    },
    content = function(file) {
      write.csv(cox(), file,row.names=FALSE)
    }
  )
  
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
  
  ###############GS-Surv page, GS-Surv image
  output$gs_surv_img2 <- renderImage({
    
    list(src = "www/gs_surv_diagram.png",
         height=1200,
         width=800)
    
  }, deleteFile = F)
  
  
}

# Create a Shiny app object ----------------------------------------------------

shinyApp(ui = ui, server = server)

