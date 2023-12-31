# Load packages ----------------------------------------------------------------

library(BiocManager)
options(repos = BiocManager::repositories())

library(shiny)
library(ggplot2)
library(data.table)
library(dplyr)
library(abind)
library(shinyscreenshot)
library(shinyjs)
library(survival)
library(survminer)
library(bslib)
library(shinydashboard)
library(shinyWidgets)
library(fst)



options(rsconnect.max.bundle.size=3145728000)
# Load data --------------------------------------------------------------------




##this is the patients categorised into high or low for each gene set using the optimal cutpoint from the surv cutpoint and surv categorise functions
##for all gene sets- each line is a different group of gene sets
#also includes survival info- time and event, and tumour type
opt1<-read_fst("h_opt.fst")
opt2<-read_fst("c2_opt.fst")
opt3<-read_fst("c3_opt.fst")
opt4<-read_fst("c6_opt.fst")
opt5<-read_fst("c7_opt.fst")
#put all the grouped gene sets together
data_opt<-list(opt1,opt2,opt3,opt4,opt5)
##name the list
names(data_opt)<-c("Hallmark","C2","C3","C6","C7")
#remove the individual ones now theyre stuck together
rm(opt1, opt2, opt3, opt4, opt5)



##these contain 2 things for each gene set
#1- for each patient high or low based upon median expression
#2 for each patient the average gene set expression
p<-read_fst("h.fst")
q<-read_fst("c2.fst")
m<-read_fst("c3.fst")
n<-read_fst("c6.fst")
o<-read_fst("c7.fst")
#put together
dataa50<-list(p,q,m,n,o)
##name the list
names(dataa50)<-c("Hallmark","C2","C3","C6","C7")
#remove originals
rm(p,q,m,n,o)

##this is a table of the cutpoint used for the optimal cut off
cuts<-readRDS("cutpoints.RDS")
names(cuts)<-c("Hallmark","C2","C3","C6","C7")



##options of groups of pathways fo later
pathway_choices<-c("Hallmark","C2","C3","C7","C6")
pathway_choices<-as.data.frame(pathway_choices)

##this is a list of all pathways and what group they go in- provides the options for later
paths<-read.csv("paths.csv")
pathws<-paths$pathway

##this is all the pathways and the genes in them
all_paths<-readRDS("all_paths.RDS")
names(all_paths)<-c("Hallmark","C2","C3","C6","C7")







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
  #CC-GSEA CC-GSEA page box
  tags$head(tags$style("
                     
                     #mybox6{height:800px !important;}
                 
                                 ")),
  
  ##GS-Corr GS-Corr page box
  tags$head(tags$style("
                     
                     #mybox7{height:775px !important;}
                 
                
                     ")),
  
  
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
                                      label = "GSEA Dataset",
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
                         ##free text for info 
                          HTML("<p>For more information about tumour type abbreviations please click <a href='https://gdc.cancer.gov/resources-tcga-users/tcga-code-tables/tcga-study-abbreviations'>here</a></p> 
         <p>For information about specifc gene sets please click <a href='https://www.gsea-msigdb.org/gsea/msigdb/genesets.jsp'>here</a></p>")
                        ),
                        
                        # main panels
                        mainPanel(
                          tabsetPanel(
                            
                            #median page
                            tabPanel("Median split",
                                     
                                     ##kaplan plot
                                     plotOutput("fifty_plot1"),
                                     ##download
                                     downloadButton('downloadPlotsurvfif', 'Download Survival Plot'),
                                     br(),
                                     
                                     ##histogram
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
                                    ##download
                                     downloadButton('downloadPlotfifhist','Download Average Expression Plot'),
                                     
                                     HTML("<p>The histogram shows the mean gene expression of each tumour after scaling. The red line shows the median cut off point used to separate high and low tumours for survival analysis.</p>"),
                                     
                                    ##genes in pathway table 
                                     fluidRow(
                                       column(
                                         dataTableOutput("Genes2"),width = 4
                                       )
                                     ),
                                     downloadButton('downloadgenefif'),
                                     
                            ),
                            
                            #optimal split
                            tabPanel("Optimal Split",
                                     
                                     #kaplan plot
                                     
                                     
                                     plotOutput("opt_plot1"),
                                     
                                     downloadButton('downloadPlotopt_surv', 'Download Survival Plot'),
                                     br(),
                                     
                                     ##histogram
                                     plotOutput("ave_hist_cut"),
                                     ##labels underneath
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
                                     
                                     ##genes in pathway
                                     fluidRow(
                                       column(
                                         dataTableOutput("Genes"),width = 4
                                       )
                                     ),
                                     downloadButton('downloadgeneopt'),
                                     
                                     
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
  
  
  

  
  ##optimal plots
  opt<-reactive({
    
    
    data=data_opt
    ##select pathways
    pathw<-as.data.frame(pathws)
    pathway<-pathw%>%
      filter(pathws==input$gene)
    pathway<-as.character(pathway)
    
    
    
    
    to_use<-data_opt[[input$surv_datset]]
    ##select tumour type
    tumour_type<-to_use%>%
      filter(tumtype==input$ttumour)
    ##gather correct data
    tumour_type<-tumour_type%>%
      select(OS.time,OS,pathway)
    colnames(tumour_type)<-c("OS.time","OS","Gene_Set")#
    
    return(tumour_type)
    
    
  })
  
  
  ##plot survival plot
  output$opt_plot1<-renderPlot({
    fit=survfit(Surv(OS.time,OS)~Gene_Set,data=opt())  
    ggsurvplot(fit,data=opt(),pval=TRUE,title = paste0(input$ttumour," ",input$gene),xlab="Time (Days)")
    
  })
  
  ###############################histogram plot optimal
  hist<-reactive({
    
    data=dataa50
    ##select pathway
    pathw<-as.data.frame(pathws)
    pathway<-pathw%>%
      filter(pathws==input$gene)
    pathway<-as.character(pathway)
    pathway2<-paste0(pathway,"_value")
    
    
    to_use<-dataa50[[input$surv_datset]]
    
    ##select tumour type
    tumour_type<-to_use%>%
      filter(tumtype==input$ttumour)%>%
      select(pathway2)
    
    colnames(tumour_type)<-"pathway"
    
    ##select correct cutpoint
    cuts_to_use<-cuts[[input$surv_datset]]
    tumour_type_cuts<-cuts_to_use%>%
      filter(tumour_type==input$ttumour)%>%
      filter(nm %in% pathway)
    
    cutpt<-tumour_type_cuts$pt
    cutpt<-as.numeric(cutpt)
    
    ##make histogram
    toplot<-ggplot(tumour_type,aes(x=pathway))+
      geom_histogram(bins=50)+
      xlab("Average Pathway Expression")+
      theme_bw()+
      ggtitle ( paste0(input$ttumour," ",input$gene))+
      geom_vline(xintercept=cutpt,linetype="dashed",color="red")
    print(toplot)
    
  })
  
  
  ##output histogram
  output$ave_hist_cut<-renderPlot({hist()})
  
  ##make genes in pathway table
  genetab<-reactive({
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
  output$Genes<-renderDataTable({
    genetab()
  })
  
  ################################################################################
  ###############################################################################
  ##########5050 plot
  fifty<-reactive({
    
    
    data=dataa50
    ##select pathway
    pathw<-as.data.frame(pathws)
    pathway<-pathw%>%
      filter(pathws==input$gene)
    pathway<-as.character(pathway)
    
    
    
    
    to_use<-dataa50[[input$surv_datset]]
    ##select tumour type
    tumour_type<-to_use%>%
      filter(tumtype==input$ttumour)
    ##gather survival things
    tumour_type<-tumour_type%>%
      select(OS.time,OS,pathway)
    colnames(tumour_type)<-c("OS.time","OS","Gene_Set")#
    
    return(tumour_type)
    
    
  })
  
  
  ##plot kaplan plot for median
  output$fifty_plot1<-renderPlot({
    fit=survfit(Surv(OS.time,OS)~Gene_Set,data=fifty())  
    ggsurvplot(fit,data=fifty(),pval=TRUE,title = paste0(input$ttumour," ",input$gene),xlab="Time (Days)")
    
  })
  
  ###############################histogram plot for median
  hist50<-reactive({
    
    data=dataa50
    ##select pathway
    pathw<-as.data.frame(pathws)
    pathway<-pathw%>%
      filter(pathws==input$gene)
    pathway<-as.character(pathway)
    pathway2<-paste0(pathway,"_value")
    
    
    to_use<-dataa50[[input$surv_datset]]
    
    ##select tumour type
    tumour_type<-to_use%>%
      filter(tumtype==input$ttumour)%>%
      select(pathway2)
    
    colnames(tumour_type)<-"pathway"
    
    ##make histogram
    toplot<-ggplot(tumour_type,aes(x=pathway))+
      geom_histogram(bins=50)+
      xlab("Average Pathway Expression")+
      theme_bw()+
      ggtitle ( paste0(input$ttumour," ",input$gene))+
      geom_vline(xintercept=(median(tumour_type$pathway)),linetype="dashed",color="red")
    print(toplot)
    
  })
  
  
  ##plot histogram
  output$ave_hist<-renderPlot({hist50()})
  
  ##genes in pathway
  genetab50<-reactive({
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
  output$Genes2<-renderDataTable({
    genetab50()
  })
  
  
  ##download plots
  
  
  
  
  output$downloadPlotopt_surv <- downloadHandler(
    filename = function(){
      paste0(input$ttumour,"_",input$gene,"_opt_cutoff_survival_plot.svg")
    },
    content = function(file) {
      svg(file)
      #  print(d())
      print(ggsurvplot((survfit(Surv(OS.time,OS)~Gene_Set,data=opt())  ),data=opt(),pval=TRUE,title =  paste0(input$ttumour," ",input$gene)))
      dev.off()
    }
  )
  
  
  
  
  output$downloadPlotsurvfif <- downloadHandler(
    filename = function(){
      paste0(input$ttumour,"_",input$gene,"_median_cutoff_survival_plot.svg")
    },
    content = function(file) {
      svg(file)
      #  print(d())
      print(ggsurvplot((survfit(Surv(OS.time,OS)~Gene_Set,data=fifty())  ),data=fifty(),pval=TRUE,title =  paste0(input$ttumour," ",input$gene)))
      dev.off()
    }
  )
  
  
  output$downloadPlotfifhist <- downloadHandler(
    filename = function(){
      paste0(input$ttumour,"_",input$gene,"_median_cutoff_histogram.svg")
    },
    content = function(file) {
      svg(file)
      print(hist50())
      dev.off()
    }
  )
  

  
  output$downloadPlotopt_hist <- downloadHandler(
    filename = function(){
      paste0(input$ttumour,"_",input$gene,"_opt_cutoff_histogram.svg")
    },
    content = function(file) {
      svg(file)
      print(hist())
      dev.off()
    }
  )
  
  
  
  
  output$downloadgeneopt <- downloadHandler(
    filename = function(){
      paste0(input$gene,"_genes.csv")
    },
    content = function(file) {
      vroom::vroom_write(genetab(), file)
    }
  )
  
  output$downloadgenefif <- downloadHandler(
    filename = function(){
      paste0(input$gene,"_genes.csv")
    },
    content = function(file) {
      vroom::vroom_write(genetab50(), file)
    }
  )
  
  
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
  
  #############GS-Corr page GS-Corr diagram
  output$gs_corr_img2 <- renderImage({
    
    list(src = "www/GS-Corr_diagram.png",
         height=700,
         width=750)
    
  }, deleteFile = F)
  
  

}

# Create a Shiny app object ----------------------------------------------------

shinyApp(ui = ui, server = server)

