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
library(fst)




#######UI

ui <- fluidPage(
  useShinydashboard(),
  useShinyjs(),


  
  
  #Custom CSS
  #make 3 boxes for my 3 figures to go in later
  ##CCGSEA home page box
  tags$head(tags$style("
                     
                     #mybox{height:800px !important;}
                 
                                 ")),
  ##GS-Surv home page box
  tags$head(tags$style("
                     
                     #mybox2{height:1300px !important;}
                 
                                     ")),
  #GS-Corr home page box
  tags$head(tags$style("
                     
                     #mybox3{height:775px !important;}
                 
                
                     ")),
  
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
  
  ##make a footer which I fill in the page section
  tags$footer(tags$style("
                     
                     #footer {
                    display: flex;
                    justify-content: space-between;
                    background:#403c3c;
                    color:white;
                    height:50px;
                    align-items:center;
                    
               }
                
                     
                     
                     ")),
  

  
  
 
  #########################home page###########################################
  navbarPage("GS-TCGA",
             theme = bs_theme(bootswatch = "sandstone"),
             
             tabPanel(title = "Home",
                      
      
   ##intro section                   
                      h4(strong("GS-TCGA: Gene Set-based Analysis of The Cancer Genome Atlas")),
                      p(style="text-align: justify; font-size = 50px",
                        "This website is a resource designed to use Gene Set Enrichment Analysis
                        to provide novel interpretation of the wealth of data contained in The Cancer Genome Atlas."),
                      p(style="text-align: justify; font-size = 50px",strong("Gene Set Survival Analysis: "),
                        "GS-Surv allows the user to investigate how the average expression of genes in a specified gene set relates to overall survival in patient data. "),
                      p(style="text-align: justify; font-size = 50px",strong("Co-Correlative Gene Set Enrichment Analysis: "),
                        "CC-GSEA allows generation of novel hypotheses of gene function through performing GSEA on co-correlated genes."),
                      p(style="text-align: justify; font-size = 50px",strong("Gene Set Correlative Analysis: "),
                        "GS-Corr calculates the average expression of a gene set and correlates this with individual genes."),
                      p(style="text-align: justify; font-size = 50px",strong("GS-Surv (Custom): "),
                        "This function allows you to upload your own gene set for GS-Surv survival analysis. "),
                      br(),
                      hr(),
                      
                      
        ##GS-Surv diagram and info     
                      box(
                        title = h4("GS-Surv"),
                        width = 12,
                        "The survival analysis tool allows the user to see the relationship of their selected gene set with survival of patients in the TCGA in a Kaplan-Meier plot. The mean expression of genes in the selected gene set and tumour type is calculated, 
               this is then used to generate a Kaplan-Meier plot either splitting data at the median or at the optimal point for discrimination between high and low gene set expression. A histogram of distributions of gene expression and the cut-off is also shown. Users can upload custom gene sets for analysis on the GS-Surv (Custom) page.",
                        br(),
                        id = "mybox2",
                        imageOutput("gs_surv_img")
                      ),
                      
     ##CC-GSEA diagram and info                 
                        box(
               title = h4("CC-GSEA"),
               width = 12,
               "Genes involved in similar cellular processes/pathways may correlate with one another. Co-correlative Gene Set Enrichment Analysis aims to suggest novel gene functions through performing GSEA on a ranked list of genes which correlate with the gene of interest.
               The user can select a tumour type and gene of interest, then the expression of the gene of interest in the TCGA is correlated with all other genes, forming the ranked list used as an input to GSEA. 
               GSEA  analysis outputs a table of gene sets along with their enrichment scores and other statistics. Gene sets with a high normalised enrichment score are enriched in genes that correlate positively with the gene of interest. The user can then select a gene set from this list to explore in enrichment plots and export summary statistics.",
               id = "mybox",
               br(),
               
               imageOutput("cc_gsea_img")
             ),
             
             
       ##GS-Corr diagram and info      
             box(
               title = h4("GS-Corr"),
               width = 12,
               "The Gene Set Correlates tool aims to suggest novel gene functions through correlating the expression of genes in a pre defined gene set with expression of all other genes. 
               The user can select a tumour type and gene set of interest, then the average expression of genes in the gene set of interest is calculated and correlated with expression of all other genes in the TCGA.",
               id = "mybox3",
               br(),
               
               imageOutput("gs_corr_img")
             ),
             
             HTML("<p><font size='2'> The data used here are in whole or part based upon data generated by the TCGA Research Network: https://www.cancer.gov/tcga <br> We acknowledge our use of the gene set enrichment analysis, GSEA software, and Molecular Signature Database (MSigDB) (Subramanian, Tamayo, et al. (2005), PNAS 102, 15545-15550, http://www.broad.mit.edu/gsea/,  Xie et al. (2005), Nature 434).</font></p>"),
             

      ##footer with names and things       
             
         HTML('<div class="row">
               <div id="footer">
               <div >Tarrion Baird</div>
               <div >Roychoudhuri Laboratory</div>
               <div >University of Cambridge</div>
               </div>
               </div>
             ')
             ),
   
   
   
     ##the other website pages    
   
   
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
  
  ###########homepage
 
  ##this is putting the images in the home page
  output$gs_surv_img <- renderImage({
    
    list(src = "www/gs_surv_diagram.png",
    height=1200,
    width=800)
    
  }, deleteFile = F)
  
  output$cc_gsea_img <- renderImage({
    
    list(src = "www/ccgsea_diagram.png",
         height=700,
         width=750)
    
  }, deleteFile = F)
  
  output$gs_corr_img <- renderImage({
    
    list(src = "www/GS-Corr_diagram.png",
         height=700,
         width=750)
    
  }, deleteFile = F)
  
  output$box_state <- renderText("CC-GSEA")
  output$box_state2 <- renderText("GS-Surv")
  output$box_state3 <- renderText("GS-Corr")
  

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
  
  #############GS-Corr page GS-Corr diagram
  output$gs_corr_img2 <- renderImage({
    
    list(src = "www/GS-Corr_diagram.png",
         height=700,
         width=750)
    
  }, deleteFile = F)
  
  
}

# Create a Shiny app object ----------------------------------------------------

shinyApp(ui = ui, server = server)

