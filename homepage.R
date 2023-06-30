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


##make customisable loading page
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
  useShinyjs(),
  inlineCSS(appCSS),

  
  
  #Custom CSS
  #make 3 boxes for my 3 figures to go in later
  tags$head(tags$style("
                     
                     #mybox{height:800px !important;}
                 
                                 ")),
  
  tags$head(tags$style("
                     
                     #mybox2{height:1300px !important;}
                 
                                     ")),
  
  tags$head(tags$style("
                     
                     #mybox3{height:775px !important;}
                 
                
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
  
  # Loading message
  div(
    id = "loading-content",
    h2("Loading... this may take a minute or two...")
  ),
  
  
 
  
  navbarPage("GS-TCGA",
             theme = bs_theme(bootswatch = "sandstone"),
             
             tabPanel(title = "Home",
                      
      
   ##intro section                   
                      h4(strong("GS-TCGA: Gene Set-based Analysis of The Cancer Genome Atlas")),
                      p(style="text-align: justify; font-size = 50px",
                        "This website is a resource designed to use Gene Set Enrichment Analysis
                        to provide novel interpretation of the wealth of data contained in The Cancer Genome Atlas."),
                      p(style="text-align: justify; font-size = 50px",strong("Gene Set Survival Analysis: "),
                        "GS-Surv allows the user to investigate how the expression of genes in their specified gene set relates to overall survival in patient data. "),
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
             tabPanel(HTML("</a></li><li><a href=\"https://gs-tcga.shinyapps.io/GS-Surv\">GS-Surv")),
             tabPanel(HTML("</a></li><li><a href=\"https://gs-tcga.shinyapps.io/CC-GSEA\">CC-GSEA")),
             tabPanel(HTML("</a></li><li><a href=\"https://gs-tcga.shinyapps.io/gene_set_correlates\">GS-Corr")),
             tabPanel(HTML("</a></li><li><a href=\"https://gs-tcga.shinyapps.io/upload\">GS-Surv (Custom)")),
             
             
             )
  
  
  
  )

  

  


server <- function(input, output, session) {
  
  ###########homepage
 
  ##this is putting the images in
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
  
  # Hide the loading message when the rest of the server function has executed
  hide(id = "loading-content", anim = TRUE, animType = "fade")    
  show("app-content")
  
}

# Create a Shiny app object ----------------------------------------------------

shinyApp(ui = ui, server = server)

