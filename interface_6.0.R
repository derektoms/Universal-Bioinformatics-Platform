# 2017-07-31

########################################
#$#$#$#$#$#$    HEADER     $#$#$#$#$#$#$
########################################
# Bioinformatics packages installed via biocLite:
#source("https://bioconductor.org/biocLite.R")
#biocLite(c('limma','annotate','genefilter','ComplexHeatmap','pheatmap','cowplot','GEOmetadb','mouse4302.db','hgu133plus2.db'))

library(limma)
library(annotate)
library(genefilter)
library(ComplexHeatmap)
library(pheatmap)
library(cowplot)
library(GEOmetadb)

#biocLite(c('MergeMaid','GEOquery','inSilicoMerging','affy','sva','Rtsne','metaArray','testthat'))
library(MergeMaid)
library(GEOquery)
library(testthat)
library(metaArray)
#library(inSilicoMerging) not available for this version of R
library(Rtsne)
library(sva)
library(affy)

# Microarray platform annotations:
# Equivalent human platform is GPL570 with 127 514 samples
# HG-U133_Plus_2] Affymetrix Human Genome U133 Plus 2.0 Array

#library(mouse4302.db) 
#library(hgu133plus2.db)  Both aren't available

# Other packages:
#install.packages(c('dplyr','dbplyr','tidyr','ggplot2','RColorBrewer','readr','stringr','shiny','shinythemes','shinyjs','DT'))

library(dplyr)
library(dbplyr)
library(tidyr)
library(ggplot2)
library(RColorBrewer)
library(readr)
library(stringr)
library(DT)
library(shiny)
library(shinythemes)
library(shinyjs)

########################################
#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$
########################################

# Trying out getting and analyzing the data from GEO

if(!file.exists('GEOmetadb.sqlite')) getSQLiteFile()

db = src_sqlite("GEOmetadb.sqlite")
src_tbls(db)
gse = tbl(db, 'gse')
gse_gpl = tbl(db, 'gse_gpl')
gpl = tbl(db, 'gpl') 
gsm = tbl(db, 'gsm')
gse_gsm = tbl(db, 'gse_gsm')

# fill each row with 'not assigned' label, until they are placed in groups
# gse_selected$category <- rep("Not yet assigned", nrow(gse_selected))


########################################
#$#$#$#$#$#$    Shiny App  $#$#$#$#$#$#$
########################################

## Javascript
jscode <- '
$(function() {
  var $els = $("[data-proxy-click]");
  $.each(
    $els,
    function(idx, el) {
      var $el = $(el);
      var $proxy = $("#" + $el.data("proxyClick"));
      $el.keydown(function (e) {
        if (e.keyCode == 13) {
          $proxy.click();
        }
      });
    }
  );
});
'

## UI
ui <- fluidPage(
  tags$head(tags$script(HTML(jscode))),
  #creation of a navigation bar and mulitple pages
  navbarPage("Bioinformatics Software",
             tabPanel("Search for GEO data series (GSE)",  
                      #search GSE, and select which to include
                      helpText("After searching, click on the second tab to proceed to the next page"),
                      radioButtons("gplSelection", "Choose species:", choices = c("Mouse (GPL1260)" = "mouse", "Human (GPL570)" = "human")),
                      textOutput("gplSelection"),
                      actionButton("Search", "Search"),
                      tagAppendAttributes(
                        textInput("Key", "Enter search terms, separated by commas", value = ""),
                        `data-proxy-click` = "Search"
                      )  
             ),
             tabPanel("Select GEO data series (GSE)", uiOutput("page1"), 
                      helpText("Highlight the desired search results (GSE) and click 'Retrieve GSM' to proceed"),
                      actionButton("GSE_GSM", "Retrieve GSM"),
                      helpText("Do not click 'finish' until all selections have been made. 
                               This button removes the unselected rows and generates a new table on the next page."),
                      DT::dataTableOutput("filteredgse"),
                      tableOutput("GSEtoGSMlist")
                      
             ),
             tabPanel("Define categories for GEO samples (GSM)", uiOutput("page2"), 
                      helpText("Define the categories that you wish to compare. 
                               After this is complete, click on the third tab to proceed to the next page"),
                      textInput("cat1", "Define Category 1"),
                      textInput("cat2", "Define Category 2"),
                      textInput("cat3", "Define Category 3")
             ),
             tabPanel("Assign samples to categories", uiOutput("page3"), 
                      helpText("Highlight the desired search results and click 'assign' to assign them to the specificed category"),
                      actionButton("Assign", "Assign Categories"),
                      actionButton("Remove", "Finalize selections and remove not included"),
                      helpText("Do not click 'finish' until all selections have been made. 
                               This button removes the unselected rows and generates a new table on the next page."),
                      DT::dataTableOutput("gsm_table")
             ),
             tabPanel("Selection details", uiOutput("page4"), 
                      tableOutput("finishedtable")
             )
  )
)

## SERVER
server <- function(input, output,session) {
  ## Change GPL
  #gplSelection <- switch(input$gplSelection,
  #                 mouse = 'GPL1261',
  #                 human = 'GPL570')
  output$gplSelection <- renderText({
    paste("You chose", input$gplSelection)
  })
  
  ## Convert SQLite to data frames
  as.data.frame.DataTable(gse) -> gse.df
  as.data.frame.DataTable(gse_gsm) -> gse_gsm.df
  #as.data.frame.DataTable(gsm) -> gsm.df
  gse_to_filter <- data.frame(gse="")
  
  ## Search functions
  Totalchar <- eventReactive(input$Search, {nchar(input$Key)})
  Commas <- eventReactive(input$Search, {which(strsplit(input$Key, "")[[1]]==",")})
  Ncommas <- eventReactive(input$Search, {length(Commas())})
  Commasstart <- eventReactive(input$Search, {Commas() + 1})
  Commasend <- eventReactive(input$Search, {Commas() - 1})
  
  Searchterms <- eventReactive(input$Search, {
    substring(input$Key, c(1, Commasstart()), c(Commasend(), Totalchar()))
  })
  
  filtered_gse <- eventReactive(input$Search, {
    dplyr::filter(gse.df, str_detect(gse.df$title, Searchterms()))
  })
  
  ## Collect samples to use (GSE - GSM)
    # List of the GSM associated with the selected GSE
  
  gse_to_keep <- eventReactive(input$GSE_GSM, {
    filtered_gse()[input$filteredgse_rows_selected,c(3,7)]
  })
   
  ## Assign categories to each sample (GSM)
    # test
  gsm_annotated <- eventReactive(input$GSE_GSM, {
    dplyr::filter(gsm.df,series_id %in% gse_to_keep()$gse)
  })
  
  #import dataframe as reactive
  #for some reason the following line is not needed:
  rows <- reactiveValues() 
  eventReactive(input$Assign, {
    #s <- eventReactive(input$Lock, input$gsm_table_rows_selected)
    if (input$Assign == 1) {
      gsm_table <- filteredgsm()
      gsm_selected[input$gsm_table_rows_selected,"category"] <- input$selection
      rows$df <- gsm_selected
      gsm_selected <<- rows$df # '<<-' is necessary to get this to the enclosing environment
    }
    else
    {
      gsm_selected[input$gsm_table_rows_selected,"category"] <- input$selection
      rows$df <- gsm_selected
      gsm_selected <<- rows$df
    }
  })
  
  finishedtable <- eventReactive(input$Remove, {
    dplyr::filter(rows$df, category %in% c(input$cat1, input$cat2, input$cat3))
  })
 
  ## Outputs
  output$page3 <- renderUI(
    fluidRow(
      column(3,
             selectInput("selection", "Select a Category",
                         c("category1" <- {input$cat1},
                           "category2" <- {input$cat2},
                           "category3" <- {input$cat3},
                           "category4" <- "Not included"))
      )
    )
  )
  
  output$filteredgse <- DT::renderDataTable({
 #   if (input$GSE_GSM == 0)
      return (filtered_gse()[,c(1,2,7)]) ###
#    else ## breaks when there is no search results
#      return (rows$df)
  }, options=list(searching=TRUE, pageLength=50))
 
  output$GSEtoGSMlist <- renderTable(
    if (input$GSE_GSM == 0)
      return ()
    else
      return (filter(gse_gsm.df,gse %in% gse_to_keep()$gse)))## works
  
 # output$gsm_table <- DT::renderDataTable({
 #   filter(gsm.df,series_id %in% gse_to_keep()$gse)}, options=list(searching=TRUE))
    output$gsm_table <- DT::renderDataTable({gsm_annotated()}, options=list(searching=TRUE))

  output$finishedtable <- renderTable({
    if (input$Remove == 0)
      return (filteredgsm()[,c(1,2,7,19)])
    else
      return (finishedtable()[,c(1,2,3,7,19)])})
  
  ## Kill shinyApp when session closes
  session$onSessionEnded(stopApp)

}

shinyApp(ui, server)

