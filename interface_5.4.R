# 2017-06-19

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

########################################
#$#$#$#$#$#$    Shiny App  $#$#$#$#$#$#$
########################################

# fill each row with 'not assigned' label, until they are placed in groups
as.data.frame.DataTable(gse) -> gsm_selected
gsm_selected$category <- rep("Not yet assigned", nrow(gsm_selected))


ui <- fluidPage(
  #creation of a navigation bar and mulitple pages
  navbarPage("Bioinformatics Software",
             tabPanel("Search for GEO data series (GSE)",  
                      #search GSE, and select which to include
                      helpText("After searching, click on the second tab to proceed to the next page"),
                      textInput("Key", "Enter search terms, separated by commas", value = ""),
                      actionButton("Search", "Search")
             ),
             tabPanel("Select GEO data series (GSE)", uiOutput("page1"), 
                      helpText("Highlight the desired search results (GSE) and click 'Retrieve GSM' to proceed"),
                      actionButton("GSE-GSM", "Retrieve GSM"),
                      helpText("Do not click 'finish' until all selections have been made. 
                               This button removes the unselected rows and generates a new table on the next page."),
                      DT::dataTableOutput("gse-gsm_table")
             ),
             tabPanel("Define categories for GEO samples (GSM)", uiOutput("page2"), 
                      helpText("Define the categories that you wish to compare. 
                               After this is complete, click on the third tab to proceed to the next page"),
                      textInput("cat1", "Define Category 1"),
                      textInput("cat2", "Define Category 2"),
                      textInput("cat3", "Define Category 3")
             ),
             
             # changed the format of this slightly to accomodate the DT package (for row selection)
             tabPanel("Assign samples to categories", uiOutput("page3"), 
                      helpText("Highlight the desired search results and click 'assign' to assign them to the specificed category"),
                      actionButton("Lock", "Assign Categories"),
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

server <- function(input, output) {
  
  Totalchar <- eventReactive(input$Search, {nchar(input$Key)})
  
  Commas <- eventReactive(input$Search, {which(strsplit(input$Key, "")[[1]]==",")})
  
  Ncommas <- eventReactive(input$Search, {length(Commas())})
  
  Commasstart <- eventReactive(input$Search, {Commas() + 1})
  
  Commasend <- eventReactive(input$Search, {Commas() - 1})
  
  Searchterms <- eventReactive(input$Search, {
    substring(input$Key, c(1, Commasstart()), c(Commasend(), Totalchar()))
  })
  
  # Searchterms2 <- eventReactive(input$Search, {
  #   paste0("//b", Searchterms(), "//b")
  # })
  # 
  # Searchterms3 <- eventReactive(input$Search, {
  #   gsub(",//b", "//b|", Searchterms2())
  # })
  
  filteredgsm <- eventReactive(input$Search, {dplyr::filter(gsm_selected, str_detect(gsm_selected$title, Searchterms()))})
  
  # List of the GSM associated with the selected GSE
  gsegsm2 <- eventReactive(input$Search, {filter(gse_gsm,gse %in% filteredgsm()$gse)}) # list of series and associated samples
  gsm2 <- eventReactive(input$Search, {filter(gsm,series_id %in% filteredgsm()$gse)}) # detailed sample information
  
  #import dataframe as reactive
  #for some reason this line is not needed
  rows <- reactiveValues() 
  observeEvent(input$Lock, {
    #s <- eventReactive(input$Lock, input$gsm_table_rows_selected)
    if (input$Lock == 1) {
      gsm_selected <- filteredgsm()
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
  
  #output$Commas <- renderText(Commas())
  #output$Ncommas <- renderText(Ncommas())
  #output$Searchterms2 <- renderText(Searchterms())
  
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
  
  output$gse-gsm_table <- DT::renderDataTable({
    if (input$GSE-GSM == 0)
      return (filteredgsm()[,c(1,2,7,19)]) ###
    else ## breaks when there is no search results
      return (rows$df[,c(1,2,7,19)])}, options=list(searching=TRUE)) ###
  
  output$gsm_table <- DT::renderDataTable({
    if (input$Lock == 0)
      return (filteredgsm()[,c(1,2,7,19)])
    else ## breaks when there is no search results
      return (rows$df[,c(1,2,7,19)])}, options=list(searching=TRUE))

  output$finishedtable <- renderTable({
    if (input$Remove == 0)
      return (filteredgsm()[,c(1,2,7,19)])
    else
      return (finishedtable()[,c(1,2,3,7,19)])})
  
}

shinyApp(ui, server)

