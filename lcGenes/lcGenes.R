# Sys.setlocale("LC_ALL", "C")

library(shiny)
library(DT)
# options(warn=2)  # warning are turned into errors
library(shinyBS)
library(yaml)
library(later)
#----------------------------------------------------------------------------------------------------
printf = function (...) print (noquote (sprintf (...)))
#----------------------------------------------------------------------------------------------------
# tooltips <- yaml.load_file("tooltips.yaml")
# for(i in 1:length(tooltips)) tooltips[[i]]$text <- paste(tooltips[[i]]$text, collapse=" ")
#printf("length of tooltips read: %d", length(tooltips))
#print(tooltips[1])
#print(tooltips[2])

# tt <- list(longevityTab=list(selector="DataTables_Table_0_wrapper > div > div:nth-child(3) > div > div.dataTables_scrollHead > div > table > thead > tr > th:nth-child(2)",
#                             text=paste0('[+-] causal evidence<br>',
#                                        '[+-]? correlational only<br>',
#                                         '+/- direction unknown, causal<br>',
#                                         '+/-? unknown, correlational')))


tbl.summary <- get(load("tbl.summary.1010x6.RData"))
tbl.summary <- tbl.summary[, c("Gene", "Longevity", "Feature", "Function", "PMID")]
colnames(tbl.summary)[3] <- "Longevity Feature"

tbl.orthologsBySpecies <- get(load("tbl.orthologsBySpecies.9999x4.RData"))
#----------------------------------------------------------------------------------------------------
# https://docs.google.com/document/d/1Qwj9-vj8Q7b0GWLCs5UmHirQqGkejMex5w11E-CXLvU/edit?usp=sharing?embedded=true"
#----------------------------------------------------------------------------------------------------
introTab <- function()
{
   tabPanel(title="Introduction", value="introductionTab", includeHTML("introduction.html"))
}
#----------------------------------------------------------------------------------------------------
mainGeneTab <- function()
{
   destinations <- list(" " = "goNowhere",
                        "All Tabs" = "allTabs",
                        "GeneCards" = "genecards",
                        "HomoloGene" = "homologene",
                        "PubMed (curated)" = "pubmedCurated",
                        "PubMed (gene + longevity)" = "pubmedLongevity",
                        "Orthologs" = "orthologs"
                        )

   tab <- tabPanel(title="Gene Table", value="byGeneTab",
      wellPanel(
          fluidRow(
             #column(width=3, offset=0, selectInput("selectDestination", label="Table selection goes to",
             #                                      choices=destinations, width="400px")),
             column(width=8, offset=0, h4("Click on table row for your gene of interest, populate all tabs")),
             column(width=1, offset=3, actionButton(label="?", inputId="geneTabHelpButton", style="float: right")),
             style="margin-bottom: 20px"),
           DTOutput("geneTable"), style="margin-top:5px;")
      ) # tabPanel
   tab

} # mainGeneTab
#----------------------------------------------------------------------------------------------------
ui <- fluidPage(

   includeCSS("lcGenes.css"),
   tags$head(tags$script(type="module", src="https://unpkg.com/x-frame-bypass")),
   titlePanel("Longevity-associated Genes & Orthologs"),

   tabsetPanel(type="tabs", id="lcGenesTabs",
               #introTab(),
               mainGeneTab(),
               tabPanel(title="GeneCards", value="geneCardsTab",
                        wellPanel(htmlOutput("geneCardsDisplay"))),
               tabPanel(title="HomoloGene", value="homoloGeneTab",
                        wellPanel(htmlOutput("homologeneDisplay"))),
               tabPanel(title="Curated PubMed", value="pubmedTab",
                        wellPanel(htmlOutput("pubmedDisplay"))),
               tabPanel(title="PubMed Gene+longevity", value="pubmedDeNovoTab",
                        wellPanel(htmlOutput("pubmedDeNovoDisplay"))),
               tabPanel(title="Orthologs", value="orthologsTab",
                        wellPanel(DTOutput("orthologsTableDisplay"), style="margin-top:5px;")),
               tabPanel(title="Notes & Comments", value="notesAndCommentsTab",
                        wellPanel(htmlOutput("notesAndCommentsDisplay")))
               # tabPanel("Help", includeHTML("help.html"))
               ),
    style="margin: 10px; margin-top: 10px; margin-bottom: 50px;"

) # fluidPage
#----------------------------------------------------------------------------------------------------
server <- function(session, input, output) {

   reactiveInputs <- reactiveValues(gene="",
                                    destination="",
                                    pmid="",
                                    homologeneQuery=FALSE,
                                    geneCardsQuery=FALSE,
                                    pubmedQuery=FALSE,
                                    pubmedDeNovoQuery=FALSE,
                                    orthologsQuery=FALSE,
                                    notesAndCommentsQuery=FALSE
                                    )

   observeEvent(input$geneTabHelpButton, {
      showModal(modalDialog(includeHTML("geneTableHelp.html"),
                            title = "Gene Table Help",
                            easyClose = TRUE,
                            footer = NULL
                            ))
      })

   output$geneTable <- DT::renderDataTable({
       DT::datatable(tbl.summary,
                     rownames=FALSE,
                     options=list(dom='<lfip<t>>',
                                  scrollX=TRUE,
                                  autoWidth=TRUE,
                                  columnDefs=list(list(width="10%", targets=c(0,1)),
                                                  list(width="40%", targets=c(2,3)),
                                                  list(width="10%", targets=4)),
                                  lengthMenu = c(3,5,10,50),
                                  pageLength = 5,
                                  paging=TRUE),
                     selection="single")

      })

   observeEvent(input$geneTable_rows_selected, {

      reactiveInputs$geneCardsQuery <- FALSE
      reactiveInputs$homologeneQuery <- FALSE
      reactiveInputs$pubmedQuery <- FALSE
      reactiveInputs$pubmedDeNovoQuery <- FALSE
      reactiveInputs$orthologsQuery <- FALSE

      selectedTableRow <- isolate(input$geneTable_rows_selected)
      gene <- tbl.summary[selectedTableRow, "Gene"]
      printf("--- geneTable_rows_selected: %s", gene);
      pubmedIds <- tbl.summary[selectedTableRow, "PMID"]
         # but see https://www.html5rocks.com/en/tutorials/security/sandboxed-iframes/
         # https://github.com/niutech/x-frame-bypass
      #if(!grepl(",", pubmedIds))
      #    pubmedIds <- sprintf("%s,%s", pubmedIds, pubmedIds)
      printf("pubmedIDs: %s", pubmedIds)
      #destination <- isolate(input$selectDestination)

      destination <- "allTabs"

      printf("destination: %s", destination)

      tabName <- "byGeneTab"

      if(destination == "allTabs"){
         printf("allTabs requested")
         reactiveInputs$geneCardsQuery <- TRUE
         reactiveInputs$homologeneQuery <- TRUE
         reactiveInputs$pubmedQuery <- TRUE
         reactiveInputs$pubmedDeNovoQuery <- TRUE
         reactiveInputs$orthologsQuery <- TRUE
         }

      if(destination == "genecards"){
         printf("destination: genecards")
         tabName <- "geneCardsTab"
         reactiveInputs$geneCardsQuery <- TRUE
         }

      if(destination == "homologene"){
         printf("destination: homologene")
         tabName <- "homoloGeneTab"
         reactiveInputs$homologeneQuery <- TRUE
         }

      if(destination == "pubmedCurated"){
         printf("destination: pubmedCurated")
         tabName <- "pubmedTab"
         reactiveInputs$pubmedQuery <- TRUE
         }

      if(destination == "pubmedLongevity"){
         printf("destination: pubmedLongevity")
         tabName <- "pubmedDeNovoTab"
         reactiveInputs$pubmedDeNovoQuery <- TRUE
         }

      if(destination == "orthologs"){
         printf("destination: orthologs")
         tabName <- "orthologsTab"
         reactiveInputs$orthologsQuery <- TRUE
         }

      if(destination == "Notes & Comments"){
         printf("destination: notes & comments")
         tabName <- "notesAndCommentsTab"
         reactiveInputs$notesAndCommentsQuery <- TRUE
         }

      updateTabsetPanel(session, "lcGenesTabs", selected=tabName)
      printf("send %s to %s", gene, destination)
      reactiveInputs$gene <- gene
      reactiveInputs$destination <- destination
      reactiveInputs$pmid <- pubmedIds
      }) # observe row selection event

    output$orthologsTableDisplay <- DT::renderDataTable({
       printf(" ==== entering orthologsTableDisplay")
       if(reactiveInputs$orthologsQuery){
         goi <- isolate(reactiveInputs$gene)
         printf("--- rendering orthologsTableDisplay")
         printf("goi: %s", goi)
         printf("subset tbl.orthologsBySpecies by %s", goi)
         tbl.orthoSub <- subset(tbl.orthologsBySpecies, HsOrtholog_Gene_Symbol==goi)
         reactiveInputs$orrhtologsQuery <- FALSE
         DT::datatable(tbl.orthoSub, rownames=FALSE, options=list(paging=FALSE))
         } # if Orthologs
      }) # orthologsDisplay

    output$geneCardsDisplay <- renderUI({
       printf(" ==== entering geneCardsDisplay")
    #   if(reactiveInputs$geneCardsQuery){
    #      goi <- isolate(reactiveInputs$gene)
    #      printf(" geneCardsDisplay on %s", goi)
    #      printf("--- rendering geneCardsDisplay %s", goi)
    #      uri <- sprintf("https://www.genecards.org/cgi-bin/carddisp.pl?gene=%s", goi)
    #      printf("uri: %s", uri)
    #      htmlText <- tags$iframe(src=uri, height=1000, width="100%")
    #      # reactiveInputs$geneCardsQuery <- FALSE
    #      htmlText
    #      }
        }) # geneCardsDisplay

    output$homologeneDisplay <- renderUI({
       printf("=== entering homologeneDisplay")
       goi <- isolate(reactiveInputs$gene)
    #   if(reactiveInputs$homologeneQuery){
    #      printf("--- rendering homologeneDisplay")
    #      doi <- reactiveInputs$destination
    #      printf("goi: %s", goi)
    #      printf("doi: %s", doi)
    #      uri <- sprintf("https://www.ncbi.nlm.nih.gov/homologene/?term=%s", goi)
    #      printf("uri: %s", uri)
    #      htmlText <- tags$iframe(src=uri, height=1000, width="100%")
    #      htmlText
    #      } # HomoloGene
       }) # homologeneDisplay


    output$pubmedDisplay <- renderUI({
       printf("=== entering pubmedDisplay")
    #   if(reactiveInputs$pubmedQuery){
    #      printf("--- rendering pubmedDisplay")
    #      pmids <- reactiveInputs$pmid
    #      doi <- reactiveInputs$destination
    #      printf("pmids: %s", pmids)
    #      printf("doi: %s", doi)
    #      uri <- sprintf("https://pubmed.ncbi.nlm.nih.gov/?term=%s", pmids)
    #      printf("uri: %s", uri)
    #      htmlText <- tags$iframe(src=uri, is="x-frame-bypass", height=1000, width="100%")
    #      htmlText
    #      }
       }) # pubmedDisplay

    output$pubmedDeNovoDisplay <- renderUI({
       printf("=== entering pubmedDeNovoDisplay")
     #  goi <- isolate(reactiveInputs$gene)
     #  if(reactiveInputs$pubmedDeNovoQuery) {
     #     printf("--- rendering pubmedDeNovoDisplay")
     #     pmids <- reactiveInputs$pmid
     #     doi <- reactiveInputs$destination
     #     printf("pmids: %s", pmids)
     #     printf("doi: %s", doi)
     #     queryString <- sprintf("?term=%s+longevity", goi)
     #     #uri <- sprintf("https://www.ncbi.nlm.nih.gov/pubmed/%s", queryString)
     #     uri <- sprintf("https://pubmed.ncbi.nlm.nih.gov/%s", queryString)
     #     printf("uri: %s", uri)
     #     htmlText <- tags$iframe(src=uri, is="x-frame-bypass", height=1000, width="100%")
     #    htmlText
     #    }
       }) # pubmedDisplay

    output$notesAndCommentsDisplay <- renderUI({
       printf("--- notesAndComment")
       uri <- sprintf("https://docs.google.com/document/d/1Qwj9-vj8Q7b0GWLCs5UmHirQqGkejMex5w11E-CXLvU/edit?usp=sharing")
       printf("uri: %s", uri)
       htmlText <- tags$iframe(src=uri, height=1000, width="100%")
       htmlText
       }) # notesAndCommentsDisplay

    outputOptions(output, "geneCardsDisplay",        suspendWhenHidden = FALSE)
    outputOptions(output, "homologeneDisplay",       suspendWhenHidden = FALSE)
    outputOptions(output, "notesAndCommentsDisplay", suspendWhenHidden = FALSE)
    outputOptions(output, "pubmedDeNovoDisplay",     suspendWhenHidden = FALSE)
    outputOptions(output, "pubmedDisplay",           suspendWhenHidden = FALSE)
    outputOptions(output, "orthologsTableDisplay",   suspendWhenHidden = FALSE)

   } # server

#----------------------------------------------------------------------------------------------------
runApp(shinyApp(ui=ui, server=server), port=9003)
#shinyApp(ui=ui, server=server)




