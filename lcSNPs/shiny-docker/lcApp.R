printf <- function(...) print(noquote(sprintf(...)))
library(shiny)
library(shinyModules)
library(later)
library(RSQLite)
library(TrenaProjectHG38.generic)
tp.hg38 <- TrenaProjectHG38.generic()
#----------------------------------------------------------------------------------------------------
tbl.gh <- get(load("tbl.gh.apoe.RData"))
colnames(tbl.gh)[grep("ghid", colnames(tbl.gh))] <- "name"

tbl.snps <- get(load("tbl.snps.598.RData"))

tbl.genes <- get(load("tbl.summary.1010x6-gene-rowNames.RData"))
tbl.genes <- tbl.genes[, c("Longevity", "Feature", "Function", "PMID")]
colnames(tbl.genes)[2] <- "Longevity Feature"
#----------------------------------------------------------------------------------------------------
driver <- dbDriver("SQLite")
dbConnection <- dbConnect(driver, dbname = "sqlite.db")
template.comment <- data.frame(author="",
                               timestamp=Sys.time(),
                               entity="",
                               tags="",
                               text="",
                               stringsAsFactors=FALSE)
if(!"comments" %in% dbListTables(dbConnection))
   dbCreateTable(dbConnection, "comments", template.comment)
#----------------------------------------------------------------------------------------------------
snpsPage <- function()
{
    fluidRow(
        conditionalPanel(
            condition = "input.genesOrSnps=='snps'",
            div(dataTableUI("snpDataTable"),
                style="margin: 10px; margin-bottom: 30px; padding: 10px; border: 3px solid black; border-radius: 10px;")
            ),
        conditionalPanel(
            condition = "input.genesOrSnps=='genes'",
            div(dataTableUI("geneDataTable"),
                style="margin: 10px; margin-bottom: 30px; padding: 10px; border: 3px solid black; border-radius: 10px;")
            ),
        div(igvUI("igv"),
            style="margin: 10px; margin-bottom: 5px; padding: 10px; border: 3px solid black; border-radius: 10px;"),
        style="margin-top:5px;"
        )

} # snpsPage
#----------------------------------------------------------------------------------------------------
mainPage <- function()
{
  fluidPage(
    titlePanel("Exploring The Genetic and Molecular Bases of Human Longevity"),
    sidebarPanel(
        radioButtons("genesOrSnps", "Upper Panel", choices=c(Genes="genes", SNPs="snps", Hide="hide")),
        verbatimTextOutput("igvSelection"),
        messageBoxUI(id="messageBox.snpTable", title="snp", fontSize=16, boxWidth=120, boxHeight=23),
        messageBoxUI(id="messageBox.geneTable", title="gene", fontSize=16, boxWidth=120, boxHeight=24),
        messageBoxUI(id="messageBox.igv", title="igv", fontSize=16, boxWidth=120, boxHeight=22),
        br(),
        selectInput("genesWithSNPsSelector", "Genes w/Centarian SNPs:",
                    c("", "ACOX1", "APOC1", "APOE", "FOXO3", "GRHL1", "KNL1", "NECTIN2",
                      "SIAH1", "TCF21", "TMTC2", "TOMM40", "USP42")),
        width=2),
    mainPanel(
        snpsPage(), width=10
        )
   ) # fluidPage

} # mainPage
#----------------------------------------------------------------------------------------------------
ui <- fluidPage(
   tabsetPanel(type="tabs", id="lcGenesTabs",
          #introTab(),
          tabPanel(title="Home", mainPage()),
          tabPanel(title="GeneCards", wellPanel(iframeSearchUI(id="geneCardsSearch",   title="GeneCards"))),
          tabPanel(title="dbSNP", wellPanel(iframeSearchUI(id="dbSNPSearch",   title="dbSNP"))),
          tabPanel(title="rVarBase", wellPanel(iframeSearchUI(id="rVarBaseSearch",   title="rVarBase"))),
          tabPanel(title="ClinVar", wellPanel(iframeSearchUI(id="clinvarSearch",   title="ClinVar"))),
          tabPanel(title="Homologene", wellPanel(iframeSearchUI(id="homologeneSearch",  title="Homologene"))),
          tabPanel(title="KEGG", wellPanel(iframeSearchUI(id="keggSearch",  title="KEGG"))),
          tabPanel(title="PubMed", wellPanel(iframeSearchUI(id="pubmedSearch", title="PubMed"))),
          tabPanel(title="Google", wellPanel(iframeSearchUI(id="googleSearch", title="Google"))),
          tabPanel(title="Notes", wellPanel(commentsUI(id="comments"), title="Notes"))
          ),
    style="margin: 10px; margin-top: 10px; margin-bottom: 50px;"
   )


#----------------------------------------------------------------------------------------------------
server <- function(input, output, session){

   table.selected.snps <- reactiveVal(c("rs5117"))
   table.selected.genes <- reactiveVal("APOE")
   igv.selected.entity <- reactiveVal("rs5117")

    observeEvent(input$genesWithSNPsSelector, ignoreInit=TRUE, {
       goi <- input$genesWithSNPsSelector
       req(goi)
       dispatchOnGene(session, goi, origin="genesWithSNPsSelector")
       })

   table.selected.snps <- callModule(dataTableServer, "snpDataTable",
                                     tbl=tbl.snps,
                                     selectionPolicy="single",
                                     pageLength=10,
                                     visibleRows = reactive("all"))

   table.selected.gene <- callModule(dataTableServer, "geneDataTable",
                                     tbl=tbl.genes,
                                     selectionPolicy="single",
                                     pageLength=10,
                                     visibleRows = reactive("all"))

   igv.selected.entity <- callModule(igvServer, "igv",
                                     genome="hg38",
                                     geneModelDisplayMode="COLLAPSED",
                                     locus="chr19:44,854,808-44,940,011") #"APOE")

   observe({
      selectedGene <- table.selected.gene()
      dispatchOnGene(session, selectedGene, origin="gene.table")
      })

    observe({
       selectedSNP <- table.selected.snps()
       dispatchOnSNP(session, selectedSNP, origin="snp.table")
       })

   observe({
      selectedEntity <- igv.selected.entity()
      printf("igv.selected.entity: %s", selectedEntity)
      callModule(messageBoxServer, "messageBox.igv", newContent=reactive(selectedEntity))
      req(selectedEntity)
      printf("about to test igv selection: %s", selectedEntity)
      if(grepl("^rs", selectedEntity))
          dispatchOnSNP(session, selectedEntity, "igv")
      else
          dispatchOnGene(session, selectedEntity, "igv")
      })

  later(function(){
      tbl.bed <- tbl.snps[, c("chrom", "start", "end")]
      tbl.bed$name <- rownames(tbl.snps)
      tbl.bed$start <- tbl.bed$start - 1
      loadBedTrack(session, trackName="snps", tbl=tbl.bed, color="red");
      tbl.bed <- tbl.gh[, c("chrom", "start", "end", "combinedscore", "name")]
      loadBedGraphTrack(session, trackName="GeneHancer", tbl=tbl.bed, color="blue",
                        autoscale=FALSE, min=1, max=50)
      }, 3)

} # server
#----------------------------------------------------------------------------------------------------
dispatchOnSNP <- function(session, selectedSNP, origin)
{
  printf("--- dispatchOnSNP")

  if(origin != "igv"){
     chrom.loc <- with(tbl.snps[selectedSNP,],  sprintf("%s:%d-%d", chrom, start-10, end+10))
     showGenomicRegion(session, chrom.loc)
     }

  if(origin != "snp.table"){
    callModule(dataTableServer, "snpDataTable",
               tbl=tbl.snps,
               selectionPolicy="single",
               pageLength=10,
               visibleRows = reactive("all"),
               searchTerm=reactive(selectedSNP))
     } # !snp.table

  callModule(messageBoxServer, "messageBox.snpTable", newContent=reactive(selectedSNP))
  callModule(iframeSearchServer, "googleSearch",
             website=reactive("google"), geneSymbol=reactive(selectedSNP))

  callModule(iframeSearchServer, "dbSNPSearch",
             website=reactive("dbsnp"), geneSymbol=reactive(selectedSNP))

  callModule(iframeSearchServer, "rVarBaseSearch",
             website=reactive("rvarbase"), geneSymbol=reactive(selectedSNP))

  callModule(iframeSearchServer, "pubmedSearch",
             website=reactive("PubMed"), geneSymbol=reactive(selectedSNP))

  callModule(commentsServer, "comments", dbConnection, entityName=reactive(selectedSNP))

} # dispatchOnSNP
#----------------------------------------------------------------------------------------------------
dispatchOnGene <- function(session, selectedGene, origin)
{
  printf("--- dispatchOnGene")

  if(origin != "igv"){
     showGenomicRegion(session, selectedGene)
     displayGeneHancerTrack(session, selectedGene)
     } # !igv

  if(origin != "gene.table"){
    callModule(dataTableServer, "geneDataTable",
               tbl=tbl.genes,
               selectionPolicy="single",
               pageLength=10,
               visibleRows = reactive("all"),
               searchTerm=reactive(selectedGene))
     } # !gene.table

   callModule(messageBoxServer, "messageBox.geneTable", newContent=reactive(selectedGene))
   showGenomicRegion(session, selectedGene)
   callModule(iframeSearchServer, "googleSearch",
              website=reactive("google"), geneSymbol=reactive(selectedGene))
   callModule(iframeSearchServer, "dbSNPSearch",
              website=reactive("dbsnp"), geneSymbol=reactive(selectedGene))
   callModule(iframeSearchServer, "homologeneSearch",
              website=reactive("Homologene"), geneSymbol=reactive(selectedGene))
   callModule(iframeSearchServer, "pubmedSearch",
              website=reactive("PubMed"), geneSymbol=reactive(selectedGene))

   callModule(iframeSearchServer, "clinvarSearch",
              website=reactive("ClinVar"), geneSymbol=reactive(selectedGene))

   callModule(iframeSearchServer, "keggSearch",
              website=reactive("KEGG"), geneSymbol=reactive(selectedGene))
   callModule(iframeSearchServer, "geneCardsSearch",
              website=reactive("GeneCards"), geneSymbol=reactive(selectedGene))
   callModule(commentsServer, "comments", dbConnection, entityName=reactive(selectedGene))



} # dispatchOnSNP
#----------------------------------------------------------------------------------------------------
displayGeneHancerTrack <- function(session, selectedGene)
{
   printf("displayGeneHancerTrack, selectedGene:")
   print(selectedGene)
   print(nchar(selectedGene))

   if(length(selectedGene) > 0){
       if(nchar(selectedGene) > 0){
           printf("entering body of displayGeneHancerTrack: %s", selectedGene)
           tbl.gh <- getEnhancers(tp.hg38, selectedGene)
           if(nrow(tbl.gh) > 0){
              printf("%d enhancers", nrow(tbl.gh))
              #print(tbl.gh)
              tbl.bed <- tbl.gh[, c("chrom", "start", "end", "combinedscore")]
              loadBedGraphTrack(session, trackName="GeneHancer", tbl=tbl.bed, color="blue",
                                autoscale=FALSE, min=1, max=50)
              } # if nrow
           } # if nchar
       } # if length

} # displayGeneHancerTrack
#----------------------------------------------------------------------------------------------------
printf("about to run shiny app on port 9044")
runApp(shinyApp(ui, server), port=3838, launch.browse=FALSE, host="0.0.0.0")
#shinyApp(ui, server)
