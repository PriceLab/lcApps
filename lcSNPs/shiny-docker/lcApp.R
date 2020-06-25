printf <- function(...) print(noquote(sprintf(...)))
library(shiny)
library(shinyModules)
library(later)
library (RSQLite)
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
# geneTab <- function()
# {
#     fluidRow(
#         div(dataTableUI("geneDataTable"),
#             style="margin: 10px; margin-bottom: 30px; padding: 10px; border: 3px solid black; border-radius: 10px;"),
#         div(igvUI("igv2"),
#             style="margin: 10px; margin-bottom: 5px; padding: 10px; border: 3px solid black; border-radius: 10px;"),
#         messageBoxUI(id="messageBox.igv2", title="igvselection", titleSpan=1, boxSpan=8),
#         style="margin-top:5px;"
#         )
#
# } # geneTab
#----------------------------------------------------------------------------------------------------
mainPage <- function()
{

  fluidPage(
    titlePanel("Exploring The Genetic and Molecular Bases of Human Longevity"),
    sidebarPanel(
        radioButtons("genesOrSnps", "Upper Panel", choices=c(Genes="genes", SNPs="snps")),
        verbatimTextOutput("igvSelection"),
        messageBoxUI(id="messageBox.snpTable", title="snp", fontSize=16, boxWidth=120),
        messageBoxUI(id="messageBox.geneTable", title="gene", fontSize=16, boxWidth=120),
        messageBoxUI(id="messageBox.igv", title="igv", fontSize=16, boxWidth=120),
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

   table.selected.snps <- callModule(dataTableServer, "snpDataTable",
                               tbl=tbl.snps,
                               selectionPolicy="single",
                               pageLength=10,
                               visibleRows = reactive("all"))

   table.selected.genes <- callModule(dataTableServer, "geneDataTable",
                               tbl=tbl.genes,
                               selectionPolicy="single",
                               pageLength=10,
                               visibleRows = reactive("all"))

   igv.selected.entity <- callModule(igvServer, "igv",
                                     genome="hg38",
                                     geneModelDisplayMode="COLLAPSED",
                                     locus="chr19:44,854,808-44,940,011") #"APOE")

   observe({
      selectedGenes <- table.selected.genes()
      callModule(messageBoxServer, "messageBox.geneTable", newContent=reactive(selectedGenes))
      showGenomicRegion(session, selectedGenes)

      callModule(iframeSearchServer, "googleSearch",
                 website=reactive("google"), geneSymbol=reactive(selectedGenes))
      callModule(iframeSearchServer, "dbSNPSearch",
                 website=reactive("dbsnp"), geneSymbol=reactive(selectedGenes))
      callModule(iframeSearchServer, "homologeneSearch",
                 website=reactive("Homologene"), geneSymbol=reactive(selectedGenes))
      callModule(iframeSearchServer, "pubmedSearch",
                 website=reactive("PubMed"), geneSymbol=reactive(selectedGenes))
      callModule(iframeSearchServer, "keggSearch",
                 website=reactive("KEGG"), geneSymbol=reactive(selectedGenes))
      callModule(iframeSearchServer, "geneCardsSearch",
                 website=reactive("GeneCards"), geneSymbol=reactive(selectedGenes))
      callModule(commentsServer, "comments", dbConnection, entityName=reactive(selectedGenes))

      })


    observe({
       selectedSNP <- table.selected.snps()
       chrom.loc <- with(tbl.snps[selectedSNP,],  sprintf("%s:%d-%d", chrom, start-10, end+10))
       showGenomicRegion(session, chrom.loc)
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
       })


   observe({
      selectedEntity <- igv.selected.entity()
      callModule(messageBoxServer, "messageBox.igv", newContent=reactive(selectedEntity))
      })

   #callModule(messageBoxServer, "messageBox.geneTable", newContent=table.selected.genes)
   #callModule(messageBoxServer, "messageBox.snpTable", newContent=table.selected.snps)
   #callModule(messageBoxServer, "messageBox.igv", newContent=igv.selected.entity)


#    observe({
#       printf("entering entity observe")
#       geneOrSnp <- selectedEntity();
#       print(geneOrSnp)
#       if(!is.null(geneOrSnp)){
#          printf("--- selected in igv:")
#          print(geneOrSnp)
#          if(nchar(geneOrSnp) > 0){
#             printf("--- selected in igv: %s, calling dbSNP", geneOrSnp)
#             callModule(iframeSearchServer, "dbSNPSearch",
#                        website=reactive("dbSNP"), geneSymbol=reactive(geneOrSnp))
#             printf("--- selected in igv: %s, calling pubmed", geneOrSnp)
#             callModule(iframeSearchServer, "pubmedSearch",
#                        website=reactive("pubmed"), geneSymbol=reactive(geneOrSnp))
#             callModule(iframeSearchServer, "googleSearch",
#                        website=reactive("google"), geneSymbol=reactive(geneOrSnp))
#             } # if nchar > 0
#          } # if !is.null
#       })
#
#
#    observe({
#       rsids <- roi()
#       tbl.sub <- tbl.snps[rsids,]
#       if(nrow(tbl.sub) > 0){
#          chrom.first.encountered <- tbl.sub$chrom[1]
#          tbl.sub <- subset(tbl.sub, chrom==chrom.first.encountered)
#          start.loc <- min(tbl.sub$start)
#          end.loc <- max(tbl.sub$end)
#          shoulder <- round(1 + 0.1 * (end.loc - start.loc))
#          region.string <- sprintf("%s:%d-%d", chrom.first.encountered,
#                                   start.loc - shoulder,
#                                   end.loc + shoulder)
#          showGenomicRegion(session, region.string)
#          } # if nrow
#       })

  later(function(){
      tbl.bed <- tbl.snps[, c("chrom", "start", "end")]
      tbl.bed$name <- rownames(tbl.snps)
      loadBedTrack(session, trackName="snps", tbl=tbl.bed, color="red");
      tbl.bed <- tbl.gh[, c("chrom", "start", "end", "combinedscore", "name")]
      loadBedGraphTrack(session, trackName="GeneHancer", tbl=tbl.bed, color="blue",
                        autoscale=FALSE, min=1, max=50)
      }, 3)

} # server
#----------------------------------------------------------------------------------------------------
printf("about to run shiny app on port 9044")
runApp(shinyApp(ui, server), port=3838, launch.browse=FALSE, host="0.0.0.0")
#shinyApp(ui, server)


