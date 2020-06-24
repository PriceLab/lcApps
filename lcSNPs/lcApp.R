library(shiny)
library(later)
source("~/github/shinyModules/igv/igvModule.R")
source("~/github/shinyModules/dataTable/dataTableModule.R")
source("~/github/shinyModules/messageBox/messageBoxModule.R")
source("~/github/shinyModules/iframeSearch/iframeSearchModule.R")

load("../genehancer/tbl.gh.apoe.RData")
colnames(tbl.gh)[grep("ghid", colnames(tbl.gh))] <- "name"

tbl.snps <- read.table("threeSnps.tsv", header=TRUE, sep="\t", as.is=TRUE, row.names=1)
#----------------------------------------------------------------------------------------------------
snpsTab <- function()
{
   #wellPanel(
       fluidRow(
           div(igvUI("igv"),
               style="margin: 10px; margin-bottom: 5px; padding: 10px; border: 3px solid black; border-radius: 10px;"),
           messageBoxUI(id="messageBox.igv", title="igvselection", titleSpan=1, boxSpan=8),
           div(dataTableUI("snpDataTable"),
               style="margin: 10px; margin-bottom: 30px; padding: 10px; border: 3px solid black; border-radius: 10px;"),
           #actionButton("addSNPsTrackButton", "Add SNP Track"),
           #actionButton("addGHTrackButton", "Add GeneHancer Track"),
           #messageBoxUI(id="messageBox.igv", title="igvselection", titleSpan=2, boxSpan=8),
           #DTOutput("geneTable"),
           style="margin-top:5px;")
    #   )

} # snpsTab
#----------------------------------------------------------------------------------------------------
ui <- fluidPage(
   tabsetPanel(type="tabs", id="lcGenesTabs",
               #introTab(),
               tabPanel(title="SNP-centric View", snpsTab()),
               tabPanel(title="dbSNP", wellPanel(iframeSearchUI(id="dbSNPSearch",   title="dbSNP"))),
               tabPanel(title="PubMed", wellPanel(iframeSearchUI(id="pubmedSearch", title="PubMed"))),
               tabPanel(title="Google", wellPanel(iframeSearchUI(id="googleSearch", title="Google")))
               ),
    style="margin: 10px; margin-top: 10px; margin-bottom: 50px;"
   )


#----------------------------------------------------------------------------------------------------
server <- function(input, output, session){

   roi <- reactiveVal("APOE")

   roi <- callModule(dataTableServer, "snpDataTable",
                     tbl=tbl.snps,
                     selectionPolicy="multiple",
                     pageLength=10,
                     visibleRows = reactive("all"))

   selectedEntity <- callModule(igvServer, "igv",
                                genome="hg38",
                                geneModelDisplayMode="COLLAPSED",
                                locus="chr19:44,854,808-44,940,011") #"APOE")

   observe({
      printf("entering entity observe")
      geneOrSnp <- selectedEntity();
      print(geneOrSnp)
      if(!is.null(geneOrSnp)){
         printf("--- selected in igv:")
         print(geneOrSnp)
         if(nchar(geneOrSnp) > 0){
            printf("--- selected in igv: %s, calling dbSNP", geneOrSnp)
            callModule(iframeSearchServer, "dbSNPSearch",
                       website=reactive("dbSNP"), geneSymbol=reactive(geneOrSnp))
            printf("--- selected in igv: %s, calling pubmed", geneOrSnp)
            callModule(iframeSearchServer, "pubmedSearch",
                       website=reactive("pubmed"), geneSymbol=reactive(geneOrSnp))
            callModule(iframeSearchServer, "googleSearch",
                       website=reactive("google"), geneSymbol=reactive(geneOrSnp))
            } # if nchar > 0
         } # if !is.null
      })

   callModule(messageBoxServer, "messageBox.igv", newContent=selectedEntity)

   observe({
      rsids <- roi()

      tbl.sub <- tbl.snps[rsids,]
      if(nrow(tbl.sub) > 0){
         chrom.first.encountered <- tbl.sub$chrom[1]
         tbl.sub <- subset(tbl.sub, chrom==chrom.first.encountered)
         start.loc <- min(tbl.sub$start)
         end.loc <- max(tbl.sub$end)
         shoulder <- round(1 + 0.1 * (end.loc - start.loc))
         region.string <- sprintf("%s:%d-%d", chrom.first.encountered,
                                  start.loc - shoulder,
                                  end.loc + shoulder)
         showGenomicRegion(session, region.string)
         } # if nrow
      })

  later(function(){
      tbl.bed <- tbl.snps[, c("chrom", "start", "end")]
      tbl.bed$name <- rownames(tbl.snps)
      loadBedTrack(session, trackName="snps", tbl=tbl.bed, color="red");
      tbl.bed <- tbl.gh[, c("chrom", "start", "end", "combinedscore", "name")]
      loadBedGraphTrack(session, trackName="GeneHancer", tbl=tbl.bed, color="blue",
                        autoscale=FALSE, min=1, max=50)
      }, 5)

} # server
#----------------------------------------------------------------------------------------------------
runApp(shinyApp(ui, server), port=9044)

