library(GenomicRanges)
library(GenomicScores)

tbl.snps <- get(load("~/github/lcApps/lcSNPs/shiny-docker/tbl.snps.598.RData"))
dim(tbl.snps)   #  598 6
head(tbl.snps)
           chrom     start       end var     maf subst
rs72903185  chr1  53037233  53037233   A 0.13780   T/A
rs1856036   chr1 101880281 101880281   A 0.45450   G/A
rs12135322  chr1 211549387 211549387   C 0.33510
rs72759967  chr1 228297043 228297043   C 0.08267   T/C
rs61825286  chr1 228301358 228301358   A 0.08107   G/A
rs61825300  chr1 228312214 228312214   C 0.07867   G/C

tbl.genes <- get(load("~/github/lcApps/lcSNPs/shiny-docker/tbl.summary.1010x6-gene-rowNames.RData"))
dim(tbl.genes)  # 1010 5

gr <- GRanges(tbl.snps)
anno <- get(load(system.file(package="TrenaMultiScore", "extdata", "genomeAnnotations.RData")))
class(anno)

gr.annoResults <- annotate_regions(regions=gr, annotations=anno, ignore.strand=TRUE, quiet=FALSE)
tbl.anno <- as.data.frame(gr.annoResults, row.names=NULL)
coi <- c(existing.colnames, "annot.symbol", "annot.type")
coi[1] <- "seqnames"  # bioc-speak for chromosome

snp.associated.genes <- names(table(tbl.anno$annot.symbol))
snp.associated.kengo.genes <- intersect(snp.associated.genes, rownames(tbl.genes))  # 12
