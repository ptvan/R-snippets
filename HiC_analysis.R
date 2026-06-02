library(GenomicRanges)
library(HiCExperiment)
library(HiContacts)
library(HiContactsData)
library(InteractionSet)
library(rtracklayer)
library(BSgenome.Hsapiens.UCSC.hg38)
library(BiocParallel)

# load an .mcool file cached by ExperimentHub
coolf <- HiContactsData('yeast_wt', 'mcool')

# load a BEDPE containing chromatin loops as an GInteraction object
loops <- system.file('extdata', 'S288C-loops.bedpe', package = 'HiCExperiment') |> 
  import() |> 
  makeGInteractionsFromGRangesPairs()

#############################
# EXPLORATORY DATA ANALYSIS
#############################

# this particular .mcool has 5 resolutions
cf <- CoolFile(coolf)

# import at 2000 resolution
hic <- import(cf, focus = 'II:300001-813184', resolution = 2000)

scores(hic)
scores(hic, "count")
scores(hic, "balanced")

# subsetting
## using subsetByOverlaps() ...
telomere <- GRanges("II:700001-813184")
subsetByOverlaps(hic, telomere) |> interactions()

## ... or directly on the HiCExperiment object
hic["II:800001-813184"]
first3chr <- hic[c("I","II","III")]

# changing resolutions while preserving metadata (aka. "zooming")
# NOTE: only works on .hic and .mcool, not .cool

length(hic)
length(zoom(hic, 1000))
length(zoom(hic, 4000))

# visualizing contact maps
## general heatmap, showing different plot options
plotMatrix(first3chr, use.scores = 'count', limits = c(-3.5, -1))

## horizontal map
plotMatrix(hic, maxDistance = 100000, use.scores = 'balanced')

## layer on the loops that we loaded previously
plotMatrix(hic, loops = loops)

# Aggregated Plot Analysis (APA)
aggr_loops <- aggregate(hic, targets = loops, flankingBins = 15)

slices(aggr_loops)

topologicalFeatures(aggr_loops, 'targets')

plotMatrix(
  aggr_loops, 
  use.scores = 'detrended', 
  scale = 'linear', 
  limits = c(-1, 1), 
  cmap = bgrColors()
)

# Finding Topological Features
## annotate compartments
human_microC <- import(CoolFile(HiContactsData('microC', 'mcool')), resolution = 250000)
hg38 <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
human_compts <- getCompartments(human_microC, genome = hg38)

topologicalFeatures(human_compts, "compartments")

metadata(human_compts)$eigens

## export eigenvector as a bigwig file 
coverage(metadata(human_compts)$eigens, weight = 'eigen') |> export('microC_eigen.bw')

## export compartments as GFF file
topologicalFeatures(human_compts, "compartments") |> export('microC_compartments.gff3')

## `observed vs. expected` interaction scores, aka "saddle plots"
plotSaddle(human_compts, nbins = 25, BPPARAM = SerialParam(progressbar = FALSE))
