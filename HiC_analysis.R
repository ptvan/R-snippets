library(GenomicRanges)
library(HiCExperiment)
library(HiContacts)
library(HiContactsData)
library(InteractionSet)

# load an .mcool file cached by ExperimentHub
coolf <- HiContactsData('yeast_wt', 'mcool')


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

# visualization
## general heatmap, using some plot options
plotMatrix(first3chr, use.scores = 'count', limits = c(-3.5, -1))

## horizontal map
plotMatrix(hic, maxDistance = 100000, use.scores = 'balanced')
