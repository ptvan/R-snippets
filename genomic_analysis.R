library(maftools)
library(GenVisR)

## load in MAFs file and associated annotations
TCGA_LAML_MAF <- system.file('extdata', 'tcga_laml.maf.gz', package = 'maftools')  
TCGA_LAML_ANNO <- system.file('extdata', 'tcga_laml_annot.tsv', package = 'maftools') 
BRCA_MAF <- system.file('extdata', 'brca.maf.gz', package = 'maftools')  

TCGA_data <- read.maf(maf = TCGA_LAML_MAF, clinicalData = TCGA_LAML_ANNO)
BRCA_data <- read.maf(BRCA_MAF)

## get data summaries
getSampleSummary(TCGA_data)
getGeneSummary(TCGA_data)
getClinicalData(TCGA_data)
getFields(TCGA_data)

##############
## PLOTTING ##
##############

# summary plot
plotmafSummary(TCGA_data, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)

# oncoplots/waterfall plots
oncoplot(TCGA_data, top = 20)

# rainfall plot with change detection
rainfallPlot(maf = BRCA_data, detectChangePoints = TRUE, pointSize = 0.4)

# plotting VAF requires a user-specified column
# this throws an error since no VAF was calculated for the BRCA MAF:
# plotVaf(TCGA_data)
plotVaf(TCGA_data, vafCol = 'i_TumorVAF_WU')

