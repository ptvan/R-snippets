library(maftools)
library(GenVisR)

# load in MAF file and associated annotations
FILE_PATH <- system.file('extdata', 'tcga_laml.maf.gz', package = 'maftools')  
ANNO_PATH <- system.file('extdata', 'tcga_laml_annot.tsv', package = 'maftools') 

maf_data <- read.maf(maf = FILE_PATH, clinicalData = ANNO_PATH)

# get data summaries
getSampleSummary(maf_data)
getGeneSummary(maf_data)
getClinicalData(maf_data)
getFields(maf_data)

# plotting
plotmafSummary(maf = maf_data, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)
oncoplot(maf_data, top = 20)
