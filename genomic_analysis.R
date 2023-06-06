library(maftools)
library(GenVisR)
library(MutationalPatterns)
library(BSgenome)
library(NMF)

## load in MAFs file and associated annotations
TCGA_LAML_MAF <- system.file('extdata', 'tcga_laml.maf.gz', package = 'maftools')  
TCGA_LAML_ANNO <- system.file('extdata', 'tcga_laml_annot.tsv', package = 'maftools') 
BRCA_MAF <- system.file('extdata', 'brca.maf.gz', package = 'maftools')  

TCGA_data <- read.maf(maf = TCGA_LAML_MAF, clinicalData = TCGA_LAML_ANNO)
BRCA_data <- read.maf(BRCA_MAF)

## load in VCF files as GRanges via MutationalPatterns
ref_genome <- "BSgenome.Hsapiens.UCSC.hg19"
library(ref_genome, character.only = TRUE)

vcf_files <- list.files(system.file("extdata", package = "MutationalPatterns"),
                        pattern = "sample.vcf", full.names = TRUE
)

## tissue type as a vector for plotting
vcf_metadata <- c(
  "colon1", "colon2", "colon3",
  "intestine1", "intestine2", "intestine3",
  "liver1", "liver2", "liver3"
)

grl <- read_vcfs_as_granges(vcf_files, vcf_metadata, ref_genome)

## get data summaries
getSampleSummary(TCGA_data)
getGeneSummary(TCGA_data)
getClinicalData(TCGA_data)
getFields(TCGA_data)

#####################
## DATA SUBSETTING ##
#####################

TCGA_single_gene <- subsetMaf(TCGA_data, genes = "FLT3")
TCGA_single_sample_df <- subsetMaf(TCGA_data, 
                                tsb = "TCGA-AB-3009",
                                query = "Variant_Classification == 'Missense_Mutation'",
                                mafObj = FALSE
                                )
TCGA_survivor_samples_df <- subsetMaf(TCGA_data, 
                                   clinQuery = "Overall_Survival_Status == '1'",
                                   mafObj = FALSE
)

## Mutation Spectra

# extract subtypes
subtype_counts <- mut_type_occurrences(grl, ref_genome) 

# calculate trinucleotide spectra
mutation_matrix <- mut_matrix(vcf_list = grl, ref_genome = ref_genome)

# can also calculate larger contexts
mutation_matrix_ext_context <- mut_matrix(vcf_list = grl, 
                                          ref_genome = ref_genome, 
                                          extension = 2)

##############
## PLOTTING ##
##############

## multi-panel summary plot
plotmafSummary(TCGA_data, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)

## simple barplot of mutation types
mafbarplot(TCGA_data, includeCN = TRUE)

## oncoplots/waterfall plots
oncoplot(TCGA_data, top = 20)

## rainfall plot with change detection
rainfallPlot(maf = BRCA_data, detectChangePoints = TRUE, pointSize = 0.4)

## plotting VAF requires a user-specified column
# this throws an error since no VAF was calculated for the BRCA MAF:
# plotVaf(TCGA_data)
plotVaf(TCGA_data, vafCol = 'i_TumorVAF_WU')

# facet spectrum plot by vector of metadata
plot_spectrum(subtype_counts, 
              by = vcf_metadata, 
              indv_points = TRUE,
              error_bars = 'none')

# plot trinucleotide spectrum, comparing across samples
plot_96_profile(mutation_matrix[, c(1, 7)])

# plot extended context spectrum as a heatmap, also comparing across samples
plot_profile_heatmap(mutation_matrix_ext_context, by = vcf_metadata)
