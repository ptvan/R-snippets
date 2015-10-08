# this example workflow reads in .fcs files, creates the associated flowCore-based structures
# performs gating from scratch (ie. we were not provided gated data, usually in the form of flowJo
# workspaces)
# Phu T. Van, Gottardo Lab, FHCRC, September 2015

library(openCyto)
library(data.table)
library(ggplot2)
library(gdata)
library(gridExtra)
library(latticeExtra)

path <- '/home/pvan/flowProject'
dataPath <- '/home/pvan/flowProject/data'
setwd(path)

# this project has each patient's data under a separate dir named for the patientID

fcs_files <- list.files(path, recursive=TRUE, pattern="fcs")
ptids <- list.files(dataPath)

# read in Excel sheet containing metadata 
meta <- read.xls("clinical_data.xlsx")
meta$controller_status <- gsub(" ", "", meta$controller_status)
meta$has_data <- gsub(" ", "", meta$has_data)
meta <- as.data.table(meta)
meta$ptid <- as.character(meta$ptid)
setkey(meta, ptid)

# read all fcs_files into an ncdfFlowset (uses NetCDF library to store flowSets on disk)
# regular flowSets are memory-resident and can be huge
fs  <- read.ncdfFlowSet(fcs_files, ncdfFile = file.path(path, "output", "flowProject.cdf"))

# get markers
markers <- pData(parameters(fs[[1]]))  
markers <- markers[!is.na(markers$desc),]
markers <- markers[!grepl("bead",markers$desc),]
markers <- markers[,c(1:2)]
markers <- data.table(markers)

cytokines <- c("TNFa", 
               "IFNg", 
               "IL-2", 
               "IL-4", 
               "IL-17", 
               "IL-21",
               "Perforin",
               "GrzA",
               "GrzB",
               "CD154",
               "CD107"
)

# make list of channels that need to be transformed
chnls <- as.vector(markers$name)
names(chnls) <- markers$desc

# estimate parameters of the logicle transformation from the data,
# then transform flowSet 
fs_trans <- fs[seq(along=fs)]
fr1 <- fs[["sample1.fcs"]]
tlist <- estimateLogicle(fr1, channels = chnls, type = "data")
fs_trans <- transform(fs_trans, tlist)

# # compare before-vs-after transformation, "after" version should be more spread out 
p0 <- densityplot(~Rh103Di, fs[1], main="CD3 marker, raw data", margin=T)
p1 <- densityplot(~Rh103Di, fs_trans[1], main="logicle transformed", margin=T)
grid.arrange(p0,p1, nrows=2)

# transformation looks good, make an empty (no gates) gatingSet from our transformed data
gs <- GatingSet(fs_trans)

# gatingSets uses phenoData structure from bioC to store metadata
pd <- pData(gs)

# clean metadata and flag samples with their stimulation, so we can facet later
# can also conceivably get this from .fcs file headers using flowCore:::read.FCSheader()
# but users rarely fill this out

pd$ptid <- pd$name
pd$ptid <- substr(pd$ptid, start=0, stop=6)
pd$antigen[grep("ENV|Env|env", pd$name, )] <- "ENV"
pd$antigen[grep("GAG|Gag|gag", pd$name, )] <- "GAG"
pd$antigen[grep("SEB|Seb|seb", pd$name, )] <- "SEB"
pd$antigen[grep("neg|unstim|NA|na", pd$name, )] <- "unstim"

pd <- merge(pd, meta[,.(ptid,controller_status,neut_status)], by="ptid")

# write the metadata to the gatingSet
pData(gs) <- pd

# save the gatingSet, breathe sigh of relief that you have survived this far
save_gs(gs, "output/gs_auto", overwrite=T)

# read gatingSet back in
gs <- load_gs("output/gs_auto/")

# load gatingTemplate 
gtFile <- "cyTOF_gt.csv"
gt <- gatingTemplate(gtFile)

# perform gating using the template. Good things come to those who wait while gating() runs
# can also optionally gate from a particular gate in the hierarchy downstream

gating(gt, gs
       , mc.cores = 4
       , parallel_type = "multicore"
       #, start ="dna"
)

# get some population statistics
getPopStats(gs)
getProp(gs[["first_fcs_file.fcs"]], "CD3")

# example of removing a gate. Afterwards when we run gating() again, first gate to be gated will be CD3
Rm("CD3", gs)

# some basic plotting, details and examples at:
# http://www.bioconductor.org/packages/release/bioc/vignettes/flowWorkspace/inst/doc/plotGate.html

plotGate(gs, "CD3", type="densityplot")
plotGate(gs[["first_fcs_file.fcs"]], "CD3", main="example of one sample's 2D gate")
useOuterStrips(plotGate(gs, "live", type="densityplot",  cond="stim+ptid", main="faceted by stimulation and patientID from pData"))
