# this example workflow reads in .fcs files, creates the associated flowCore-based structures
# performs gating from scratch (ie. we were not provided gated data, usually in the form of flowJo
# workspaces)
# Phu T. Van, Gottardo Lab, FHCRC, September 2015

library(openCyto)
library(data.table)
library(ggplot2)
library(gdata)

path <- '/home/pvan'
dataPath <- '/home/pvan/data/flowProject/'
setwd(path)

# this project has each patient's data under a separate dir named for the patientID

fcs_files <- list.files(path, recursive=TRUE, pattern="fcs")
ptids <- list.files(dataPath)

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

# # compare before-vs-after transformation, after should 
p0 <- densityplot(~Rh103Di, fs[1], main="CD3 marker, raw data", margin=T)
p1 <- densityplot(~Rh103Di, fs_trans[1], main="logicle transformed", margin=T)

# transformation looks good, make an empty (ungated) gatingSet from our transform data
gs <- GatingSet(fs_trans)

# gatingSets ues phenoData structure from bioC to store metadata
pd <- pData(gs)

# clean things up and flag samples with their stimulation
# can also conceivably get this from .fcs file headers using flowCore:::read.FCSheader()
# but users rarely fill this out

pd$ptid <- pd$name
pd$ptid <- substr(pd$ptid, start=0, stop=6)
pd$antigen[grep("ENV|Env|env", pd$name, )] <- "ENV"
pd$antigen[grep("GAG|Gag|gag", pd$name, )] <- "GAG"
pd$antigen[grep("SEB|Seb|seb", pd$name, )] <- "SEB"
pd$antigen[grep("neg|unstim|NA|na", pd$name, )] <- "unstim"

pd <- merge(pd, meta[,.(ptid,controller_status,neut_status)], by="ptid")

pData(gs) <- pd

# save the gatingSet, breathe sigh of relief that you have survived this far
save_gs(gs, "output/gs_auto", overwrite=T)

# read gatingSet back in
gs <- load_gs("output/gs_auto/")

# load gatingTemplate 
gtFile <- "cyTOF_gt.csv"
gt <- gatingTemplate(gtFile)

# perform gating using the template. Good things come to those who wait while gating() runs
gating(gt, gs
       , mc.cores = 4
       , parallel_type = "multicore"
       #, start ="dna"
)

# optionally remove a gate
Rm("CD3", gs)

# consult plotGate's extensive help for further plotting shenanigans
plotGate(gs, "CD3", type="densityplot")
