# this example workflow reads in .fcs files, creates the associated flowCore-based structures
# performs gating from scratch (ie. we were not provided gated data, normally flowJo workspaces)
# the data used is cyTOF, so no compensation is needed, only transformation (using logicle transforms)

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
# this example has clinical variables "controller_status", among others

meta <- read.xls("clinical_data.xlsx")
meta$controller_status <- gsub(" ", "", meta$controller_status)
meta$has_data <- gsub(" ", "", meta$has_data)
meta <- as.data.table(meta)
meta$ptid <- as.character(meta$ptid)
setkey(meta, ptid)

# read all fcs_files into an ncdfFlowset (uses NetCDF library to store flowSets on disk)
# regular flowSets are memory-resident and can be huge
# DO NOT MOVE OR DELETE THIS .cdf FILE OR THINGS WILL BREAK
fs  <- read.ncdfFlowSet(fcs_files, ncdfFile = file.path(path, "output", "flowProject.cdf"))

# get marker/channel information from the flowSet
markers <- pData(parameters(fs[[1]]))  

# make list of channels that need to be transformed
chnls <- as.vector(markers$name)
names(chnls) <- markers$desc

# for non-cyTOF data you would also need to compensate here...
# flowCore provides compensate() function, ?flowCore::compensate

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

# example of cleaning metadata and flagging samples with their stimulation, so we can facet later
# can also conceivably get this from .fcs file headers using flowCore:::read.FCSheader()
# but users rarely fill this out

pd$ptid <- pd$name
pd$ptid <- substr(pd$ptid, start=0, stop=6)
pd$antigen[grep("ENV|Env|env", pd$name )] <- "ENV"
pd$antigen[grep("GAG|Gag|gag", pd$name )] <- "GAG"
pd$antigen[grep("SEB|Seb|seb", pd$name )] <- "SEB"
pd$antigen[grep("neg|unstim|NA|na", pd$name )] <- "unstim"

pd <- merge(pd, meta[,.(ptid,controller_status,neut_status)], by="ptid")

# write the metadata to the gatingSet
pData(gs) <- pd

# save the gatingSet, breathe sigh of relief that you have survived this far
# note that at this point you have created a gatingSet, but since 
# you have not gated your data, this gatingSet is empty

save_gs(gs, "output/gs_auto")

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

# plots the gating hierarchy, which now shouldn't be empty
plot(gs)

# get some population statistics
getPopStats(gs)
getProp(gs[["first_fcs_file.fcs"]], "CD3")

# gatingSets only contain the gates, the flow event data are stored in flowSets
# but thankfully, flowSets are attached to gatingSets, so you can get the event data:
# note that this data is in the state you used to create the gatingSet (in our case, this is transformed data)

# this gets you a flowFrame of flow events for the first .fcs file, *DOWNSTREAM* of the `CD4+` gate
fr <- getData(gs[["first_fcs_file.fcs"]], "CD4+")

# openCyto has various functions for manipulating flowFrames, so safer than others
fr <- openCyto:::.truncate_flowframe(fr, channels = "Rh103Di", min = 1)

# this plots a histogram of the data using base R
x <- as.vector(exprs(fr)[, "Rh103Di"])
hist(x)

# you can set various parameters used by flowWorkspace
flowWorkspace.par.set("plotGate", list("default.y" = "Rh103Di"))
flowWorkspace.par.set("plotGate", list("type" = "histogram"))

# example of removing a gate. Afterwards when we run gating() again, first gate to be gated will be CD3
# removing a gate will also remove all downstream gates
Rm("CD4+", gs)

# some basic plotting, details and examples at:
# http://www.bioconductor.org/packages/release/bioc/vignettes/flowWorkspace/inst/doc/plotGate.html
plotGate(gs, "CD3", type="densityplot")
plotGate(gs[["first_fcs_file.fcs"]], "CD3", main="example of one sample's 2D gate")
useOuterStrips(plotGate(gs, "live", type="densityplot",  cond="stim+ptid", main="density plot faceted by stimulation and patientID from pData"))

# you can also exclude samples from the gatingSet by subsetting
# this removes all unstimulated (control) samples, as determined by the gatingSet's phenoData
# one more reason to fill out your pData(gs) correctly !!!
gs <- subset(gs, !stim %in% c("unstim") )

# changes to the gatingSet are only in memory, so you need to explicitly save.
# NOTE: you will not a warning when overwrite=TRUE, so be careful !!!
save_gs(gs, "output/gs_auto", overwrite=TRUE)
