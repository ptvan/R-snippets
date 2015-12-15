library(COMPASS)
library(data.table)

# extract cell counts from the COMPASS container (CR[[antigen]]$data$count_s)
#

COMPASS_counts_to_SPICE <- function(CR){
  
  antigens <- names(CR)
  
  # get cytokines (last column is "Counts", exclude it)
  cytokines <- colnames(CR[[1]]$data$categories)
  cytokines <- cytokines[-length(cytokines)]
  
  # create the big data.frame to be returned later
  big <- data.frame(matrix(vector(),0, length(cytokines)+4))
  colnames(big) <- c("ptid","time",cytokines,"proportion","antigen")
  
  # for each antigen in the COMPASSResult structure...  
  for (i in 1:length(antigens)){
    antigen <- antigens[i]
    cat(antigen, "\n")  
    p <- as.data.frame(CR[[antigen]]$data$n_s)
    # get total cell counts
    totals <- CR[[antigen]]$data$counts_s
    
    p_new <- data.frame(matrix(vector(),0, length(cytokines)+4))
    colnames(p_new) <- c("ptid","time",cytokines,"proportion","antigen")
    
    k <- 0
    for (i in 1:nrow(p)){
      rown <- rownames(p)[i]
      ptid <- unlist(strsplit(rown, "_"))[1]
      time <- unlist(strsplit(rown, "_"))[2]
      
      # calculate proportion
      for (j in 1:ncol(p)){
        coln <- colnames(p)[j]
        x <- unlist(strsplit(coln, "&"))
        vec <- rep("+", 5)
        vec[grep("\\!", x)] = "-"
        vec <- c(ptid, time, vec, p[i,j] / totals[i], antigen)
        k <- k+1
        # cat("row", k, ":  ", vec, "\n")
        p_new[k,] <- vec 
      }
    }
    big <- rbind(big,p_new)
  }
  big$time <- as.factor(big$time)
  big$proportion <- as.numeric(big$proportion)
  big <- as.data.table(big)
  return(big)
}