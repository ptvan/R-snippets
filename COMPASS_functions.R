library(COMPASS)
library(data.table)

# extract cell counts from the COMPASSResult container (CR[[antigen]]$data$count_s)
# convert it to proportions and export to a format (mostly) readable by SPICE

COMPASS_counts_to_SPICE <- function(CR){
  
  if(!is.list(CR) || class(CR[[1]]) != "COMPASSResult") {
    stop("this function requires a list containing at least 1 COMPASSResult !!!")
  }
  
  antigens <- names(CR)
  
  # get cytokines (last column is "Counts", exclude it)
  cytokines <- colnames(CR[[1]]$data$categories)
  cytokines <- cytokines[-length(cytokines)]
  
  # create the big data.frame to be returned later
  big <- data.frame(matrix(vector(),0, length(cytokines)+6))
  colnames(big) <- c("ptid","time",cytokines,"unstim_proportion" ,"stim_proportion","antigen")
  
  # for each antigen in the COMPASSResult structure...  
  for (i in 1:length(antigens)){
    antigen <- antigens[i]
    cat(antigen, "\n")  
    p <- as.data.frame(CR[[antigen]]$data$n_s)
    p_u <- as.data.frame(CR[[antigen]]$data$n_u)
    
    # get total cell counts
    totals <- CR[[antigen]]$data$counts_s
    totals_u <- CR[[antigen]]$data$counts_s
    
    p_new <- data.frame(matrix(vector(),0, length(cytokines)+5))
    colnames(p_new) <- c("ptid","time",cytokines,"unstim_proportion", "stim_proportion","antigen")
    
    k <- 0
    for (i in 1:nrow(p)){
      rown <- rownames(p)[i]
      # since this is time-series, need to split jointID into ptid and time
      # eg. "1_0", ptid == 1, time == 0
      
      ptid <- unlist(strsplit(rown, "_"))[1]
      time <- unlist(strsplit(rown, "_"))[2]
      
      for (j in 1:ncol(p)){ # this assumes count_s & count_u have same dimensions
        coln <- colnames(p)[j]
        x <- unlist(strsplit(coln, "&"))
        vec <- rep("+", 5)
        vec[grep("\\!", x)] = "-"
        vec <- c(ptid, time, vec, p[i,j] / totals[i], p_u[i,j] / totals_u[i], antigen)
        k <- k+1
        # cat("row", k, ":  ", vec, "\n")
        p_new[k,] <- vec 
      }
    }
    big <- rbind(big,p_new)
  }
  # append
  big$time <- as.factor(big$time)
  big$unstim_proportion <- as.numeric(big$unstim_proportion)
  big$stim_proportion <- as.numeric(big$stim_proportion)
  big <- as.data.table(big)
  return(big)
}