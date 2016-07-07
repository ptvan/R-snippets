library(COMPASS)
library(data.table)


COMPASS_cell_proportions <- function(CR){
# extract cell counts from the COMPASSResult container (CR[[antigen]]$data$count_s)
# !!! IMPORTANT: assumes data comes from time-course experiment and therefore has
# joint_id (<ptid>_<time>)
  
# convert it to proportions and export to a format (mostly) readable by SPICE
# null category is kept, so proportions should be very small since cytokines are rare

  
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
    totals_u <- CR[[antigen]]$data$counts_u
    
    p_new <- data.frame(matrix(vector(),0, length(cytokines)+5))
    colnames(p_new) <- c("ptid","time",cytokines,"unstim_proportion", "stim_proportion","antigen")
    
    k <- 0
    for (i in 1:nrow(p)){
      rown <- rownames(p)[i]
      # split jointID into ptid and time
      # eg. "1_0", ptid == 1, time == 0
      
      ptid <- unlist(strsplit(rown, "_"))[1]
      time <- unlist(strsplit(rown, "_"))[2]
      
      for (j in 1:ncol(p)){ # NOTE: this assumes count_s & count_u have same dimensions
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


COMPASS_cell_proportions_nonull <- function(CR){
## !!! IMPORTANT !!!
## This function calculates percentages EXCLUDING THE NULL CATEGORY !!!
## eg. you have 10000 cells 
## 5 cells express IL2 only, 3 cells express IFNg only, 2 cells express both (10 cells express any)
## the percentages returned by this fx are IL2 = 5/10 = 0.5, IFNg = 3/10 = 0.3, IL2+IFNg = 2/10 = 0.2
## this is different than when the calculation is over ALL cells (10000 in this example)
  
  if(!is.list(CR) || class(CR[[1]]) != "COMPASSResult") {
    stop("this function requires a list containing at least 1 COMPASSResult !!!")
  }
  
  antigens <- names(CR)
  
  # create the big data.frame to be returned later
  big <- data.frame(matrix(vector(),0, length(cytokines)+6))
  colnames(big) <- c("ptid","time",cytokines,"unstim_proportion" ,"stim_proportion","antigen")
  
  # get cytokines (last column is "Counts", exclude it)
  cytokines <- colnames(CR[[1]]$data$categories)
  cytokines <- cytokines[-length(cytokines)]  
  
  for (i in 1:length(antigens)){
    antigen <- antigens[i]
    cat(antigen, "\n")  
    
    # get cell counts, REMOVE NULL CATEGORY
    p <- as.data.frame(CR[[antigen]]$data$n_s)
    p <- p[,1:ncol(p)-1]
    p_u <- as.data.frame(CR[[antigen]]$data$n_u)
    p_u <- p_u[,1:ncol(p_u)-1]
    
    # RECALCULATE TOTAL
    totals <- apply(p, 1, sum)
    totals_u <- apply(p_u, 1, sum)
    
    p_new <- data.frame(matrix(vector(),0, length(cytokines)+5))
    colnames(p_new) <- c("ptid","time",cytokines,"unstim_proportion", "stim_proportion","antigen")
    k <- 0
    for (i in 1:nrow(p)){
      rown <- rownames(p)[i]
      # since this is time-series, need to split jointID into ptid and time
      # eg. "1_0", ptid == 1, time == 0
      
      ptid <- unlist(strsplit(rown, "_"))[1]
      time <- unlist(strsplit(rown, "_"))[2]
      
      for (j in 1:ncol(p)){ # NOTE: this assumes count_s & count_u have same dimensions
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
 
