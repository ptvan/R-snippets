# this is a convenience wrapper to make plotting PCA biplots a bit nicer
# it removes data columns that have zero variance which does not contribute to
# the PCs but interferes with plotting
# it also colors the groups with user-supplied labels

plotpca <- function (df, title=NULL, labels=NULL, varname.adjust=2) {
  require(ggbiplot)  
  require(prcomp)
  # remove data columns that have zero variance
  df <- df[,apply(df,2,var, na.rm=TRUE) != 0]
  # run PCA
  pcobj <- prcomp(df)
  if (length(labels) != nrow(df)){
    cat("length of `labels` did not match rows of input, groups will not be colored in biplot !\n")
    labels <- NULL
  }
  
  if (is.null(labels)){
    # no labels
    out <- ggbiplot(pcobj,...)
  } else {
  out <- ggbiplot(pcobj, groups=labels) +
    ggtitle(title) +
    geom_text(aes(label=labels))
  }
  out  
  return(out)
}