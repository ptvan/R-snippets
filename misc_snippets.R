# toggle scientific notation off for session
options(scipen=999)

# disable scientific notation for a single value
format("1000000000", scientific=F)

# bind a list of same-rank dfs into one 
# much faster than base R's do.call(rbind, list_of_df)
big_df <- data.table::rbindlist(list_of_dfs)

# globbing files from within R
Sys.glob("/home/pvan/working/data/*/data.txt")

# initialize an empty df with same columns as another df
new_df <- df[FALSE, ]

# remap vector based on another vector
old <- c("x", "y", "x", "z")
mapping <- c("x" = "a", "y" = "b", "z" = "c")
new <- mapping[old]

# categorize numeric values into discrete factors
cut(1:100, breaks=3, labels = c("low", "medium", "high"))

# convert a vector of strings to CamelCase
re_from <- "\\b([[:lower:]])([[:lower:]]+)"
strings <- c("first phrase", "another phrase to convert",
             "and here's another one", "last-one")
gsub(re_from, "\\U\\1\\L\\2" ,strings, perl=TRUE)
