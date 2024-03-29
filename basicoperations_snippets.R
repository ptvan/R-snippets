### Matrix operations

# trace of a covariance matrix
V <- matrix(c(10,-5,10,-5,20,10,0,30), nrow=3))
sum(diag(V))

W <- matrix(c(5,-10,10,-5,10,20,30,0), nrow=3))

# matrix multiplication
W %*% V

# outer product
W %o% V

# W'B
crossprod(W,V)

# Eigen vectors & values
# eigen(V)$vectors[,i] is weights of principal component i
# eigen(V)$values[i]/trace is %variation explained by principal component i
eigen(V)

# toggle scientific notation off for session
options(scipen=999)

# disable scientific notation for a single value
format("1000000000", scientific=F)

# rank() assigns ranks to each element in a vector
scores <- c(80, 60, 90, 75)
ranking <- rank(scores) # returns "3 1 4 2"

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
