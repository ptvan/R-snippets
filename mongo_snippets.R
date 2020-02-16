library(mongolite)
m <- mongo(collection = "test", db = "test", url = "mongodb://localhost",
           verbose = FALSE, options = ssl_options())