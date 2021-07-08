# much faster than base R's do.call(rbind, list_of_df)
big_df <- data.table::rbindlist(list_of_dfs)

# globbing files from within R
Sys.glob("/home/pvan/working/data/*/data.txt")
