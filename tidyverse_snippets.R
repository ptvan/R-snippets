library(tidyverse)

# create new columns dynamically using bang-bang and walrus operators
new_col_name <- "total_score"
my_tibble <- tibble(score1 = 10, score2 = 20)
my_tibble %>% 
  mutate(!!new_col_name := score1 + score2)

# same process but providing a vector of column names
purrr::map_dfc(COLUMN_NAMES, ~ tibble(!!.x := rep(0, nrows(my_df)))) 

