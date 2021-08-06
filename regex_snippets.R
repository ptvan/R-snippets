library(rex)
library(magrittr)

# character classes
valid_chars <- rex(except_some_of(".", "/", " ", "-"))

# combine with built-in rex definitions
re <- rex(start, valid_chars, maybe(":", digit %>% between(2, 5)), end)

re