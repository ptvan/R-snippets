library(dplyr) 
library(tidytext) 
library(topicmodels) 
library(tm) 
library(SnowballC)
library(janeaustenr)
library(quanteda)

## generate a corpus of Jane Austen's Mansfield Park which is a character vector
# if called on data.frame, needs `text_field` param 
mp_corpus <- corpus(mansfieldpark)

## tokenize and view tokens in context
mp_tokens <- tokens(mp_corpus, remove_punct = TRUE)
mp_love <- kwic(mp_tokens, pattern = "love")
mp_your_love <- kwic(mp_tokens, pattern = phrase(c("my love", "our love")))

## generate N-grams
mp_ngrams <-  tokens_ngrams(mp_tokens, n = 2:4)

## create Document Feature Matrix (DFM)
mp_dfm <- dfm(mp_tokens)
