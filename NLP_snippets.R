library(dplyr) 
library(tidytext) 
library(topicmodels) 
library(tm) 
library(SnowballC)
library(janeaustenr)
library(quanteda)
library(quanteda.textstats)
library(quanteda.corpora)

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
topfeatures(mp_dfm)

## create Feature Co-occurrence Matrix (FCM)
mp_fcm <- fcm(mp_dfm)
topfeatures(mp_fcm)

## descriptive stats

# word frequency
mp_freq <- textstat_frequency(mp_dfm, n = 5)

# lexical diversity
mp_lexdiv <- textstat_lexdiv(mp_dfm)

## use US Presidents State of the Union address corpus
# tokenize
dfmat_sotu <- data_corpus_sotu %>% 
  tokens(remove_punct = TRUE, remove_url = TRUE, remove_symbols = TRUE) %>% 
  dfm() 

# create DFM
ndoc(dfmat_sotu)
dfm_group(dfmat_sotu, groups = "docs")
dfm_speakers <- dfm_group(dfmat_sotu, groups =  dfmat_sotu@docvars$President)

# extract and filter speakers
dfm_speakers <- dfm_speakers %>% 
  dfm_select(min_nchar = 2) %>% 
  dfm_trim(min_termfreq = 10) 

# cluster speakers
tstat_dist <- as.dist(textstat_dist(dfm_speakers))
speaker_clust <- hclust(tstat_dist)
plot(speaker_clust)
