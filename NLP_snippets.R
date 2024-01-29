library(dplyr) 
library(tidytext)
library(textdata)
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

# topic models
mp_dtm <- DocumentTermMatrix(mp_corpus, 
                   control = list(stemming = TRUE, stopwords = TRUE, minWordLength = 3, removeNumbers = TRUE, removePunctuation = TRUE))



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

# tokenize
sotu_tokens <- tokens(data_corpus_sotu, remove_punct = TRUE)
phrases <- c("united", "union")
sotu_tokens_keep <- tokens_keep(sotu_tokens, phrases)

# apply dictionary
data_dictionary_LSD2015_pos_neg <- data_dictionary_LSD2015[1:2]
sotu_tokens_keep_posneg <- tokens_lookup(sotu_tokens_keep, dictionary = data_dictionary_LSD2015_pos_neg)

# group by speaker
sotu_posneg_dfm <- dfm(sotu_tokens_keep_posneg) %>% 
  dfm_group(groups = dfmat_sotu@docvars$President)

as.data.frame(sotu_posneg_dfm)

## modelling topics on AssociatedPress data
data("AssociatedPress")
ap_lda <- LDA(AssociatedPress, 
              k = 3, 
              control = list(seed = 31415))

# extract beta, per-topic-per-word probabilities
ap_topics <- tidy(ap_lda, matrix = "beta")

# extract top terms, group them by topics and plot
ap_top_terms <- ap_topics %>%
  group_by(topic) %>%
  slice_max(beta, n = 10) %>% 
  ungroup() %>%
  arrange(topic, -beta)

ap_top_terms %>%
  mutate(term = reorder_within(term, beta, topic)) %>%
  ggplot(aes(beta, term, fill = factor(topic))) +
  geom_col(show.legend = FALSE) +
  facet_wrap(~ topic, scales = "free") +
  scale_y_reordered()
