require(tidyverse)

# extract mutational signature info from 
# https://cancer.sanger.ac.uk/cosmic/signatures
# Correct as of date: 2019-01-27
sig_info <- function(data, col, sig_i) {

  # each field has a specific h4 index
  description_col <- c(cancer_type = 3,
                       proposed_aetiology = 4,
                       additional_features = 5,
                       comments = 6)
  # don't ask - probably better ways but this works fine
  xml_child(xml_child(xml_child(data, sig_i), description_col[col]), 2) %>% 
    rvest::html_text() %>% 
    trimws()
}

# we're interested in the second div of class panel
x <-
  xml2::read_html("https://cancer.sanger.ac.uk/cosmic/signatures") %>% 
  rvest::html_nodes(".panel") %>% 
  `[[`(2)

sig_table <- lapply(1:30, function(i) {
  tibble::tibble(
    cancer_type = sig_info(x, "cancer_type", i),
    proposed_aetiology = sig_info(x, "proposed_aetiology", i),
    additional_features = sig_info(x, "additional_features", i),
    comments = sig_info(x, "comments", i)
  )
})

sig_table %>% 
  dplyr::bind_rows() %>%
  readr::write_tsv(path = "sig_description.tsv")

