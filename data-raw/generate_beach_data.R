library(dplyr)

# Load the full data
load("./data-raw/Beaches_data.Rdata")

# We need the data object
beach_preferences <- data %>%
  as_tibble() %>%
  transmute(
    assessor = as.integer(assessors),
    bottom_item = as.integer(gsub("[[:alpha:]]", "", firstItem)),
    top_item = as.integer(gsub("[[:alpha:]]", "", secondItem))
  )

devtools::use_data(beach_preferences, overwrite = TRUE)
