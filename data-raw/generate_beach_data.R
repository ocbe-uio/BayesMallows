library(BayesMallows)

# Load the full data
load("./data-raw/Beaches_data.Rdata")

# We need the data object
data <- as.data.frame(data)
beach_preferences <- data.frame(
  assessor = as.numeric(data$assessors),
  bottom_item = as.numeric(gsub("[[:alpha:]]", "", data$firstItem)),
  top_item = as.numeric(gsub("[[:alpha:]]", "", data$secondItem))
)

usethis::use_data(beach_preferences, overwrite = TRUE)
