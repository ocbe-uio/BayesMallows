library(dplyr)
library(BayesMallows)

# Load the full data
load("./data-raw/Beaches_data.Rdata")

# We need the data object
beach_preferences <- data %>%
  as_tibble() %>%
  transmute(
    assessor = as.numeric(assessors),
    bottom_item = as.numeric(gsub("[[:alpha:]]", "", firstItem)),
    top_item = as.numeric(gsub("[[:alpha:]]", "", secondItem))
  )

# Now generate the transitive closure
beach_tc <- generate_transitive_closure(beach_preferences)

# Find which assessor have given non-transitive relations
df <- beach_tc %>%
  mutate(
    lowest = pmin(bottom_item, top_item),
    highest = pmax(bottom_item, top_item)
  ) %>%
  group_by(assessor, lowest, highest) %>%
  count()



devtools::use_data(beach_preferences, overwrite = TRUE)
