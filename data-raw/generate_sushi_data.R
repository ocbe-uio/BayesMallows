# Generate sushi data

# Load it
load("./data-raw/sushiData.RData")

# Rename and save
sushi_rankings <- Ro
colnames(sushi_rankings) <- c(
  "shrimp", "sea eel", "tuna", "squid", "sea urchin",
  "salmon roe", "egg", "fatty tuna", "tuna roll", "cucumber roll"
)

devtools::use_data(sushi_rankings, overwrite = TRUE)
