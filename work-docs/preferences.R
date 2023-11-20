devtools::load_all()
dd <- beach_preferences
dd$assessor <- factor(dd$assessor)
levels(dd$assessor) <- sample(levels(dd$assessor))
data <- setup_rank_data(
  preferences = dd
)

