old_wd <- getwd()

setwd("vignettes-raw/")
knitr::knit("BayesMallows.Rmd", output = "../vignettes/BayesMallows.Rmd")
knitr::knit("SMC-Mallows.Rmd", output = "../vignettes/SMC-Mallows.Rmd")
knitr::knit("parallel_chains.Rmd", output = "../vignettes/parallel_chains.Rmd")

imgs <- list.files(pattern = "\\.png$")
imgs_new <- file.path("..", "vignettes", imgs)
file.rename(imgs, imgs_new)
setwd(old_wd)
beepr::beep()
