---
output:
  pdf_document: default
  html_document: default
---
``` r
library(BayesMallows)
sessionInfo()
#> R version 4.3.2 (2023-10-31)
#> Platform: aarch64-apple-darwin20 (64-bit)
#> Running under: macOS Sonoma 14.1.2
#>
#> Matrix products: default
#> BLAS:   /Library/Frameworks/R.framework/Versions/4.3-arm64/Resources/lib/libRblas.0.dylib
#> LAPACK: /Library/Frameworks/R.framework/Versions/4.3-arm64/Resources/lib/libRlapack.dylib;  LAPACK version 3.11.0
#>
#> locale:
#> [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
#>
#> time zone: Europe/Oslo
#> tzcode source: internal
#>
#> attached base packages:
#> [1] stats     graphics  grDevices utils     datasets  methods   base
#>
#> other attached packages:
#> [1] BayesMallows_1.5.0.9000
#>
#> loaded via a namespace (and not attached):
#>  [1] styler_1.10.2     digest_0.6.33     fastmap_1.1.1     xfun_0.41
#>  [5] magrittr_2.0.3    glue_1.6.2        R.utils_2.12.3    knitr_1.45
#>  [9] htmltools_0.5.7   rmarkdown_2.25    lifecycle_1.0.4   Rdpack_2.6
#> [13] cli_3.6.1         R.methodsS3_1.8.2 vctrs_0.6.5       reprex_2.0.2
#> [17] withr_2.5.2       compiler_4.3.2    R.oo_1.25.0       R.cache_0.16.0
#> [21] rbibutils_2.2.16  purrr_1.0.2       rstudioapi_0.15.0 tools_4.3.2
#> [25] evaluate_0.23     Rcpp_1.0.11.2     yaml_2.3.7        rlang_1.1.2
#> [29] fs_1.6.3
mod <- compute_mallows(
  data = setup_rank_data(cluster_data),
  model_options = set_model_options(n_clusters = 3)
)
assess_convergence(mod)
```

![](https://i.imgur.com/rIyFDo9.png)<!-- -->

``` r
mod$burnin <- 100
plot(mod)
```

![](https://i.imgur.com/E7lQdNk.png)<!-- -->

``` r
compute_posterior_intervals(mod, parameter = "rho")
#>    parameter   cluster   item mean median hpdi central_interval
#> 1        rho Cluster 1 Item 1    1      1  [1]              [1]
#> 2        rho Cluster 1 Item 2    3      3  [3]              [3]
#> 3        rho Cluster 1 Item 3    4      4  [4]              [4]
#> 4        rho Cluster 1 Item 4    2      2  [2]              [2]
#> 5        rho Cluster 1 Item 5    5      5  [5]              [5]
#> 6        rho Cluster 2 Item 1    1      1  [1]              [1]
#> 7        rho Cluster 2 Item 2    2      2  [2]              [2]
#> 8        rho Cluster 2 Item 3    3      3  [3]              [3]
#> 9        rho Cluster 2 Item 4    4      4  [4]              [4]
#> 10       rho Cluster 2 Item 5    5      5  [5]              [5]
#> 11       rho Cluster 3 Item 1    1      1  [1]            [1,2]
#> 12       rho Cluster 3 Item 2    2      2  [2]            [1,2]
#> 13       rho Cluster 3 Item 3    5      5  [5]              [5]
#> 14       rho Cluster 3 Item 4    4      4  [4]            [3,4]
#> 15       rho Cluster 3 Item 5    3      3  [3]            [3,4]
compute_consensus(mod)
#>      cluster ranking   item   cumprob
#> 1  Cluster 1       1 Item 1 1.0000000
#> 2  Cluster 1       2 Item 4 1.0000000
#> 3  Cluster 1       3 Item 2 1.0000000
#> 4  Cluster 1       4 Item 3 1.0000000
#> 5  Cluster 1       5 Item 5 1.0000000
#> 6  Cluster 2       1 Item 1 1.0000000
#> 7  Cluster 2       2 Item 2 1.0000000
#> 8  Cluster 2       3 Item 3 1.0000000
#> 9  Cluster 2       4 Item 4 1.0000000
#> 10 Cluster 2       5 Item 5 1.0000000
#> 11 Cluster 3       1 Item 1 0.9573684
#> 12 Cluster 3       2 Item 2 1.0000000
#> 13 Cluster 3       3 Item 5 0.9563158
#> 14 Cluster 3       4 Item 4 0.9968421
#> 15 Cluster 3       5 Item 3 1.0000000
```

<sup>Created on 2023-12-11 with [reprex v2.0.2](https://reprex.tidyverse.org)</sup>
