Transitive Closures
================

Here I investigate the results of `generate_transitive_closure` for the beach data.

``` r
library(dplyr)
library(BayesMallows)
```

First, I generate the transitive closure.

``` r
beach_tc <- generate_transitive_closure(beach_preferences)
```

Assessor 1
----------

### Item 2

Let us look at how assessor 1 has ranked item 2.

``` r
beach_preferences %>% 
  filter(assessor == 1, 
         bottom_item == 2 | top_item == 2) %>% 
  knitr::kable()
```

|  assessor|  bottom\_item|  top\_item|
|---------:|-------------:|----------:|
|         1|             2|         15|

Assessor 1 has stated that she prefers item 15 to item 2.

Next, let us look at all rows involving item 2 in `beach_tc`:

``` r
beach_tc %>% 
  filter(assessor == 1, 
         bottom_item == 2 | top_item == 2) %>% 
  knitr::kable()
```

|  assessor|  bottom\_item|  top\_item|
|---------:|-------------:|----------:|
|         1|             2|          6|
|         1|             2|         15|

The ordering {*A*<sub>2</sub> ≺ *A*<sub>6</sub>} is implied by the other orderings.

Why is it so? Let us look at items 2, 6, and 15 in the original stated preferences (not the transitive closure):

``` r
beach_preferences %>% 
  filter(assessor == 1, 
         bottom_item %in% c(2, 6, 15) | top_item %in% c(2, 6, 15)) %>% 
  knitr::kable()
```

|  assessor|  bottom\_item|  top\_item|
|---------:|-------------:|----------:|
|         1|             2|         15|
|         1|             5|         15|
|         1|            12|          6|
|         1|            14|          6|
|         1|             3|          6|
|         1|            15|          6|
|         1|             9|          6|

We see that items 3, 9, 12, 14, and 15 are all disfavored to item 6. Since item 15 is preferred to item 2, as we saw above, it follows that item 6 is favored to item 2.

So the logic here is that {*A*<sub>2</sub> ≺ *A*<sub>15</sub>, *A*<sub>15</sub> ≺ *A*<sub>6</sub>} implies {*A*<sub>2</sub> ≺ *A*<sub>6</sub>}.

### Item 4

Let us look at how assessor 1 has ranked item 4.

``` r
beach_preferences %>% 
  filter(assessor == 1, 
         bottom_item == 4 | top_item == 4) %>% 
  knitr::kable()
```

|  assessor|  bottom\_item|  top\_item|
|---------:|-------------:|----------:|
|         1|             4|          7|

Assessor 1 has stated that she prefers item 7 to item 4.

Next, let us look at all rows involving item 4 in `beach_tc`:

``` r
beach_tc %>% 
  filter(assessor == 1, 
         bottom_item == 4 | top_item == 4) %>% 
  knitr::kable()
```

|  assessor|  bottom\_item|  top\_item|
|---------:|-------------:|----------:|
|         1|             4|          1|
|         1|             4|          3|
|         1|             4|          6|
|         1|             4|          7|
|         1|             4|          9|
|         1|             4|         10|
|         1|             4|         11|

The orderings {*A*<sub>4</sub> ≺ *A*<sub>1</sub>}, {*A*<sub>4</sub> ≺ *A*<sub>3</sub>}, {*A*<sub>4</sub> ≺ *A*<sub>6</sub>}, {*A*<sub>4</sub> ≺ *A*<sub>9</sub>}, {*A*<sub>4</sub> ≺ *A*<sub>10</sub>}, and {*A*<sub>4</sub> ≺ *A*<sub>11</sub>} are implied by the other orderings.

Why is it so? Let us look at items 4, 7, 1, 3, 6, 9, 10, and 11 in the original stated preferences (not the transitive closure):

``` r
items <- c(1, 3, 4, 6, 7, 9, 10, 11)
beach_preferences %>% 
  filter(assessor == 1, 
         bottom_item %in% items | top_item %in% items) %>% 
  arrange(bottom_item, top_item) %>% 
  knitr::kable()
```

|  assessor|  bottom\_item|  top\_item|
|---------:|-------------:|----------:|
|         1|             1|         10|
|         1|             3|          6|
|         1|             3|         11|
|         1|             4|          7|
|         1|             5|          3|
|         1|             5|          7|
|         1|             5|         10|
|         1|             7|          1|
|         1|             7|          9|
|         1|             8|         10|
|         1|             9|          3|
|         1|             9|          6|
|         1|             9|         11|
|         1|            12|          6|
|         1|            13|          3|
|         1|            13|         11|
|         1|            14|          6|
|         1|            14|          9|
|         1|            14|         11|
|         1|            15|          6|

Here is the logic:

We have {*A*<sub>4</sub> ≺ *A*<sub>7</sub>}.

We also have {*A*<sub>7</sub> ≺ *A*<sub>1</sub>} and {*A*<sub>7</sub> ≺ *A*<sub>9</sub>}, which imply {*A*<sub>4</sub> ≺ *A*<sub>1</sub>} and {*A*<sub>4</sub> ≺ *A*<sub>9</sub>}.

Further, we have {*A*<sub>1</sub> ≺ *A*<sub>10</sub>}, which implies {*A*<sub>4</sub> ≺ *A*<sub>10</sub>} because {*A*<sub>4</sub> ≺ *A*<sub>1</sub>}.

We also have {*A*<sub>9</sub> ≺ *A*<sub>3</sub>}, {*A*<sub>9</sub> ≺ *A*<sub>6</sub>}, and {*A*<sub>9</sub> ≺ *A*<sub>11</sub>}, which impliy {*A*<sub>4</sub> ≺ *A*<sub>3</sub>}, {*A*<sub>4</sub> ≺ *A*<sub>6</sub>}, and {*A*<sub>4</sub> ≺ *A*<sub>11</sub>}, because {*A*<sub>4</sub> ≺ *A*<sub>9</sub>}.

That was it!
