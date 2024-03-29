% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/compute_exact_partition_function.R
\name{compute_exact_partition_function}
\alias{compute_exact_partition_function}
\title{Compute exact partition function}
\usage{
compute_exact_partition_function(
  alpha,
  n_items,
  metric = c("cayley", "hamming", "kendall")
)
}
\arguments{
\item{alpha}{Dispersion parameter.}

\item{n_items}{Number of items.}

\item{metric}{Distance function, one of "cayley", "hamming", or "kendall".}
}
\value{
The logarithm of the partition function.
}
\description{
For Cayley, Hamming, and Kendall distances, computationally
tractable functions are available for the exact partition function.
}
\examples{
compute_exact_partition_function(
  alpha = 3.4, n_items = 34, metric = "cayley"
)

compute_exact_partition_function(
  alpha = 3.4, n_items = 34, metric = "hamming"
)

compute_exact_partition_function(
  alpha = 3.4, n_items = 34, metric = "kendall"
)
}
\references{
\insertAllCited{}
}
\seealso{
Other partition function: 
\code{\link{estimate_partition_function}()},
\code{\link{get_cardinalities}()}
}
\concept{partition function}
