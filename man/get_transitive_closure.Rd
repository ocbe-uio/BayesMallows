% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_transitive_closure.R
\name{get_transitive_closure}
\alias{get_transitive_closure}
\title{Get transitive closure}
\usage{
get_transitive_closure(rank_data)
}
\arguments{
\item{rank_data}{An object of class \code{"BayesMallowsData"} returned from
\link{setup_rank_data}.}
}
\value{
A dataframe with transitive closure, if there is any.
}
\description{
A simple method for showing any transitive closure computed by
\code{\link[=setup_rank_data]{setup_rank_data()}}.
}
\examples{
# Original beach preferences
head(beach_preferences)
dim(beach_preferences)
# We then create a rank data object
dat <- setup_rank_data(preferences = beach_preferences)
# The transitive closure contains additional filled-in preferences implied
# by the stated preferences.
head(get_transitive_closure(dat))
dim(get_transitive_closure(dat))

}
\seealso{
Other preprocessing: 
\code{\link{set_compute_options}()},
\code{\link{set_initial_values}()},
\code{\link{set_model_options}()},
\code{\link{set_priors}()},
\code{\link{set_progress_report}()},
\code{\link{set_smc_options}()},
\code{\link{setup_rank_data}()}
}
\concept{preprocessing}
