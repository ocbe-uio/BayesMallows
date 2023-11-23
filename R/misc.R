#' @importFrom Rdpack reprompt
#' @importFrom Rcpp sourceCpp
#' @importFrom stats aggregate
#' @importFrom utils head
#' @useDynLib BayesMallows, .registration = TRUE
NULL


.onUnload <- function(libpath) {
  library.dynam.unload("BayesMallows", libpath)
}

#' Check if a vector is a permutation
#' @param vec a vector
#' @return TRUE if vec is a permutation
#' @noRd
validate_permutation <- function(vec) {
  if (!any(is.na(vec))) {
    return(all(sort(vec) == seq_along(vec)))
  } else if (all(is.na(vec))) {
    return(TRUE)
  } else {
    return(all(vec[!is.na(vec)] <= length(vec)) &&
      all(vec[!is.na(vec)] >= 1) && !any(duplicated(vec[!is.na(vec)])))
  }
}

#' Prepare partition functions
#'
#' Utility function for estimating partition function of the Mallows model.
#'
#' @param logz_estimate Optional argument containing the result of calling
#'   \code{\link{estimate_partition_function}}.
#' @param metric Metric to be used.
#' @param n_items Number of items.
#'
#' @return List with two elements, \code{cardinalities} and \code{logz_estimate},
#' one of which is \code{NULL} and the other of which contains partition function
#' estimates.
#'
#' @export
#' @family preprocessing
#'
prepare_partition_function <- function(logz_estimate = NULL, metric, n_items) {
  # First, has the user supplied an estimate?
  if (!is.null(logz_estimate)) {
    return(list(cardinalities = NULL, logz_estimate = logz_estimate))
  }

  # Second, do we have a sequence?
  relevant_params <- partition_function_data[partition_function_data$n_items == n_items &
    partition_function_data$metric == metric &
    partition_function_data$type == "cardinalities", , drop = FALSE]

  if (nrow(relevant_params) == 1) {
    return(list(cardinalities = unlist(relevant_params$values), logz_estimate = NULL))
  }

  # Third, do we have an importance sampling estimate?
  relevant_params <- partition_function_data[partition_function_data$n_items == n_items &
    partition_function_data$metric == metric &
    partition_function_data$type == "importance_sampling", , drop = FALSE]

  if (nrow(relevant_params) == 1) {
    return(list(cardinalities = NULL, logz_estimate = unlist(relevant_params$values)))
  }

  # Fourth, is it the Ulam distance?
  if (metric == "ulam") {
    message("Exact partition function no longer available for Ulam distance with >95 items.")
  }

  # Fifth, can we compute the partition function in our C++ code?
  if (metric %in% c("cayley", "hamming", "kendall")) {
    return(list(cardinalities = NULL, logz_estimate = NULL))
  }

  stop("Partition function not available. Please compute an estimate using estimate_partition_function().")
}

# function taken from PLMIX package:
# Copyright Cristina Mollica and Luca Tardella
unit_to_freq <- function(data) {
  data <- fill_single_entries(data = data)
  K <- ncol(data)
  freq <- table(apply(data, 1, paste, collapse = "-"))
  obs_seq <- matrix(as.numeric(unlist(strsplit(names(freq),
    split = "-"
  ))), nrow = length(freq), ncol = K, byrow = TRUE)
  rownames(obs_seq) <- NULL
  out <- cbind(obs_seq, freq = freq, deparse.level = 0)
  rownames(out) <- NULL
  return(out)
}

# function taken from PLMIX package:
# Copyright Cristina Mollica and Luca Tardella
fill_single_entries <- function(data) {
  if (is.vector(data)) {
    data <- t(data)
  }
  K <- ncol(data)
  r_single_miss <- (rowSums(data == 0) == 1)
  if (any(r_single_miss)) {
    w_row <- which(r_single_miss)
    w_col <- apply(data[w_row, , drop = FALSE], 1, function(x) {
      which(x ==
        0)
    })
    w_item <- apply(data[w_row, , drop = FALSE], 1, setdiff,
      x = 1:K
    )
    data[cbind(w_row, w_col)] <- w_item
    warning(paste(paste0("Top-", K - 1, ""), "sequencies correspond to full orderings. Single missing entries filled."),
      call. = FALSE
    )
  }
  return(data)
}

## Source: https://stackoverflow.com/questions/11095992/generating-all-distinct-permutations-of-a-list-in-r
permutations <- function(n) {
  if (n == 1) {
    return(matrix(1))
  } else {
    sp <- permutations(n - 1)
    p <- nrow(sp)
    A <- matrix(nrow = n * p, ncol = n)
    for (i in 1:n) {
      A[(i - 1) * p + 1:p, ] <- cbind(i, sp + (sp >= i))
    }
    return(A)
  }
}
