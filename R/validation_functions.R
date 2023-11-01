validate_class <- function(argument, class) {
  if(!inherits(argument, class)) {
    stop(paste(deparse(substitute(argument)), "must be an object of class",
               class, "."))
  }
}

validate_integer <- function(argument) {
  if(!is.numeric(argument) || argument < 1 || (round(argument) != argument)) {
    stop(paste(deparse(substitute(argument)), "must be a positive integer"))
  }
}

validate_positive <- function(argument) {
  if(argument <= 0 || !is.numeric(argument)) {
    stop(paste(deparse(substitute(argument)),
               "must be a strictly positive number of length one"))
  }
}

validate_positive_vector <- function(argument) {
  if(any(argument <= 0) || !is.numeric(argument)) {
    stop(paste(deparse(substitute(argument)),
               "must be a vector of strictly positive numbers"))
  }
}


validate_logical <- function(argument) {
  if(!is.logical(argument) || length(argument) != 1) {
    stop(paste(deparse(substitute(argument)),
               "must be a logical value of length one"))
  }
}

check_larger <- function(larger, smaller) {
  if(larger <= smaller) {
    stop(paste(deparse(substitute(larger)), "must be strictly larger than",
               deparse(substitute(smaller))))
  }
}

validate_preferences <- function(data, model) {
  if(inherits(data$preferences, "BayesMallowsIntransitive") &&
     model$error_model == "none") {
    stop("Intransitive pairwise comparisons. Please specify an error model.")
  }
}

validate_initial_values <- function(initial_values, data) {

  if(!is.null(initial_values$rho)) {
    if(length(unique(initial_values$rho)) != length(initial_values$rho)) {
      stop("initial value rho must be a ranking")
    }
    if(length(initial_values$rho) != data$n_items) {
      stop("initial value rho must have one value per item")
    }
  }


}
