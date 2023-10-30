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
    stop(paste(deparse(substitute(argument)), "must be a positive number"))
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
