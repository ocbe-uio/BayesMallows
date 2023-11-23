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
