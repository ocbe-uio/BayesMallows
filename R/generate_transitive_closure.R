generate_transitive_closure <- function(preferences, cl = NULL) {
  if (is.null(preferences)) {
    return(NULL)
  }
  stopifnot(is.null(cl) || inherits(cl, "cluster"))

  if (!is.numeric(preferences$assessor)) {
    stop("assessor column in preferences must be numeric")
  }

  prefs <- splitpref(preferences)
  if (is.null(cl)) {
    lapplyfun <- lapply
  } else {
    lapplyfun <- function(X, FUN, ...) {
      parallel::parLapply(cl = cl, X = X, fun = FUN, ...)
    }
  }
  prefs <- lapplyfun(seq_along(prefs), function(i) {
    cbind(
      assessor = as.numeric(names(prefs)[[i]]),
      .generate_transitive_closure(prefs[[i]])
    )
  })

  prefs <- do.call(rbind.data.frame, prefs)

  # Check if there are any inconsistencies
  check <- merge(prefs, prefs,
    by.x = c("assessor", "bottom_item", "top_item"),
    by.y = c("assessor", "top_item", "bottom_item")
  )
  if (nrow(check) > 0) {
    class(preferences) <- c("BayesMallowsIntransitive", class(preferences))
    return(preferences)
  } else {
    class(prefs) <- c("BayesMallowsTransitiveClosure", class(prefs))
    return(prefs)
  }
}


.generate_transitive_closure <- function(mat) {
  # This line was an answer to StackOverflow question 51794127
  my_set <- do.call(sets::set, apply(mat, 1, sets::as.tuple))
  r <- relations::endorelation(graph = my_set)
  tc <- relations::transitive_closure(r)
  incidence <- relations::relation_incidence(tc)

  new_mat <- which(incidence == 1, arr.ind = TRUE)
  row_inds <- as.numeric(gsub("[^0-9]+", "", rownames(incidence)))
  result <- data.frame(
    bottom_item = row_inds[new_mat[, 1, drop = FALSE]],
    top_item = row_inds[new_mat[, 2, drop = FALSE]]
  )
  result
}
