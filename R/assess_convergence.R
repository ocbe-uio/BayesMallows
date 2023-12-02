#' Trace Plots from Metropolis-Hastings Algorithm
#'
#' `assess_convergence` provides trace plots for the parameters of the Mallows
#' Rank model, in order to study the convergence of the Metropolis-Hastings
#' algorithm.
#'
#' @param model_fit A fitted model object of class `BayesMallows` returned from
#'   [compute_mallows()] or an object of class `BayesMallowsMixtures` returned
#'   from [compute_mallows_mixtures()].
#'
#' @param parameter Character string specifying which parameter to plot.
#'   Available options are `"alpha"`, `"rho"`, `"Rtilde"`, `"cluster_probs"`, or
#'   `"theta"`.
#'
#' @param items The items to study in the diagnostic plot for `rho`. Either a
#'   vector of item names, corresponding to `model_fit$items` or a vector of
#'   indices. If NULL, five items are selected randomly. Only used when
#'   `parameter = "rho"` or `parameter = "Rtilde"`.
#'
#' @param assessors Numeric vector specifying the assessors to study in the
#'   diagnostic plot for `"Rtilde"`.
#'
#' @param ... Other arguments passed on to other methods. Currently not used.
#'
#' @export
#' @family diagnostics
#'
#' @example /inst/examples/assess_convergence_example.R
assess_convergence <- function(model_fit, ...) {
  UseMethod("assess_convergence")
}

#' @export
#' @rdname assess_convergence
assess_convergence.BayesMallows <- function(
    model_fit,
    parameter = c("alpha", "rho", "Rtilde", "cluster_probs", "theta"),
    items = NULL,
    assessors = NULL,
    ...) {
  parameter <- match.arg(
    parameter,
    c("alpha", "rho", "Rtilde", "cluster_probs", "theta")
  )

  if (parameter == "alpha") {
    trace_alpha(model_fit$alpha, FALSE)
  } else if (parameter == "rho") {
    trace_rho(model_fit, items)
  } else if (parameter == "Rtilde") {
    trace_rtilde(model_fit, items, assessors)
  } else if (parameter == "cluster_probs") {
    m <- model_fit$cluster_probs
    m$n_clusters <- model_fit$n_clusters
    trace_cluster_probs(m)
  } else if (parameter == "theta") {
    trace_theta(model_fit)
  }
}

#' @export
#' @rdname assess_convergence
assess_convergence.BayesMallowsMixtures <- function(
    model_fit,
    parameter = c("alpha", "cluster_probs"),
    items = NULL,
    assessors = NULL,
    ...) {
  parameter <- match.arg(parameter, c("alpha", "cluster_probs"))

  if (parameter == "alpha") {
    m <- do.call(rbind, lapply(model_fit, function(x) {
      x$alpha$cluster <- as.character(x$alpha$cluster)
      x$alpha$n_clusters <- x$n_clusters
      x$alpha
    }))
    trace_alpha(m, TRUE)
  } else if (parameter == "cluster_probs") {
    m <- do.call(rbind, lapply(model_fit, function(x) {
      x$cluster_probs$cluster <- as.character(x$cluster_probs$cluster)
      x$cluster_probs$n_clusters <- x$n_clusters
      x$cluster_probs
    }))
    trace_cluster_probs(m)
  }
}

trace_alpha <- function(m, clusters) {
  p <- ggplot2::ggplot(m, ggplot2::aes(
    x = .data$iteration, y = .data$value,
    group = interaction(.data$chain, .data$cluster),
    color = .data$cluster, linetype = .data$chain
  )) +
    ggplot2::geom_line() +
    ggplot2::xlab("Iteration") +
    ggplot2::ylab(expression(alpha)) +
    ggplot2::labs(color = "Cluster") +
    ggplot2::labs(linetype = "Chain")

  if (clusters) {
    p <- p +
      ggplot2::theme(legend.position = "none") +
      ggplot2::facet_wrap(ggplot2::vars(.data$n_clusters),
        labeller = ggplot2::as_labeller(cluster_labeler_function), scales = "free_y"
      )
  }
  return(p)
}

trace_rho <- function(model_fit, items, clusters = model_fit$n_clusters > 1) {
  if (is.null(items) && model_fit$n_items > 5) {
    message("Items not provided by user. Picking 5 at random.")
    items <- sample.int(model_fit$n_items, 5)
  } else if (is.null(items) && model_fit$n_items > 0) {
    items <- seq.int(from = 1, to = model_fit$n_items)
  } else if (!is.null(items)) {
    if(is.numeric(items) &&
       length(setdiff(items, seq_len(model_fit$n_item))) > 0) {
      stop("numeric items vector must contain indices between 1 and the number of items")
    }
    if(is.character(items) && length(setdiff(items, model_fit$items) > 0)) {
      stop("unknown items provided")
    }
  }

  if (!is.character(items)) {
    items <- model_fit$items[items]
  }

  df <- model_fit$rho[model_fit$rho$item %in% items, , drop = FALSE]

  p <- ggplot2::ggplot(
    df, ggplot2::aes(
      x = .data$iteration, y = .data$value, color = .data$item)) +
    ggplot2::geom_line() +
    ggplot2::theme(legend.title = ggplot2::element_blank()) +
    ggplot2::xlab("Iteration") +
    ggplot2::ylab(expression(rho))

  if (clusters) {
    p <- p + ggplot2::facet_wrap(ggplot2::vars(.data$cluster))
  } else {
    p <- p +
      ggplot2::facet_wrap(
        ggplot2::vars(.data$chain),
        labeller = ggplot2::as_labeller(function(x) paste("Chain", x))
      )
  }
  return(p)
}

trace_rtilde <- function(model_fit, items, assessors, ...) {
  if (!model_fit$save_aug) {
    stop("Please rerun with compute_mallows with save_aug = TRUE")
  }

  if (is.null(items) && model_fit$n_items > 5) {
    message("Items not provided by user. Picking 5 at random.")
    items <- sample.int(model_fit$n_items, 5)
  } else if (is.null(items) && model_fit$n_items > 0) {
    items <- seq.int(from = 1, to = model_fit$n_items)
  }

  if (is.null(assessors) && model_fit$n_assessors > 5) {
    message("Assessors not provided by user. Picking 5 at random.")
    assessors <- sample.int(model_fit$n_assessors, 5)
  } else if (is.null(assessors) && model_fit$n_assessors > 0) {
    assessors <- seq.int(from = 1, to = model_fit$n_assessors)
  } else if (!is.null(assessors)) {
    if (length(setdiff(assessors, seq(1, model_fit$n_assessors, 1))) > 0) {
      stop("assessors vector must contain numeric indices between 1 and the number of assessors")
    }
  }

  if (is.factor(model_fit$augmented_data$item) && is.numeric(items)) {
    items <- levels(model_fit$augmented_data$item)[items]
  }

  df <- model_fit$augmented_data[
    model_fit$augmented_data$assessor %in% assessors &
      model_fit$augmented_data$item %in% items, ,
    drop = FALSE
  ]

  df$assessor <- as.factor(df$assessor)
  levels(df$assessor) <- paste("Assessor", levels(df$assessor))
  df$chain <- as.factor(df$chain)
  levels(df$chain) <- paste("Chain", levels(df$chain))

  ggplot2::ggplot(df, ggplot2::aes(x = .data$iteration, y = .data$value, color = .data$item)) +
    ggplot2::geom_line() +
    ggplot2::facet_wrap(ggplot2::vars(.data$assessor, .data$chain)) +
    ggplot2::theme(legend.title = ggplot2::element_blank()) +
    ggplot2::xlab("Iteration") +
    ggplot2::ylab("Rtilde")
}

trace_cluster_probs <- function(m) {
  ggplot2::ggplot(m, ggplot2::aes(
    x = .data$iteration, y = .data$value,
    color = .data$cluster
  )) +
    ggplot2::geom_line() +
    ggplot2::theme(legend.position = "none") +
    ggplot2::xlab("Iteration") +
    ggplot2::ylab(expression(tau[c])) +
    ggplot2::facet_wrap(ggplot2::vars(.data$n_clusters),
      labeller = ggplot2::as_labeller(cluster_labeler_function), scales = "free_y"
    )
}

trace_theta <- function(model_fit) {
  if (is.null(model_fit$theta) || length(model_fit$theta) == 0) {
    stop("Theta not available. Run compute_mallows with error_model = 'bernoulli'.")
  }
  p <- ggplot2::ggplot(model_fit$theta, ggplot2::aes(x = .data$iteration, y = .data$value)) +
    ggplot2::xlab("Iteration") +
    ggplot2::ylab(expression(theta)) +
    ggplot2::geom_line()

  return(p)
}

cluster_labeler_function <- function(n_clusters) {
  paste(n_clusters, ifelse(n_clusters == 1, "cluster", "clusters"))
}
