#' @title SMC Processing
#' @param output a subset of an SMCMallows object (though technically any matrix will do)
#' @param colnames colnames
#' @return A processed file of the SMCMallows class
#' @seealso [smc_mallows_new_item_rank()] and
#' [smc_mallows_new_users()], which are functions generating objects
#' of SMCMallows class.
#' @noRd
#' @importFrom methods is
smc_processing <- function(output, colnames = NULL) {
  # Recasting of input for proper handling below
  df <- data.frame(data = output)

  # if colnames are specified, then incorporate them
  if (is.null(colnames)) {
    n_items <- ncol(df)
    cletters <- rep("Item", times = n_items)
    cindexes <- seq_len(n_items)
    cnames <- c(paste(cletters, cindexes, sep = " "))
    colnames(df) <- cnames
  } else {
    colnames(df) <- colnames
  }
  new_df <- stats::reshape(
    df,
    direction = "long",
    varying = names(df),
    new.row.names = seq_len(prod(dim(df))),
    v.names = "value",
    timevar = "item",
    times = names(df)
  )
  new_df$id <- NULL # the "id" should not be part of the SMCMallows object
  class(new_df) <- c("SMCMallows", "data.frame")
  return(new_df)
}


