#' Summarizing Generalized Meta Analysis
#'
#' This function prints the summary of GMeta results.
#' @param object an object of class "GMeta"
#' @param signi_digits an optional numeric indicating the number of significant digits to be shown in the summary. Default is 3.
#' @examples
#' # This example shows how to obtain the summary of GMeta object.
#' data(reference_data)
#' data(study_info)
#' model <- "logistic"
#' result_diff <- GMeta(study_info, reference_data, model, variable_intercepts = TRUE)
#' GMeta.summary(result_diff, signi_digits = 4)
#' result_same <- GMeta(study_info, reference_data, model)
#' GMeta.summary(result_same)
#' @export
# summary.GMeta <-function(object, ...){
#   UseMethod("summary")
#   NextMethod("generic = NULL, object = NULL", ...)
# }
GMeta.summary <- function(object, signi_digits = 3)
{
  x <- object
  GMeta_opt_estimate <- as.vector(x[[1]])
  GMeta_opt_std_error <- sqrt(as.vector(diag(x[[2]])))
  Var_opt_GMeta <- diag(GMeta_opt_std_error)
  z_stat_opt_GMeta <- GMeta_opt_estimate/GMeta_opt_std_error
  p_val_opt_GMeta <- 1 - pnorm(abs(z_stat_opt_GMeta))
  summary_data_frame_opt <- data.frame(cbind(GMeta_opt_estimate, GMeta_opt_std_error, z_stat_opt_GMeta, p_val_opt_GMeta))
  colnames(summary_data_frame_opt) <- c("Estimate", "Std.Error", "z value", "Pr(>|z|)")
  rownames(summary_data_frame_opt) <- names(x[[1]])
  summary_data_frame_opt <- signif(summary_data_frame_opt, signi_digits)
  summary_data_frame_opt[,4] <- noquote(sapply(summary_data_frame_opt[, 4], sign.star))
  cat("Call:\n")
  cat("GMeta(study_info, reference_data, model)\n")
  cat("\n")
  cat("Coefficients: \n")
  print(summary_data_frame_opt)
  cat("\n---\n")
  cat("Significant codes:\n")
  cat("0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 \n")
  cat("\n")
  cat("Total number of iterations: \n")
  print(x[[4]])

}
