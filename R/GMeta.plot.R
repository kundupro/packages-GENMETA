#' Generalized Meta-analysis(forest plot)
#'
#' This function plots the confidence intervals with boxes as the study specific estimates and diamond as the GMeta estimate. For the current version, it assumes that the estimate of the variance-covariance matrix in each of the studies is provided.
#' It is demonstrated using a different dataset, "study_info_plot", which meets the assumption.
#' @param x an object of class "GMeta"
#' @param study_info_plot a list of lists containing information about the studies(similar to the study_info argument used in GMeta function.)
#' @examples
#' # This example shows how to obtain the forest plot of GMeta object.
#' data(study_info_plot)
#' data(reference_data)
#' model <- "logistic"
#' result_diff <- GMeta(study_info_plot, reference_data, model, variable_intercepts = TRUE)
#' GMeta.plot(result_diff, study_info_plot)
#' @export

GMeta.plot <- function(x, study_info_plot)
{
  no_of_studies <- length(study_info_plot)
  row_names <- c()
  tot_var_names <- c()
  for(i in 1: no_of_studies)
  {
    tot_var_names <- union(tot_var_names, names(study_info_plot[[i]][[1]]))
  }
  tot_var_names_wo_intercept <- tot_var_names[-1]

  y_indices <- seq((no_of_studies+1), 1, -1)

  for(i in 1:no_of_studies)
  {
    row_names <- c(row_names, paste0("Study",i))
  }
  row_names <- c(row_names, "GMeta")


  #op <- par(ask=TRUE)
  for(i in 1:length(tot_var_names_wo_intercept))
  {

    data_plot <- matrix(NA, (no_of_studies+1), 4)
    data_plot[,1] <- y_indices

    for(j in 1: no_of_studies)
    {
      if(tot_var_names_wo_intercept[i] %in% names(study_info_plot[[j]][[1]]) == T)
      {
        index <- which(names(study_info_plot[[j]][[1]]) %in% tot_var_names_wo_intercept[i] == T)
        data_plot[j,2] <- as.vector(study_info_plot[[j]][[1]])[index]
        data_plot[j,3] <- as.vector(study_info_plot[[j]][[1]])[index] - 1.96*sqrt(as.vector(diag(study_info_plot[[j]][[2]]))[index])
        data_plot[j,4] <- as.vector(study_info_plot[[j]][[1]])[index] + 1.96*sqrt(as.vector(diag(study_info_plot[[j]][[2]]))[index])
      }

    }
    index <- which(names(x[[1]]) %in% tot_var_names_wo_intercept[i] == T)
    data_plot[(no_of_studies+1),2] <- as.vector(x[[1]])[index]
    data_plot[(no_of_studies+1),3] <- as.vector(x[[1]])[index] - 1.96*sqrt(as.vector(diag(x[[2]])))[index]
    data_plot[(no_of_studies+1),4] <- as.vector(x[[1]])[index] + 1.96*sqrt(as.vector(diag(x[[2]])))[index]



    data_plot <- as.data.frame(data_plot)
    colnames(data_plot) <- c("Y", "coef","low","high")
    rownames(data_plot) <- row_names
    #row_names <- list(list(row_names))
    #Enter the forestplot command
    #forestplot(row_names, data_plot$coef, data_plot$low, data_plot$high, zero = 1, cex  = 2, lineheight = "auto", xlab = "log odds ratio")
    #plot.new()

    plot(data_plot[,2], data_plot[,1], axes=F, xlab = "", ylab = "", pch = 15, xlim = c(-1.5,1.5))
    for(k in 1:no_of_studies)
    {
      segments(data_plot[k,3], data_plot[k,1], x1=data_plot[k,4], y1=data_plot[k,1], col = "black")
      #points(data_plot[k,2], data_plot[k,1], pch = 18, cex=2)
    }
    segments(data_plot[(no_of_studies+1),3], data_plot[(no_of_studies+1),1], x1=data_plot[(no_of_studies+1),4], y1=data_plot[(no_of_studies+1),1], col = "black")
    points(data_plot[(no_of_studies+1),2], data_plot[(no_of_studies+1),1], pch = 18, cex=2)
    mtext("Log Odds Ratio", side=1, las=0, col="black", outer = F, line = 2.5)
    abline(v=0, lty = 3)
    axis(1, at = seq(-2,2,0.5), labels = seq(-2,2,0.5))
    axis(2, at = data_plot[,1], labels = rownames(data_plot), las = 2)
    title(main = tot_var_names_wo_intercept[i])

    #legend("topright", legend = c("Study specific estimates", "GMeta estimate", "95% CI"), lty = c(NA,NA,1), pch = c(15,18,NA))
    readline(prompt="Press [enter] for the forestplot of next covariate:")
  }
}
