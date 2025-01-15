#' Perform OLS regression
#'
#' This function performs an ordinary least squares (OLS) regression on a dataset
#' that contains a response variable and one or more explanatory variables.
#'
#' @param data A data frame containing the response variable and explanatory variables as columns.
#' @param response_variable The column name (as a string) in `data` that represents the response variable.
#' @param explanatory_variables A character vector specifying the names of the columns in `data` that represent the explanatory variables (e.g. c("x", "z")).
#' @return A list containing:
#'   \item{summary_table}{A data frame with the regression results, including:
#'     \itemize{
#'       \item \code{Estimate}: The estimated coefficients.
#'       \item \code{Std. Error}: The standard errors for the coefficients.
#'       \item \code{t-value}: The t-statistics for each coefficient.
#'       \item \code{p-value}: The p-values for each coefficient.
#'     }}
#'   \item{R_squared}{The R-squared value.}
#'   \item{adjusted_R_squared}{The adjusted R-squared value.}
#'   \item{model_p_value}{The p-value for the overall model significance.}
#' @examples
#' # Example usage:
#' data <- data.frame(x = rnorm(100, 30000, 5000), y = rnorm(100, 45, 10))
#' ols_model_function(data, response_variable = "y", explanatory_variables = c("x"))
#'
#' @export

ols_model_function <- function(data, response_variable, explanatory_variables) {
  
  # Checking that data contains all input variables 
  if (!all(c(response_variable, explanatory_variables) %in% colnames(data))) {
    stop("The data and the entered variables do not match.")
  }
  
  # Creating the X matrix and the y vector
  X_matrix <- as.matrix(data[, explanatory_variables, drop = FALSE])  
  X_matrix <- cbind(1, X_matrix) 
  y_vector <- as.matrix(data[[response_variable]], ncol = 1) 
  
  # Calculating OLS coefficients
  XtX <- t(X_matrix) %*% X_matrix
  Xty <- t(X_matrix) %*% y_vector 
  coefficients <- solve(XtX) %*% Xty 
  
  # Calculating residuals
  predictions <- X_matrix %*% coefficients
  residuals <- y_vector - predictions
  
  # Calculating R², adjusted R², and p-value for the model
  sst <- sum((y_vector - mean(y_vector))^2) 
  sse <- sum(residuals^2)                   
  R_squared <- 1 - (sse / sst)              
  n <- nrow(X_matrix)  
  k <- ncol(X_matrix)  
  adjusted_R_squared <- 1 - ((1 - R_squared) * (n - 1) / (n - k))
  f_stat <- ((sst - sse) / (k - 1)) / (sse / (n - k))
  model_p_value <- pf(f_stat, df1 = k - 1, df2 = n - k, lower.tail = FALSE)
  
  # Calculating standard errors, t-values, and p-values for explanatory variables
  residual_variance <- sse / (n - k)  
  coefficient_variances <- diag(residual_variance * solve(XtX))
  standard_errors <- sqrt(coefficient_variances)
  t_values <- as.vector(coefficients / standard_errors)
  p_values <- 2 * pt(-abs(t_values), df = n - k)
  
  # Rounding values
  R_squared <- round(R_squared, 4)
  adjusted_R_squared <- round(adjusted_R_squared, 4)
  model_p_value <- round(model_p_value, 4)
  coefficients <- round(coefficients, 4)
  standard_errors <- round(standard_errors, 4)
  t_values <- round(t_values, 4)
  p_values <- round(p_values, 4)
  
  # Combining results into a summary table
  summary_table <- data.frame(
    table_coefficients = as.vector(coefficients),
    standard_errors = standard_errors,
    t_values = t_values,
    p_values = p_values,
    row.names = c("Intercept", explanatory_variables))
  colnames(summary_table) <- c("Estimate", "Std. Error", "t-value", "p-value")
  
  # Printing model summary statistics and returning the summary table
  results <- list(summary_table = summary_table)
  cat(sprintf("Model: R-squared: %.4f | Adjusted R-squared: %.4f | p-value: %.4g\n", 
              R_squared, adjusted_R_squared, model_p_value))
  return(results)
}