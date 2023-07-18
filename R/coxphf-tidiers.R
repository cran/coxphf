
#' Tidy a coxphf object
#' 
#' @param x A `coxphf` object.
#' @param conf.int Logical indicating whether or not to include 
#'   a confidence interval in the tidied output. Defaults to FALSE.
#' @param conf.level The confidence level to use for the confidence 
#'   interval if conf.int = TRUE. Must be strictly greater than 0 
#'   and less than 1. Defaults to 0.95, which corresponds to a 
#'   95 percent confidence interval.
#' @param exponentiate Logical indicating whether or not to display coefficient
#'   estimates on an exponential scale.
#' @param ... Unused, included for generic consistency only.
#' @return A tidy [tibble::tibble()] summarizing component-level
#'   information about the model
#' 
#' import tibble
#' @export

tidy.coxphf <- function(x, conf.int = FALSE, conf.level = 0.95, exponentiate = FALSE, ...){
  
  result <- with(
    x,
    tibble::tibble(
      term = names(coefficients),
      estimate = coefficients,
      std.error = sqrt(diag(var)),
      statistic = qchisq(1 - prob, 1),
      p.value = prob
    )
  )
  
  if(exponentiate){
    result$estimate <- exp(result$estimate)
  }
  
  if(conf.int){
    ci <- stats::confint(x, level = conf.level)
    colnames(ci) <- c("conf.low", "conf.high")
    if(exponentiate) ci <- exp(ci)
    result <- cbind(
      result,
      ci
    ) |> tibble::tibble()
  }
  
  
  return(result)
}

#' Glance at a coxphf object
#' 
#' @param x A `coxphf` object.
#' @param ... Unused, included for generic consistency only.
#' @return A tidy [tibble::tibble()] summarizing component-level
#'   information about the model
#' 
#' @export
glance.coxphf <- function(x, ...){
  
  result <- with(
    x,
    tibble::tibble(
      n = n,
      nevent = sum(y[, "status"], na.rm = TRUE),
      logLik = as.numeric(stats::logLik(x)),
      AIC = stats::AIC(x)
    )
  )
  
  return(result)
  
}

#' Augment a coxphf object
#' 
#' @param x A `coxphf` object.
#' @param data the dataset used to fit the model.
#' @param ... Unused, included for generic consistency only.
#' @return A tidy [tibble::tibble()] summarizing component-level
#'   information about the model
#' 
#' @export
augment.coxphf <- function(x, data = x$y, ...){
  
  result <- cbind(
    data,
    ".linear.predictor" = x$linear.predictor
  ) |> tibble::as_tibble()
  
  return(result)
  
}

#' @importFrom generics glance 
#' @export
generics::glance






