##' Find index of maximum, breaking ties at random and ignoring NAs.
##' @title Index of maximum with random ties
##' @param x Input vector.
##' @return Index of the maximum value.
##' @author Marvin N. Wright
which.max.random <- function(x) {
  if (all(is.na(x))) {
    return(NA)
  }
  which(rank(x, ties.method = "random", na.last = FALSE) == length(x))
}

##' Find index of minimum, breaking ties at random and ignoring NAs.
##' @title Index of minimum with random ties
##' @param x Input vector.
##' @return Index of the minimum value.
##' @author Marvin N. Wright
which.min.random <- function(x) {
  if (all(is.na(x))) {
    return(NA)
  }
  which(rank(x, ties.method = "random", na.last = TRUE) == 1)
}

##' Compute median survival if available or largest quantile available in all strata if median not available.
##' @title Compute largest available quantile
##' @param formula Formula for survival model.
##' @return Ordered factor levels
##' @author Marvin N. Wright
largest.quantile <- function(formula) {
  ## Fit survival model
  fit <- survfit(formula)
  smry <- summary(fit)
  
  ## Use median survival if available or largest quantile available in all strata if median not available
  max_quant <- max(aggregate(smry$surv ~ smry$strata, FUN = min)[, "smry$surv"])
  quantiles <- quantile(fit, conf.int = FALSE, prob = min(0.5, 1 - max_quant))[, 1]
  names(quantiles) <- gsub(".+=", "", names(quantiles))
  
  ## Return ordered levels
  names(sort(quantiles))
}

##' Reorder factor columns. Use mean for continuous response, class counts for factors and mean survival for survival response.
##' @title Reorder factor columns
##' @param data Data with factor columns.
##' @return Data with reordered factor columns.
##' @author Marvin N. Wright
##' @importFrom coin logrank_trafo
reorder.factor.columns <- function(data, multiclass_mode, survsort_mode) {
  ## Recode characters and unordered factors
  character.idx <- sapply(data[, -1], is.character)
  ordered.idx <- sapply(data[, -1], is.ordered)
  factor.idx <- sapply(data[, -1], is.factor)
  recode.idx <- character.idx | (factor.idx & !ordered.idx)
  
  ## Numeric response
  response <- data[, 1]
  if (is.factor(response)) {
    num.response <- as.numeric(response)
  } else if ("Surv" %in% class(response)) {
    num.response <- response[, 1]
  } else {
    num.response <- response
  }
  
  ## Recode each column
  data[, -1][, recode.idx] <- lapply(data[, -1][, recode.idx, drop = FALSE], function(x) {
    if ("Surv" %in% class(response)) {
      if (survsort_mode == "median") {
        ## Use median survival if available or largest quantile available in all strata if median not available
        levels.ordered <- largest.quantile(response ~ x)
      } else if (survsort_mode == "logrank") {
        ## Use log-rank scores to sort
        scores <- aggregate(response~x, FUN=logrank_trafo)
        levels.ordered <- scores$x[order(scores$response)]
      }
    } else if (is.factor(response) & nlevels(response) > 2) {
      if (multiclass_mode == "old") {
        means <- aggregate(num.response~x, FUN=mean)
        levels.ordered <- means$x[order(means$num.response)]
      } else if (multiclass_mode == "largest_class") {
        largest_class <- names(which.max.random(table(response)))
        means <- aggregate((response == largest_class) ~ x, FUN=mean)
        levels.ordered <- means$x[order(means[, 2])]
      } else {
        stop("Unknown multiclass mode.")
      }
    } else {
      ## Order factor levels by num.response
      means <- aggregate(num.response~x, FUN=mean)
      levels.ordered <- means$x[order(means$num.response)]
    }
    
    ## Return reordered factor
    factor(x, levels = levels.ordered, ordered = TRUE)
  })
  
  ## Return data
  data
}