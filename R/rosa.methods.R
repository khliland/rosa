
#' Predictions from ROSA models
#'
#' @param object A \code{rosa} object.
#' @param newdata Optional new data with the same types of predictor blocks as the ones used for fitting the object.
#' @param ncomp An \code{integer} vector giving the number of components to use for prediction.
#' @param comps An \code{integer} vector giving the exact components to use for predictions.
#' @param type A \code{character} indicating if responses or scores should be predicted (default = "response", or "scores")
#' @param na.action Function determining what to do with missing values in \code{newdata}.
#' @param ... Additional arguments. Currently not implemented.
#'
#' @return Returns predicted responses or scores depending on inputs.
#' @export
#'
#' @examples
predict.rosa <- function(object, newdata, ncomp = 1:object$ncomp, comps,
                         type = c("response", "scores"), na.action = na.pass, ...){
  if (missing(newdata) || is.null(newdata)){
    newX <- object$X.concat
    # newX <- model.matrix(object)
  } else {
    if(is.matrix(newdata)){
      if (ncol(newdata) != length(object$Xmeans))
        stop("'newdata' does not have the correct number of columns")
      newX <- newdata
    } else { # Assume newdata is a list
      newX <- do.call(cbind,newdata)
      if(ncol(newX) != length(object$Xmeans))
        stop("'newdata' does not have the correct number of columns")
    }
  }

  nobs <- dim(newX)[1]

  if (!is.null(object$scale)) newX <- newX / rep(object$scale, each = nobs)
  type <- match.arg(type)

  if (type == "response") {
    if (missing(comps) || is.null(comps)) {
      ## Predict with models containing ncomp[1] components,
      ## ncomp[2] components, etc.
      if (missing(newdata)) return(fitted(object)[,,ncomp, drop=FALSE])
      B <- coef(object, ncomp = ncomp, intercept = TRUE)
      dPred <- dim(B)
      dPred[1] <- dim(newX)[1]
      dnPred <- dimnames(B)
      dnPred[1] <-
        if(is.null(dimnames(newX))) list(NULL) else dimnames(newX)[1]
      pred <- array(dim = dPred, dimnames = dnPred)
      for (i in seq(along = ncomp))
        pred[,,i] <- newX %*% B[-1,,i] + rep(B[1,,i], each = nobs)
      return(pred)
    } else {
      ## Predict with a model containing the components `comps'
      B <- rowSums(coef(object, comps = comps), dims = 2)
      B0 <- object$Ymeans# - object$Xmeans %*% B
      pred <- newX %*% B + rep(c(B0), each = nobs)
      if (missing(newdata) && !is.null(object$na.action))
        pred <- napredict(object$na.action, pred)
      return(pred)
    }
  } else {
    ## Return predicted scores (for scores, `cumulative' has no meaning)
    ## When predicting scores, we allow ncomp as an alias for comps:
    if (missing(comps) || is.null(comps)) comps <- ncomp
    if (missing(newdata)) {
      TT <- object$scores[,comps]
      if (!is.null(object$na.action))  TT <- napredict(object$na.action, TT)
    } else {
      if (is.null(object$projection))
        stop("`object' has no `projection' component.  Maybe it was fitted with `stripped = TRUE'.")
      TT <- (newX - rep(object$Xmeans, each = nobs)) %*%
        object$projection[,comps]
    }
    return(TT)
  }
}

#' Regression coefficients from ROSA model.
#'
#' @param object A \code{rosa} object.
#' @param ncomp An \code{integer} vector giving the number of components to use for prediction.
#' @param comps An \code{integer} vector giving the exact components to use for predictions.
#' @param intercept A \code{logical} indicating if coefficients for the intercept should be included (default = FALSE).
#' @param ... Additional arguments. Currently not implemented.
#'
#' @return
#' @export
#'
#' @examples
coef.rosa <- function(object, ncomp = object$ncomp, comps, intercept = FALSE,
                     ...)
{
  if (missing(comps) || is.null(comps)) {
    ## Cumulative coefficients:
    B <- object$coefficients[,,ncomp, drop=FALSE]
    if (intercept == TRUE) {      # Intercept has only meaning for
      # cumulative coefficients
      dB <- dim(B)
      dB[1] <- dB[1] + 1
      dnB <- dimnames(B)
      dnB[[1]] <- c("(Intercept)", dnB[[1]])
      BInt <- array(dim = dB, dimnames = dnB)
      BInt[-1,,] <- B
      for (i in seq(along = ncomp))
        BInt[1,,i] <- object$Ymeans - object$Xmeans %*% B[,,i]
      B <- BInt
    }
  } else {
    ## Individual coefficients:
    B <- object$coefficients[,,comps, drop=FALSE]
    g1 <- which(comps > 1)
    ## Indiv. coef. must be calculated since object$coefficients is
    ## cumulative coefs.
    B[,,g1] <- B[,,g1, drop=FALSE] -
      object$coefficients[,,comps[g1] - 1, drop=FALSE]
    dimnames(B)[[3]] <- paste("Comp", comps)
  }
  return(B)
}

