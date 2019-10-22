#' Print method for ROSA
#'
#' @param x An object fitted by \code{rosa}
#' @param ... Unused.
#'
#' @return Returns \code{x} invisibly.
#' @export
print.rosa <- function(x, ...) {
     ana <- "Response Orinented Sequential Alternation"
     alg <- "CPPLS"
  cat(ana, ", fitted with the", alg, "algorithm.")
  if (!is.null(x$validation))
    cat("\nCross-validated using", length(x$validation$segments),
        attr(x$validation$segments, "type"), "segments.")
  cat("\nCall:\n", deparse(x$call), "\n", sep = "")
  invisible(x)
}

#' Summary method for ROSA objects
#'
#' @param object A \code{rosa} object.
#' @param what A \code{character} indicating if summary should include all, validation or training.
#' @param digits The number of digits used for printing.
#' @param print.gap Gap between columns when printing.
#' @param ...
#'
#' @return No return
#' @export
#' @importFrom pls RMSEP MSEP mvrValstats R2
#'
#' @examples
summary.rosa <- function(object, what = c("all", "validation", "training"),
                        digits = 4, print.gap = 2, ...)
{
  what <- match.arg(what)
  if (what == "all") what <- c("validation", "training")
  if (is.null(object$validation)) what <- "training"

  nobj <- nrow(object$scores)
  nresp <- length(object$Ymeans)
  yvarnames <- dimnames(fitted(object))[[2]]
  cat("Data: \tX dimension:", nobj, length(object$Xmeans),
      "\n\tY dimension:", nobj, nresp)
  cat("\nFit method:", object$method)
  cat("\nNumber of components considered:", object$ncomp)

  for (wh in what) {
    if (wh == "training") {
      cat("\nTRAINING: % variance explained\n")
      xve <- explvar(object)
      yve <- 100 * drop(R2(object, estimate = "train",
                           intercept = FALSE)$val)
      tbl <- rbind(cumsum(xve), yve)
      dimnames(tbl) <- list(c("X", yvarnames),
                            paste(1:object$ncomp, "comps"))
      print(tbl, digits = digits, print.gap = print.gap, ...)
    } else {
      cat("\n\nVALIDATION: RMSEP")
      cat("\nCross-validated using", length(object$validation$segments),
          attr(object$validation$segments, "type"), "segments.\n")
      print(RMSEP(object), digits = digits, print.gap = print.gap, ...)
    }
  }
}

#' Explained variance (block-wise and component-wise)
#'
#' @param object A \code{rosa} object.
#' @param ncomp Integer to control the number of components to display (if fewer than the fitted number of components).
#' @param type Character indicating which type of explained variance to compute (default = "train", alternative = "CV").
#'
#' @return An object of class \code{rosaexpl} containing block-wise and component-wise explained variance contained in a matrix with attributes.
#' @export
#'
#' @examples
blockexpl <- function(object, ncomp = object$ncomp, type = c("train","CV")){
  nblock  <- length(object$X)
  X       <- object$X.orig
  categ   <- !is.null(object$classes)
  Y <- model.response(object$model)
  nobj    <- ifelse(is.null(dim(Y)), length(Y),nrow(Y))
  if(categ){
    Y <- model.matrix(~Y-1, data.frame(Y = factor(Y)))
  }
  Ytotvar <- sum(c(Y-rep(object$Ymeans, each = nobj))^2)
  Xtotvar <- object$Xtotvar
  Xvar    <- object$Xvar

  # Component-wise explained variance
  if(type[1] == "train"){
    fits <- object$fitted.values
  } else {
    fits <- object$validation$pred
    ncomp <- min(ncomp, dim(fits)[3])
  }
  resc <- matrix(0.0, 2, ncomp+1)
  resc[1,1] <- Xvar[1]/Xtotvar
  resc[2,1] <- 1-sum(c(fits[,,1]-Y)^2)/Ytotvar
  for(i in 2:ncomp){
    resc[1,i] <- Xvar[i]/Xtotvar
    resc[2,i] <- 1-sum(c(fits[,,i]-Y)^2)/Ytotvar - sum(resc[2,])
    # resc[2,i] <- sum(c(fits[,,i]-fits[,,i-1])^2)/Ytotvar
  }
  resc[,ncomp+1] <- sum(c(fits[,,i]-Y)^2)/Ytotvar #1-rowSums(resc) # res[,nblock+1]
  nameVec <- character(ncomp)
  for(i in 1:ncomp){
    nameVec[i] <- paste(names(object$X)[object$order[[i]]], collapse = ":")
  }
  dimnames(resc) <- list(c("X","Y"), c(paste(paste("comp.", 1:ncomp, sep=""), paste(' (',nameVec,')',sep=""),sep=""),"residual"))

  cblocks <- unlist(lapply(object$order, function(i)length(i)>1)) # Common blocks
  cnblock <- length(un <- unique(object$order[cblocks]))+nblock
  corder  <- as.list(c(1:nblock, un))

  # Block-wise explained variance
  res     <- matrix(0.0, 2, cnblock+1)
  for(i in 1:cnblock){
    ids <- unlist(lapply(object$order[1:ncomp], function(j)identical(j,corder[[i]])))
    if(sum(ids) > 0){
      res[1,i] <- sum(Xvar[ids]) / Xtotvar #- res[1,1]
      res[2,i] <- sum(resc[,-(ncomp+1)][2,ids])
    }
  }
  if(any(res[!is.na(res)] < 0)){
    warning('Negative block-wise explained variance encountered.')
  }
  res[,cnblock+1] <- 1-rowSums(res, na.rm = TRUE)
  nameVec <- character(cnblock)
  for(i in 1:cnblock){
    nameVec[i] <- paste(names(object$X)[corder[[i]]], collapse=":")
  }
  dimnames(res) <- list(c("X","Y"), c(nameVec,"residual"))
  if(length(remove <- nblock + which(colSums(res[,-c(1:nblock,cnblock+1), drop=FALSE])==0)) > 0){
    res <- res[,-remove]
    corder <- corder[-remove]
  }

  if(type[1] != "train"){
    res  <- res[2,,drop=FALSE]
    resc <- resc[2,,drop=FALSE]
  }

  attr(res,'compwise') <- resc
  attr(res,'index') <- corder
  class(res) <- c("rosaexpl","matrix")
  res
}

#' Print function for explained variance
#'
#' @param object A \code{rosaexpl} object.
#' @param digits Integer number of digits to print.
#' @param compwise Logical indicating if block-wise (default/FALSE) or component-wise (TRUE) explained variance should be printed.
#'
#' @return No return
#' @export
#'
#' @examples
print.rosaexpl <- function(object, digits = 3, compwise = FALSE){
  if(compwise){
    cat("Component-wise explained variance\n\n")
    print(round(attr(object, "compwise"), digits))
  } else {
    attr(object, "compwise") <- NULL
    attr(object, "index") <- NULL
    class(object) <- "matrix"
    cat("Block-wise explained variance\n\n")
    print(round(object, digits))
  }
}
#
# ## Print method for mvrVal objects:
# print.mvrVal <- function(x, digits = 4, print.gap = 2, ...) {
#   nresp <- dim(x$val)[2]
#   yvarnames <- dimnames(x$val)[[2]]
#   names(dimnames(x$val)) <- NULL
#   for (i in 1:nresp) {
#     if (nresp > 1) cat("\nResponse:", yvarnames[i], "\n")
#     print(x$val[,i,], digits = digits, print.gap = print.gap, ...)
#   }
#   invisible(x)
# }
