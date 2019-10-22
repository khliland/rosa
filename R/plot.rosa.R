
#' Image plot for ROSA component interpretation.
#'
#' @param x A \code{rosa} object
#' @param type An optional \code{character} for selecting the plot type: "correlation" (default), "residual" or "order".
#' @param ncomp Integer to control the number of components to plot (if fewer than the fitted number of components).
#' @param col Colors used for the image plot, defaulting to mcolors(128).
#' @param legend Logical indicating if a legend should be included (default = TRUE).
#' @param mar Figure margins, default = c(5,6,4,7).
#' @param las Axis text direction, default = 1.
#' @param zlim A two element \code{numeric} limiting the color axis.
#' @param main The main title of the plot.
#' @param ... Additional parameters passed to \code{image}
#'
#' @return No return.
#' @importFrom graphics image
#' @importFrom plotrix color.legend
#' @export
#'
#' @examples
image.rosa <- function(x, type = c("correlation","residual","order"), ncomp = x$ncomp,
                       col = mcolors(128), legend = TRUE, mar = c(5,6,4,7), las = 1, zlim, main, ...){
  if(type[1] == "correlation"){
    im <- t(x$candidate.correlation)[1:ncomp,,drop=FALSE]
    if(missing(zlim))
      zlim <- c(-1,1)
    leg <- seq(-1,1,length.out = length(col))
  } else {
    if(type[1] == "residual"){
      im <- -t(x$candidate.RMSE)[1:ncomp,,drop=FALSE]
      if(missing(zlim))
        zlim <- c(min(im),max(im))
      leg <- seq(min(-im),max(-im),length.out = length(col))
    } else {
      if(type[1] == "order"){
        im <- t(x$candidate.RMSE)[1:ncomp,,drop=FALSE]
        for(i in 1:nrow(im)){
          im[i,] <- im[i,]-min(im[i,])
          im[i,] <- 1-im[i,]/max(im[i,])
          if(missing(zlim))
            zlim <- c(min(im),max(im))
        }
        leg <- seq(min(im),max(im),length.out = length(col))
      } else {
        stop("Unsupported plot type")
      }
    }
  }
  legw <- whichMins(legp <- pretty(leg),leg)
  legs <- character(length(col)); legs[legw[!is.na(legw)]] <- legp[!is.na(legw)]

  nresp <- ncol(im)
  if(missing(main))
    main <- ifelse(type[1] == "correlation", "Candidate score correlations", ifelse(type[1] == "residual","Candidate component RMSE","Candidate component residual order"))

  pars <- par(mar = mar, las = las)
  image(im, axes = FALSE, col = col, zlim = zlim, main = main, ...)
  axis(1, at = (0:(ncomp-1))/(ncomp-1), 1:ncomp)
  axis(2, at = (0:(nresp-1))/(nresp-1), colnames(im))
  box()
  for(i in 1:ncomp){
    w <- x$order[[i]]
    points(rep((i-1)/(ncomp-1), length(w)), (w-1)/(nresp-1), pch = 16, col = "white")
    points(rep((i-1)/(ncomp-1), length(w)), (w-1)/(nresp-1))
  }
  if(type[1] == "residual") col <- rev(col)
  color.legend(1.14,0,1.19,1, rect.col = col, gradient = TRUE, legend = legs)
  par(pars)
}


#' Color palette generation from matrix of RGB values
#'
#' @param n Integer number of colors to produce.
#' @param colmatrix A numeric \code{matrix} of three columns (R,G,B) to generate color palette from.
#'
#' @return A vector of n colors.
#' @export
#'
#' @examples
mcolors <- function(n, colmatrix = matrix(c(0,0,1, 1,1,1, 1,0,0), 3,3, byrow = TRUE)){
  cols <- character(n)
  colR <- approx(seq(0,1,length.out=nrow(colmatrix)), colmatrix[,1], (0:(n-1))/(n-1))$y
  colG <- approx(seq(0,1,length.out=nrow(colmatrix)), colmatrix[,2], (0:(n-1))/(n-1))$y
  colB <- approx(seq(0,1,length.out=nrow(colmatrix)), colmatrix[,3], (0:(n-1))/(n-1))$y
  for(i in 1:n){
    cols[i] <- rgb(colR[i],colG[i],colB[i])
  }
  cols
}

whichMins <- function(short, long){
  out <- integer(length(short))
  for(i in 1:length(short)){
    out[i] <- which.min((short[i]-long)^2)
  }
  out[short < min(long)] <- NA; out[short > max(long)] <- NA
  out
}

#' Barplot of block and component explained variances
#'
#' @param height A \code{rosa} object.
#' @param type Character indicating if explained variance should be based on training data ("train") or cross-validation ("CV").
#' @param ncomp Integer to control the number of components to plot (if fewer than the fitted number of components).
#' @param col Colors to use for the components.
#' @param ... Additional parameters passed to \code{barplot}
#'
#' @return No return
#' @export
#'
#' @examples
barplot.rosa <- function(height, type = c("train","CV"), ncomp = height$ncomp, col = mcolors(ncomp),  ...){
  nums <- attr(numsb <- blockexpl(height, type = type[1], ncomp = ncomp), "compwise")
  nums <- nums[nrow(nums),,drop=FALSE]
  maxCount <- max(height$count)
  nam <- colnames(numsb)[-ncol(numsb)]
  mat <- matrix(0, ncomp, length(nam))
  counts <- rep(1, length(nam))
  corder <- attr(numsb, 'index')
  for(i in 1:ncomp){
    ids <- unlist(lapply(corder, function(j)identical(j,height$order[[i]])))
    mat[i, ids] <- nums[i]
    counts[ids] <- counts[ids] + 1
  }
  colnames(mat) <- nam
  barplot(mat, col = col, ...)
}


#' Loading weights plot
#'
#' @param object A \code{rosa} object.
#' @param ylab Character label to show on the y axis (default = "loading weight")
#' @param ... Additional parameters passed to \code{loadingplot}
#'
#' @return
#' @export
#'
#' @examples
loadingweightplot <- function(object, ylab = "loading weight", ...){
  object$loadings <- object$loading.weights
  loadingplot(object, ylab = ylab, ...)
}
