explvar <- function(object)
  switch(class(object)[1],
         mvr = 100 * object$Xvar / object$Xtotvar,
         princomp =,
         prcomp = 100 * object$sdev^2 / sum(object$sdev^2),
         scores =,
         loadings = attr(object, "explvar"),
         rosa = 100 * object$Xvar / object$Xtotvar
  )

dummycode <- function(Y, n){
  nlev <- nlevels(Y)
  lev  <- levels(Y)
  X    <- model.matrix(~y-1,data.frame(y=Y))
  ref  <- X[,nlev,drop=FALSE]
  X    <- X[,-nlev,drop=FALSE]
  attributes(X) <- list(dim = attributes(X)$dim)
  X    <- X-ref%*%colSums(X)/sum(ref)
  X
}

model.matrix.rosa <- function (object, ...) {
  if (n_match <- match("X.concat", names(object), 0))
    object[[n_match]]
  else {
    data <- model.frame(object, ...)
    mm <- NextMethod("model.matrix", data = data)
    mm <- delete.intercept(mm)
    mt <- terms(object)
    if (length(attr(mt, "term.labels")) == 1 && !is.null(colnames(data[[attr(mt,
                                                                             "term.labels")]])))
      colnames(mm) <- sub(attr(mt, "term.labels"), "",
                          colnames(mm))
    return(mm)
  }
}
