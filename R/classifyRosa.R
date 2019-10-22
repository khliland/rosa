#' @importFrom MASS lda qda
classifyRosa <- function(model, grouping, newdata, ncomp, LQ){
  if(LQ == "max"){
    labels  <- names(table(grouping))
    predVal <- predict(model, newdata = newdata, ncomp = 1:ncomp)
    class   <- apply(predVal,c(1,3),which.max)
    for(i in 1:ncol(class)){
      class[[i]]   <- labels[class[[i]]]
    }
    colnames(class) <- paste("Comp.", 1:ncomp, sep="")
    return(class)

  } else { # LDA or QDA
    # Extract and predict scores
    scoresCal <- scores(model)
    scoresVal <- predict(model, newdata = newdata, type = "scores")

    # Prepare for storage
    N <- dim(scoresVal)
    class <- matrix(0, N[1],ncomp)

    # Create ncomp lda models and predict classes
    for(i in 1:ncomp){
      if(LQ == "lda"){
        ldai <- lda(scoresCal[, 1:i, drop = FALSE], grouping, tol = 1.0e-10)
      }
      if(LQ == "qda"){
        ldai <- qda(scoresCal[, 1:i, drop = FALSE], grouping, tol = 1.0e-10)
      }
      class[, i] <- predict(ldai, scoresVal[, 1:i, drop = FALSE])$class
    }
    colnames(class) <- paste("Comp.", 1:ncomp, sep="")
    return(class)
  }
}
