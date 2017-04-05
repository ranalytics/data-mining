predict.fuse <-
function(object, newdata, ...)
  {
    library(foreach)

    if(is.function(object$method))
      {
        if(is.null(object$args))
          {
            object$args <- list(levels = object$levels)
          } else {
            object$args$levels <- object$levels
          }
      }
    
    if(!object$probs)
      {
        if(is.null(object$predict))
          {
            preds <- foreach(i = seq(along = object$models), .combine = cbind) %do% as.character(predict(object$models[[i]], newdata))
          } else {
            preds <- foreach(i = seq(along = object$models), .combine = cbind) %do% as.character(object$predict[[i]](object$models[[i]], newdata))
          }
        if(is.function(object$method))
          {          
            preds <- do.call(object$method, object$args)
            preds <- factor(as.character(preds), levels = object$levels)
          } else {
            if(object$method == "vote")
              {
                preds <- factor(apply(preds, 1, getMajVote), levels = object$levels)
              }
          }
      } else {
        probs <- foreach(i = seq(along = object$models), .combine = c) %do% predict(object$models[[i]], newdata, type = "prob")
        probs <- array(unlist(probs),
                       dim = c(
                         nrow(newdata),          # dim1 = rows are samples
                         length(object$levels),  # dim2 = cols are classes
                         length(object$models))) # dim3 = models
        if(is.function(object$method))
          {
            preds <- do.call(object$method, object$args)

          } else {
            if(object$method == "meanProb")
              {
                probs <- sweep(probs, 3L, object$weights, "*")
                probs <- apply(probs, 1:2, mean)
                preds <- factor(object$levels[apply(probs, 1, which.max)],
                                levels = object$levels)
              }          
            if(object$method == "prod")
              {
                probs <- apply(probs, 1:2, prod)
                preds <- factor(object$levels[apply(probs, 1, which.max)],
                                levels = object$levels)
              }
          }
      }
    preds
  }

