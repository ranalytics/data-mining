## todo:
## custom predict methods
## test many options
## implement more methods
## add/test dopar functionality


"fuse" <- function(mods, ...) UseMethod("fuse")

fuse.default <-
  function(mods, ## named list of models
           classes = NULL, ## class levels
           probs = TRUE, ## 
           predict = NULL,
           weights = rep(1, length(mods)),
           method = "vote",
           methodArgs = NULL,
           ...)
  {
    call <- match.call()
    if(!is.function(method))
      {
        if(length(method) > 1) stop("please speficy a single combination method")
        if(!(method %in% c("vote", "meanProb", "prod")))
          stop("method should be either a function of 'vote', 'meanProb', 'prod'")
      }
    if(!is.character(classes)) stop("classes should be a character vector")
    if(length(classes) < 2) stop("there should be at least two classes")
    
    if(!is.function(method))
      {
        if(probs & method %in% c("vote")) stop("this method is inconsistent with class probabilities")
      }

    if(is.null(names(mods))) names(mods) <- paste("model", seq(along = mods), sep = "")
    if(length(mods) == 1) stop("more than one model should be used")

    weights  <- weights/sum(weights)
    
    out <- list(models = mods,
                levels = classes,
                probs = probs,
                call = call,
                weights = weights,
                args = methodArgs,
                predict = predict,
                method = method)
    class(out) <- "fuse"
    out
  }

