print.fuse <-
function(x, ...)
  {
    cat("\nCall:\n", deparse(x$call), "\n\n", sep = "")
    meth <- if(is.function(x$method)) "custom algorithm" else x$method
    cat("Combination of", length(x$models), "classifiers using", meth, "\n")


  }

