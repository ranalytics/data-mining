getMajVote <-
function(x)
  {
    x <- table(x)
    sample(names(x)[x == max(x)], 1)

  }

