surf.pco <- function (ord, var, ax = 1, ay = 2, thinplate = TRUE, col = 2, 
    labcex = 0.8, family = gaussian, grid = 50, gamma = 1, ...) 
{
    if (class(ord) != "pco") 
        stop("You must supply an object of class pco from pco()")
    if (missing(var)) {
        stop("You must specify a variable to surface")
    }
    x <- ord$points[, ax]
    y <- ord$points[, ay]
    if (any(is.na(var))) {
        cat("Omitting plots with missing values \n")
        x <- x[!is.na(var)]
        y <- y[!is.na(var)]
        var <- var[!is.na(var)]
    }
    if (is.logical(var)) {
        if (thinplate) 
            tmp <- gam(var ~ s(x, y), gamma = gamma, family = binomial)
        else tmp <- gam(var ~ s(x) + s(y), family = binomial, 
            gamma = gamma)
    }
    else {
        if (thinplate) 
            tmp <- gam(var ~ s(x, y), gamma = gamma, family = family)
        else tmp <- gam(var ~ s(x) + s(y), family = family, gamma = gamma)
    }
    new.x <- seq(min(x), max(x), len = grid)
    new.y <- seq(min(y), max(y), len = grid)
    xy.hull <- chull(x, y)
    xy.hull <- c(xy.hull, xy.hull[1])
    new.xy <- expand.grid(x = new.x, y = new.y)
    inside <- as.logical(pip(new.xy$x, new.xy$y, x[xy.hull], 
        y[xy.hull]))
    fit <- predict(tmp, type = "response", newdata = as.data.frame(new.xy))
    fit[!inside] <- NA
    contour(x = new.x, y = new.y, z = matrix(fit, nrow = grid), 
        add = TRUE, col = col)
    print(tmp)
    d2 <- (tmp$null.deviance - tmp$deviance)/tmp$null.deviance
    cat(paste("D^2 = ", formatC(d2, width = 4), "\n"))
}
