
n <- length(x$value)
x$value <- sort(x$value)
lower <- x$value[1 : (n - floor(n * level))]
upper <- x$value[(floor(n * level) + 1) : n]
ind <- which.min(upper - lower)
c(lower[ind], upper[ind])
