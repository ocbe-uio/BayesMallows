L <- 10
n <- 100
nmc <- 1e5
type1 <- replicate(nmc, {
  l <- sample(L, 1)
  u <- sample(n - l, 1)
  u
})

type2 <- replicate(nmc, {
  sample(n - L, 1)
})

hist(type1)
hist(type2)
