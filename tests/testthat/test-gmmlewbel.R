set.seed(510)
N <- 100
b1 <- c(1,1,1)
b2 <- c(1,1,1)
b <- 2
g <- 0.5
true_coefs <- c(b2, g, b1)
d2012 <- data.frame(
  x1 = rep(1, N),
  x2 = rnorm(N, 0, 1),
  x3 = rnorm(N, 0, 1),
  u = rnorm(N, 0, 1),
  s1 = rnorm(N, 0, 1),
  s2 = rnorm(N, 0, 1),
  ov = rnorm(N, 0, 1)
) %>%
  dplyr::mutate(e1 = u + exp(x2)*s1 + exp(x3)*s1,
         e2 = u + exp(-x2)*s2 + exp(-x3)*s2,
         y2 = b2[1]*x1 + b2[2]*x2 + b2[3]*x3 + ov + e2,
         y1 = b1[1]*x1 + b1[2]*x2 + b1[3]*x3 + g*y2 + b*ov + e1)

test_that("gmmlewbel() for heteroskedastic errors has correct output", {
  f <- y1 | y2 ~ x2 + x3 + y2 | x2 + x3 | x2 + x3
  m <- gmmlewbel(formula = f, data = d2012)
  expect_true((coef(m) - rep(0, 7))^2 < 1e-10)
})

test_that("gmmlewbel() for non-normal errors has correct output", {
  data("AcemogluJohnson2007")
  f <- loggdppcmadd | loglifeexpect ~ loggdppcmadd301980 + loglifeexpect | loggdppcmadd301980 | compsjmhatit
  m <- gmmlewbel(formula = f, data = AcemogluJohnson2007)
  expect_true((coef(m)[3] - -0.378)^2 < 1e-10)
})
