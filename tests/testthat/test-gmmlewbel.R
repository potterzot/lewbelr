
library(lfe)
library(gmm)

# Generate fake data
set.seed(1234)
N <- 1000
b1 <- c(1,1,1)
b2 <- c(1,1,1)
b <- 2
g <- 1
d <- 2
true_coefs <- c(b2, g, b1)
d2012 <- data.frame(
  x1 = rnorm(N, 0, 1),
  x2 = rnorm(N, 0, 1),
  u = rnorm(N, 0, 1),
  s1 = rnorm(N, 0, 1),
  s2 = rnorm(N, 0, 1),
  ov = rnorm(N, 0, 1),
  z1 = rnorm(N, 0, 1)
) %>%
  dplyr::mutate(
    e1 = u + exp(x1)*s1 + exp(x2)*s1,
    e2 = u + exp(-x1)*s2 + exp(-x2)*s2,
    y2 = b2[1] + b2[2]*x1 + b2[3]*x2 + ov + e2 + d*z1,
    y1 = b1[1] + b1[2]*x1 + b1[3]*x2 + g*y2 + b*ov + e1)


test_that("gmmlewbel() for heteroskedastic errors has correct output", {
  # Use engeldat for testing
  url <- "http://fmwww.bc.edu/ec-p/data/wooldridge/engeldat.dta"
  tmp <- tempfile("engeldat")
  download.file(url, tmp)
  engel <- haven::read_dta(tmp)
  unlink(tmp)

  # First estimate the OLS of the first stage with no constructed instruments


  ## Match to Lewbel 2012 Table 1:
  f_x <- formula(~age + age2 + agesp + agesp2 + spwork + s1 + s2 + s3 + washer + gasheat + onecar + twocars)
  str_x <- paste(all.vars(f_x), collapse = " + ")
  X <- as.matrix(cbind(1, engel[, all.vars(f_x)]))
  X2 <- X
  Z2 <- X2

  # Column 1
  m_col1 <- lm(update(f_x, foodshare ~ lrtotexp + .), data = engel)
  expect_true(round(coef(m_col1)[2], 3) + 0.127 == 0)
  expect_true(round(mean(X %*% coef(m_col1)[-2]), 3) == 0.361)
  expect_true(round(sd(colMeans(X) * coef(m_col1)[-2]), 4) == 0.0056)   # not working

  # Column 2
  f_col2 <- formula(foodshare ~ age + age2 + agesp + agesp2 + spwork + s1 + s2 + s3 + washer + gasheat + onecar + twocars | 0 | (lrtotexp ~ lrinc))
  m_col2 <- lfe::felm(f_col2, data = engel)
  expect_true(round(coef(m_col2)[14], 3) + 0.086 == 0)
  expect_true(round(mean(X %*% coef(m_col2)[-14]), 3) == 0.336)
  expect_true(round(sd(colMeans(X) %*% coef(m_col2)[-14]), 3) == 0.012)   # not working

  # Column 3
  # We construct instruments Z = X[-1] - colMeans(X[-1])) * e2
  e2 <- residuals(lm(update(f_x, lrtotexp ~ .), data = engel))
  Z1 <- scale(X[,-1], scale = FALSE) * e2
  colnames(Z1) <- paste0("z_", colnames(Z1))
  str_z <- paste(colnames(Z1), collapse = " + ")
  engelz <- cbind(engel, Z1)

  f_col3 <- formula(paste0("foodshare ~ age + age2 + agesp + agesp2 + spwork + s1 + s2 + s3 + washer + gasheat + onecar + twocars | 0 | (lrtotexp ~ ",
                           str_z, ")"))
  m_col3 <- lfe::felm(f_col3, data = engelz)
  expect_true(round(coef(m_col3)[14], 3) + 0.055 == 0)
  expect_true(round(mean(X %*% coef(m_col3)[-14]), 3) == 0.318)
  expect_true(round(sd(colMeans(X) * coef(m_col3)[-14]), 3) == 0.035)   # not working

  # Col 3 using REndo package
  library(REndo)
  f_endo5 <- formula(paste0("foodshare ~ age + age2 + agesp + agesp2 + spwork + s1 + s2 + s3 + washer + gasheat + onecar + twocars + lrtotexp | lrtotexp | IIV(",
                            paste(all.vars(f_x), collapse = ", "),
                            ")"))
  m_endo5 <- hetErrorsIV(f_endo5, data = engelz)


  # Column 4
  f_col4 <- formula(paste0("foodshare ~ age + age2 + agesp + agesp2 + spwork + s1 + s2 + s3 + washer + gasheat + onecar + twocars | 0 | (lrtotexp ~ lrinc + ",
                           str_z, ")"))
  m_col4 <- lfe::felm(f_col4, data = engelz)
  expect_true(round(coef(m_col4)[14], 4) + 0.0846 == 0) # This should be 0.086 in the paper
  expect_true(round(mean(X %*% coef(m_col4)[-14]), 3) == 0.336)
  expect_true(round(sd(colMeans(X) %*% coef(m_col4)[-14]), 3) == 0.011)   # not working

  ### Column 5
  # Now we get into GMM
  init <- coef(m_col2)
  f_gmm1 <- update(f_x, foodshare ~ . + lrtotexp)
  f_z5 <- formula(paste0("~ lrinc + ", str_x, ""))

  # Straight twostep gmm
  m_col5a <- gmm(f_gmm1, f_z5, data = engelz,
                vcov = "iid",
                t0 = init)
  coef(m_col5a)

  # Efficient gmm as described in stata doc
  Z <- as.matrix(cbind(engel[, "lrinc"], Z2))
  u <- residuals(m_col2)
  S <- t(Z * diag(u^2)) %*% Z / nrow(Z)
  W <- chol2inv(chol(S))

  m_col5b <- gmm(f_gmm1, f_z5, data = engelz,
                 weightsMatrix = W,
                 t0 = init,
                 vcov = "iid")
  coef(m_col5b)
  expect_true(round(coef(m_col5)[14], 3) + 0.086 == 0)
  expect_true(round(mean(X %*% coef(m_col5)[-14]), 3) == 0.336)
  expect_true(round(sd(colMeans(X) %*% coef(m_col5)[-14]), 3) == 0.012)   # not working

  # Using the new "modelfit" package
  mod5 <- momentModel(update(f_x, foodshare ~ . + lrtotexp), f_z5,
                      data = engelz, vcov = "iid")
  w0 <- evalWeights(mod5, w = "ident")
  m_col5a <- gmmFit(mod5, type = "onestep", t0 = init)
  w1 <- evalWeights(mod5, solveGmm(mod5, w0)$theta)
  m_col5b <- gmmFit(mod5, weights = w1, t0 = init)
  coef(m_col5a)
  coef(m_col5b)

  ### Column 6
  # ref: https://ageconsearch.umn.edu/record/116029/files/sjart_st0030.pdf
  # First we create the weighting matrix based on the IV estimation
  u <- residuals(m_col3)
  S <- t(Z * diag(u^2)) %*% Z / nrow(Z)
  W <- chol2inv(chol(S))

  f_z6 <- formula(paste0("~ ", str_z, " + ", str_x, ""))
  m_col6 <- gmm(update(f_x, foodshare ~ . + lrtotexp), f_z6, data = engelz,
                vcov = "iid",
                type = "iterative",
                t0 = coef(m_col3))
  m_col6

  m_col6 <- gmm(update(f_x, foodshare ~ . + lrtotexp), f_z6, data = engelz,
                t0 = coef(m_col3),
                weightsMatrix = W)
  m_col6

  expect_true(round(coef(m_col6)[14], 3) + 0.055 == 0) # This should be 0.078 in the paper
  expect_true(round(mean(X %*% coef(m_col6)[-14]), 3) == 0.318) # Should be 0.332 in paper
  expect_true(round(sd(colMeans(X) %*% coef(m_col6)[-14]), 3) == 0.028)   # not working

  mod6 <- momentModel(update(f_x, foodshare ~ . + lrtotexp), f_z6, data = engelz, vcov = "iid")
  m_col6a <- gmmFit(mod6, type = "onestep")
  w0 <- evalWeights(mod6, w = "ident")
  w1 <- evalWeights(mod6, solveGmm(mod6, w0)$theta)
  m_col6b <- gmmFit(mod6, weights = w1)
  m_col6c <- gmmFit(mod6, type = "twostep")

  coef(m_col6a)
  coef(m_col6b)
  coef(m_col6c)



  # Column 7
  f_z7 <- formula(paste0("~ lrinc + ", str_z, " + ", str_x))
  m_col7 <- gmm(update(f_x, foodshare ~ . + lrtotexp), f_z7, data = engelz, vcov = "iid")
  expect_true(round(coef(m_col7)[14], 4) + 0.087 == 0) # This should be 0.087 in the paper
  expect_true(round(mean(X %*% coef(m_col7)[-14]), 3) == 0.337)
  expect_true(round(sd(colMeans(X) %*% coef(m_col7)[-14]), 3) == 0.011)   # not working




  # Column 5 using our own moments
  engelz2 <- cbind(engelz[, "foodshare"],
                   engelz[, "lrtotexp"],
                   1,
                   engelz[, c(all.vars(f_x), colnames(Z))])

  g5 <- function(theta, df) {
    k <- 12
    Y1 <- data.matrix(df[, 1])
    Y2 <- data.matrix(df[, 2])
    X <- data.matrix(df[, 2 + 1:(k+1)])
    Z <- data.matrix(df[, 2 + (k+1) + 1:k])
    beta1 <- theta[1:(k+1)]
    gamma <- theta[k+1 + 1]
    beta2 <- theta[1 + (k+1) + 1:(k+1)]
    mu <- theta[1 + 2*(k+1) + 1:k]

    # moments
    e1 <- Y1 - X %*% beta1 - Y2 %*% gamma
    e2 <- Y2 - X %*% beta2

    m1 <- X * as.vector(e1)
    m2 <- X * as.vector(e2)
    m3 <- Z - mu
    m4 <- (Z - mu) * e2[,1] * e1[,1]
    M <- cbind(m1, m2, m3, m4)
    M
  }

  # First using the identity matrix
  init <- c(coef(m_col2), coef(lm(update(f_x, lrtotexp ~ .), data = engelz)), colMeans(Z))
  k <- length(init)
  w0 <- diag(1, 50, 50)
  m1 <- g5(init, engelz2)
  v1 <- crossprod(m1)
  w1 <- solve(v1)
  chol(w1)
  w1 <- chol2inv(chol(S))

  m_col5z <- gmm(g5, engelz2, t0 = init, vcov = "iid", weightsMatrix = w1)
      vcov = "TrueFixed",
      weightsMatrix = w0)
  coef(m_col5z)


  m_col5a <- gmm(g0, engelz2, vcov = "iid", t0 = rep(0, 39))
  expect_true(round(coef(m_col6)[14], 3) + 0.055 == 0) # This should be 0.078 in the paper
  expect_true(round(mean(X %*% coef(m_col6)[-14]), 3) == 0.318) # Should be 0.332 in paper
  expect_true(round(sd(colMeans(X) %*% coef(m_col6)[-14]), 3) == 0.028)   # not working






      resid1 <- residuals(ols1)

  # Constructed instruments
  d2012a <- d2012 %>%
    dplyr::mutate(
      hetz1 = resid1 * (x1 - mean(x1)),
      hetz2 = resid1 * (x2 - mean(x2))
    )

  # 2SLS
  tsls1 <- lm(y2 ~ x1 + x2 + hetz1 + hetz2 + z1, data = d2012a)
  d2012a$y2hat <- fitted(tsls1)
  tsls2 <- lm(y1 ~ x1 + x2 + y2hat, data = d2012a)

  tsls <- lfe::felm(y1 ~ 1 + x1 + x2 | 0 | (y2 ~ hetz1 + hetz2 + z1) | 0, data = d2012a)

  # GMM
  gmm1 = gmm(y1 ~ y2 + x1 + x2, ~ x1 + x2 + hetz1 + hetz2 + z1, data = d2012a,
             vcov="iid")

  # AER
  aer <- AER::ivreg(y1 ~ x1 + x2 + y2 | x1 + x2 + hetz1 + hetz2 + z1, data = d2012a)

  # REndo package
  library(REndo)
  data("dataHetIV")
  het <- hetErrorsIV(y1 ~ x1 + x2 + y2 | y2 | IIV(x1, x2) | z1, data = d2012)

  # These should match, but aren't the real test
  expect_true(all(abs(coef(tsls2) - coef(tsls)) < 1e-10))
  expect_true(all(abs(coef(tsls) - coef(gmm1)[c(1,3,4,2)]) < 1e-10))

  # Now run gmmlewbel to compare
  init <- c(coef(tsls), coef(ols1), rep(1, 4))
  names(init)[length(init) - 3:0] <- c("t_beta", "t_U", "t_V", "t_R")
  formula <- y1 ~ y2 | x1 + x2 | x1 + x2 + z1
  m <- gmmlewbel(formula = formula, data = d2012,
                 error_type = "hetero", init = init, type = "iterative")
  expect_true(all(abs(coef(tsls) - coef(m)[1:4]) < 1e-2))
})

test_that("gmmlewbel() for non-normal errors has correct output with no instruments and no covariates", {
  data("AcemogluJohnson2007")
  ols1 <- lm(loglifeexpect ~ compsjmhatit, data = AcemogluJohnson2007)
  AcemogluJohnson2007$loglifeexpect_hat <- fitted(ols1)
  ols2 <- lm(loggdppcmadd ~ loglifeexpect_hat, data = AcemogluJohnson2007)

  # Exponential variances
  init <- c(coef(ols2)[2], log(c(2.7, 0.0248, 0.006, 0.08)))
  names(init) <- c("loglifeexpect", "t_beta", "t_s2U", "t_s2V", "t_s2R")
  formula <- loggdppcmadd ~ loglifeexpect
  #init <- init[c(1,3)]
  m <- gmmlewbel(formula = formula,
                 data = AcemogluJohnson2007,
                 init = init,
                 error_type = "nonnormal",
                 exp_variances = TRUE,
                 use_overid = FALSE,
                 type = "iterative",
                 vcov = "iid",
                 itermax = 1000,
                 #method = "L-BFGS-B",
                 control=list(fnscale=1e-12))
  coefs <- c(coef(m)[1], exp(coef(m)[-1]))
  coefs
  summary(m)
  coefs - c(-3.046, 2.791, 0.0248, 0.00587, 0.0838)
  expect_true(all(abs(coef(m) - c(-3.046, 2.791, 0.0248, 0.00587, 0.0838)) < 1e-10))

  ## Overidentified
  init <- c(init, log(0.146))
  names(init)[6] <- c("t_uWW")

  # Get the moment
  mfn <- moments_nonnormal(0, 0, TRUE, TRUE)
  g0 <- mfn(init, AcemogluJohnson2007[, c("loggdppcmadd", "loglifeexpect")])
  g0 <- scale(g0, scale = FALSE)
  class(g0) <- "gmmFct"
  v0 <- sandwich::vcovHAC(g0, sandwich = FALSE)
  w0 <- solve(v0)

  q <- qr(g0/sqrt(nrow(g0)))
  w1 <- matrix(NA, ncol(g0), ncol(g0))
  w1[q$pivot, q$pivot] <- chol2inv(q$qr)
  gmm(mfn, AcemogluJohnson2007[, c("loggdppcmadd", "loglifeexpect")], t0 = init,
      type = "iterative",
      vcov = "iid",
      weightsMatrix = w1,
      method = "BFGS",
      control = list(fnscale = 1e-8))

  m <- gmmlewbel(formula = formula,
                 data = AcemogluJohnson2007,
                 init = init,
                 error_type = "nonnormal",
                 use_overid = TRUE,
                 exp_variances = TRUE,
                 type = "iterative",
                 vcov = "TrueFixed",
                 weightsMatrix = w1,
                 method = "BFGS",
                 control = list(fnscale = 1e-8))



  m <- gmmlewbel(formula = formula,
                 data = AcemogluJohnson2007,
                 init = init,
                 error_type = "nonnormal",
                 use_overid = TRUE,
                 exp_variances = TRUE,
                 type = "iterative",
                 #method = "BFGS",
                 control = list(fnscale = 1e-8))
  coef(m)
  expect_true(all(abs(coef(m) - 1) < 1e-10))

  # Variances not exponentiated
  init <- c(coef(ols2)[2], c(c(2.7, 0.0248, 0.006, 0.08)))
  names(init) <- c("loglifeexpect", "beta", "s2U", "s2V", "s2R")
  formula <- loggdppcmadd ~ loglifeexpect
  m <- gmmlewbel(formula = formula,
                 data = AcemogluJohnson2007,
                 init = init,
                 error_type = "nonnormal",
                 exp_variances = FALSE,
                 use_overid = FALSE,
                 type = "iterative",
                 method = "L-BFGS-B",
                 control=list(fnscale=1e-8),
                 lower = c(-Inf, rep(1e-6, 4)),
                 upper = c(rep(Inf, 5)))
  coef(m)
  summary(m)

  # Overidentified
  init <- c(init, 0.146)
  names(init)[6] <- c("uWW")
  m <- gmmlewbel(formula = formula,
                 data = AcemogluJohnson2007,
                 init = init,
                 error_type = "nonnormal",
                 use_overid = TRUE,
                 exp_variances = FALSE,
                 type = "iterative",
                 method = "L-BFGS-B",
                 control=list(fnscale=1e-8),
                 lower = c(-Inf, rep(1e-5, 5)),
                 upper = c(rep(Inf, 6)))
  coef(m)
  expect_true(all(abs(coef(m) - 1) < 1e-10))
})

test_that("gmmlewbel() for non-normal errors has correct output with covariates but no exogenous instruments", {
  # With covariates
  ols <- lm(loglifeexpect ~ loggdppcmadd301980, data = AcemogluJohnson2007)
  tsls <- lfe::felm(loggdppcmadd ~ loglifeexpect + loggdppcmadd301980 | 0 | 0 | 0,
                    data = AcemogluJohnson2007)
  init <- c(coef(tsls)[2], coef(tsls)[c(1,3)], coef(ols), rep(1, 4))
  names(init)[length(init) - 3:0] <- c("Beta", "s2U", "s2V", "s2R")
  formula <- loggdppcmadd ~ loglifeexpect | loggdppcmadd301980
  m <- gmmlewbel(formula = formula,
                 data = AcemogluJohnson2007,
                 init = init,
                 error_type = "nonnormal",
                 type = "iterative",
                 method = "L-BFGS-B",
                 control=list(fnscale=1e-8),
                 lower = c(rep(-Inf,5), rep(1e-5, 4)),
                 upper = c(rep(Inf, 9)))

  coef(m)
  expect_true(all(abs(coef(m) - c(-0.378, -0.139, 0.150, 1.758, -0.184, 2.090, 0, 0.014, 0.122)) < 1e-3))

  # With an instrument and no covariates
  ols <- lm(loglifeexpect ~ compsjmhatit - 1, data = AcemogluJohnson2007)
  tsls <- lfe::felm(loggdppcmadd ~ -1 | 0 | (loglifeexpect ~ compsjmhatit),
                    data = AcemogluJohnson2007)

  init <- c(coef(tsls)[1], 0.65, log(c(2.147, 0.00836, 0.0217, 0.0846)
  names(init)[length(init) - 4:0] <- c("Delta", "Beta", "s2U", "s2V", "s2R")
  formula <- loggdppcmadd ~ loglifeexpect | 0 | compsjmhatit
  m <- gmmlewbel(formula = formula,
                 data = AcemogluJohnson2007,
                 init = init,
                 error_type = "nonnormal",
                 type = "iterative",
                 method = "L-BFGS-B",
                 control=list(fnscale=1e-8),
                 lower = c(rep(-Inf,2), rep(1e-5, 4)),
                 upper = c(rep(Inf, 6)))
  coef(m)
  expect_true(all(abs(coef(m) - 1) < 1e-10))

  init <- c(init, 1)
  names(init)[6] <- "uWW"
  m <- gmmlewbel(formula = formula,
                 data = AcemogluJohnson2007,
                 init = init,
                 error_type = "nonnormal",
                 use_overid = TRUE,
                 exp_variances = FALSE,
                 type = "iterative",
                 method = "L-BFGS-B",
                 control=list(fnscale=1e-8),
                 lower = c(rep(-Inf,1), rep(1e-5, 4)),
                 upper = c(rep(Inf, 9)))
  coef(m)
  expect_true(all(abs(coef(m) - 1) < 1e-10))


  # With an instrument and covariates
  ols <- lm(loglifeexpect ~ loggdppcmadd301980, data = AcemogluJohnson2007)
  tsls <- lfe::felm(loggdppcmadd ~ loggdppcmadd301980 | 0 | (loglifeexpect ~ compsjmhatit),
                    data = AcemogluJohnson2007)

  init <- c(coef(tsls)[3], coef(tsls)[c(1,2)], coef(ols), coef(tsls$stage1)[3], rep(1, 4))
  names(init)[length(init) - 3:0] <- c("Beta", "s2U", "s2V", "s2R")
  formula <- loggdppcmadd ~ loglifeexpect | loggdppcmadd301980 | compsjmhatit
  m <- gmmlewbel(formula = formula,
                 data = AcemogluJohnson2007,
                 init = init,
                 error_type = "nonnormal",
                 type = "iterative",
                 method = "L-BFGS-B",
                 control=list(fnscale=1e-8),
                 lower = c(rep(-Inf,5), rep(1e-5, 4)),
                 upper = c(rep(Inf, 9)))

  coef(m)
  expect_true(all(abs(coef(m) - 1) < 1e-10))

  init <- c(init, 1)
  m <- gmmlewbel(formula = formula,
                 data = AcemogluJohnson2007,
                 init = init,
                 error_type = "nonnormal",
                 use_overid = TRUE,
                 type = "iterative",
                 method = "L-BFGS-B",
                 control=list(fnscale=1e-8),
                 lower = c(rep(-Inf,5), rep(1e-5, 5)),
                 upper = c(rep(Inf, 10)))

  coef(m)
  expect_true(all(abs(coef(m) - 1) < 1e-10))
})


test_that("gmmlewbel() for non-normal errors has correct output with covariates", {
  data("AcemogluJohnson2007")
})

