library(gmm)
library(lewbelr)
library(lmtest)

data("AcemogluJohnson2007")
AcemogluJohnson2007$constant <- 1

# Simple regression with specified moments
f1 <- loggdppcmadd ~ loglifeexpect
ols <- lm(f1, data = AcemogluJohnson2007, y = TRUE, x = TRUE)
dt <- cbind(AcemogluJohnson2007[, "loggdppcmadd"], model.matrix(f1, data = AcemogluJohnson2007))
g <- function(theta, dt) {
  Y <- data.matrix(dt[,1])
  X <- data.matrix(dt[, 2:3])
  b <- theta
  e <- Y - X %*% b
  m1 <- e[,1] * X
  M <- cbind(m1)
  return(M)
}
theta <- c(b1 = 0, b2 = 0)
m <- gmm(g, dt, t0 = theta)
cbind(coeftest(m)[,1:2], coeftest(ols)[,1:2])


# Simple regression with an instrument and a covariate
ols1 <- lm(loglifeexpect ~ compsjmhatit + loggdppcmadd301980,
           data = AcemogluJohnson2007, y = TRUE, x = TRUE)
AcemogluJohnson2007$loglifeexpect_hat <- fitted(ols1)
ols2 <- lm(loggdppcmadd ~ loglifeexpect_hat + loggdppcmadd301980,
          data = AcemogluJohnson2007, y = TRUE, x = TRUE)

g <- function(theta, dt) {
  Y <- data.matrix(dt[,1])
  X <- data.matrix(dt[, 2:4])
  Z <- data.matrix(dt[, c(2,3,5)])
  b <- theta
  e <- Y - X %*% b
  m1 <- e[,1] * Z
  M <- cbind(m1)
  return(M)
}
dt <- AcemogluJohnson2007[, c("loggdppcmadd", "constant", "loggdppcmadd301980", "loglifeexpect", "compsjmhatit")]
gmm(g, dt, t0 = c(b1 = 0, b2 = 0, b3 = 0),
    method = "BFGS", control=list(fnscale=1e-8), type = "iterative")
gmm(loggdppcmadd ~ loggdppcmadd301980 + loglifeexpect, ~ loggdppcmadd301980 + compsjmhatit ,
    data = dt, t0 = c(b1 = 0, b2 = 0, b3 = 0),
    method = "BFGS", control=list(fnscale=1e-8), type = "iterative")
coef(ols2)

##################
# Using lewbel2020 with no covariates and demeaning / using just an intercept

g <- function(theta, dt) {
  W <- data.matrix(dt[,1])
  Y <- data.matrix(dt[,2])
  gamma <- theta[1]

  idx <- 1
  B <- theta[idx+ 1]
  s2U <- theta[idx + 2]
  s2V <- theta[idx + 3]
  s2R <- theta[idx + 4]
  #uWW <- theta[idx + 5]

  uYY <- s2U + s2V
  uYW <- B * s2U + gamma * uYY

  a <- gamma + B
  Q <- W - Y * gamma
  P <- W - Y * a
  Y2 <- Y^2


  # Moments. See Appendix 2 of Lewbel et al. (2020), eq. 74-77.
  m1 <- Y * W - uYW
  m2 <- Y2 - uYY
  m3 <- Q * P * Y
  m4 <- Q * P * (Y2 - uYY) - 2 * B * s2U * P * Y
  m5 <- Q^2 - B^2 * s2U - s2R
  M <- cbind(m1, m2, m3, m4, m5)
  M
}
dt <- scale(AcemogluJohnson2007[, c("loggdppcmadd", "loglifeexpect")], scale = FALSE)
theta <- c(gamma = -3, beta = 2, s2U = 0.025, s2V = 0.005, s2R = 0.01)
m <- gmm(g, dt, t0 = theta,
         control=list(
           fnscale=1e-4
           ),
         type = "iterative", method = "L-BFGS-B",
         kernel = "Bartlett",
         vcov = "iid",
         itermax = 500,
         lower = c(-Inf, rep(1e-3, 4)),
         upper = rep(Inf, 5))
m

g0 <- g(theta, dt)
g0 <- scale(g0, scale = FALSE)
class(g0) <- "gmmFct"
v0 <- vcovHAC(g0, sandwich = FALSE)
w0 <- corpcor::pseudoinverse(v0)
m <- gmm(g, dt, t0 = theta, weightsMatrix = w0,
         control=list(
           fnscale=1e-4
         ),
         type = "iterative", method = "L-BFGS-B",
         kernel = "Bartlett",
         itermax = 500,
         vcov  = "TrueFixed")
m

## Overidentified
g <- function(theta, dt) {
  W <- data.matrix(dt[,1])
  Y <- data.matrix(dt[,2])
  gamma <- theta[1]

  idx <- 1
  B <- theta[idx+ 1]
  s2U <- theta[idx + 2]
  s2V <- theta[idx + 3]
  s2R <- theta[idx + 4]
  uWW <- theta[idx + 5]

  uYY <- s2U + s2V
  uYW <- B * s2U + gamma * uYY

  a <- gamma + B
  Q <- W - Y * gamma
  P <- W - Y * a
  Y2 <- Y^2
  W2 <- W^2
  Y3 <- Y^3
  Y5a <- Y^5 - 10 * Y3 * uYY

  # Moments. See Appendix 2 of Lewbel et al. (2020), eq. 74-77.
  m1 <- Y * W - uYW
  m2 <- Y2 - uYY
  m3 <- Q * P * Y
  m4 <- Q * P * (Y2 - uYY) - 2 * B * s2U * P * Y
  m5 <- Q^2 - B^2 * s2U - s2R
  m6 <- W2 - uWW
  m7 <- -3 * Y * W2 * uYY - 6 * uYW * Y2 * W - uWW * Y3 +
    Y3 * W2 - a^2 * Y5a -
    (gamma + a) * (-6 * uYY * Y2 * W - 4 * uYW * Y3 + Y^4 * W - a * Y5a)
  M <- cbind(m1, m2, m3, m4, m5, m6,m7)
  M
}
dt <- scale(AcemogluJohnson2007[, c("loggdppcmadd", "loglifeexpect")], scale = FALSE)
theta <- c(gamma = -3, beta = 2, s2U = 0.025, s2V = 0.005, s2R = 0.01, uWW = 0.14)
m <- gmm(g, dt, t0 = theta,
    control=list(fnscale=1e-8),
    type = "twoStep",
    method = "L-BFGS-B",
    vcov = "iid",
    lower = c(-Inf, rep(1e-3, 5)),
    upper = rep(Inf, 6))

g0 <- g(theta, dt)
g0 <- scale(g0, scale = FALSE)
class(g0) <- "gmmFct"
v0 <- vcovHAC(g0, sandwich = FALSE)
w0 <- corpcor::pseudoinverse(v0)
m <- gmm(g, dt, t0 = theta, weightsMatrix = w0,
         control=list(
           fnscale=1e-4
         ),
         type = "iterative", method = "L-BFGS-B",
         kernel = "Bartlett",
         itermax = 500,
         vcov  = "TrueFixed", prewhite = 0)
m
summary(m)

q <- qr(g0/sqrt(nrow(g0)))
w1 <- matrix(NA, ncol(g0),ncol(g0))
w1[q$pivot, q$pivot] <- chol2inv(q$qr)
m <- gmm(g, dt, t0 = theta,
         weightsMatrix = w1,
         vcov  = "TrueFixed")
m


### An intercept
g <- function(theta, dt) {
  W <- data.matrix(dt[,1])
  Y <- data.matrix(dt[,2])
  X <- data.matrix(dt[,3:ncol(dt)])

  gamma <- theta[1]
  bW <- theta[2]
  bY <- theta[3]

  idx <- 3
  B <- theta[idx + 1]
  s2U <- theta[idx + 2]
  s2V <- theta[idx + 3]
  s2R <- theta[idx + 4]


  a <- gamma + B
  eW <- W - X %*% (gamma * bY + bW)
  eY <- Y - X %*% bY
  Q <- W - Y %*% gamma - X %*% bW
  P <- W - Y %*% a + X %*% (B * bY - bW)
  eY2 <- eY^2
  uYY <- s2U + s2V #u_{yy}
  uYW <- B * s2U + gamma * uYY #u_{yw}

  # Moments. See Appendix 2 of Lewbel et al. (2020), eq. 74-77.
  m1 <- eY * eW - uYW
  m2 <- eY2 - uYY
  m3 <- Q^2 - B^2 * s2U - s2R
  m4 <- Q * P * eY
  m5 <- Q * P * (eY2 - uYY) - 2 * B * s2U * P * eY
  m6 <- eW[, 1] * X
  m7 <- eY[, 1] * X
  M <- cbind(m1, m2, m3, m4, m5, m6, m7)
  M
}
AcemogluJohnson2007$constant <- 1
dt <- AcemogluJohnson2007[, c("loggdppcmadd", "loglifeexpect", "constant")]
theta <- c(gamma = -1.3, b1 = 1.3, d1 = 1.75, beta = 2, s2U = 0.02, s2V = 0.005, s2R = 0.08)
m <- gmm(g, dt, t0 = theta,
         control=list(fnscale=1e-8), type = "iterative", method = "L-BFGS-B",
         lower = c(-Inf, rep(1e-3, 4)),
         upper = rep(Inf, 5))
m








library(momentfit)
g <- list(Y = loglifeexpect ~ 1 + loggdppcmadd301980,
          W = loggdppcmadd ~ 1 + loggdppcmadd301980 + loglifeexpect)
h <- list(~ compsjmhatit)
smod1 <- sysMomentModel(g, h, vcov = "iid", data = AcemogluJohnson2007)
gmmFit(smod1[1])



## Lewbel2020 with covariates
# Table 7, column 6
g <- function(theta, dt) {
  W <- data.matrix(dt[,1])
  Y <- data.matrix(dt[,2])
  X <- data.matrix(dt[, -c(1:2)])
  k <- ncol(X)

  bW <- theta[1:k]
  gamma <- theta[k + 1]
  bY <- theta[k + 1 + 1:k]

  idx <- 2*k + 1
  B <- theta[idx + 1]
  s2U <- theta[idx + 2]
  s2V <- theta[idx + 3]
  s2R <- theta[idx + 4]

  # Conveinience variables
  a <- gamma + B
  eW <- W - X %*% (gamma * bY + bW)
  eY <- Y - X %*% bY
  Q <- W - Y %*% gamma - X %*% bW
  P <- W - Y %*% a + X %*% (B * bY - bW)
  eY2 <- eY^2
  uYY <- s2U + s2V #u_{yy}
  uYW <- B * s2U + gamma * uYY #u_{yw}

  # Moments. See Appendix 2 of Lewbel et al. (2020), eq. 74-77.
  m1 <- eY * eW - uYW
  m2 <- eY2 - uYY
  m3 <- Q^2 - B^2 * s2U - s2R
  m4 <- Q * P * eY
  m5 <- Q * P * (eY2 - uYY) - 2 * B * s2U * P * eY
  m6 <- eW[, 1] * X
  m7 <- eY[, 1] * X
  M <- cbind(m1, m2, m3, m4, m5, m6, m7)
  M
}
AcemogluJohnson2007$constant <- 1
dt <- AcemogluJohnson2007[, c("loggdppcmadd", "loglifeexpect", "constant", "loggdppcmadd301980")]
theta <- c(bWa  = -0.139, bWb = 0.150, gamma = -1.376, bYa = 1.758, bYb = -0.184, beta = 2.090, s2U = 0, s2V = 0.014, s2R = 0.122)
m <- gmm(g, dt, t0 = theta,
         vcov = "iid",
        control=list(fnscale=1e-8), type = "iterative", method = "L-BFGS-B",
        lower = c(-Inf, rep(1e-3, 4)),
        upper = rep(Inf, 5))
coef(m)
m <- gmm(g, dt, t0 = theta, vcov = "iid")
summary(m)


g <- function(theta, dt) {
  W <- data.matrix(dt[,1])
  Y <- data.matrix(dt[,2])
  X <- data.matrix(dt[,3:ncol(dt)])

  gamma <- theta[1]
  bW <- theta[2:3]
  bY <- theta[4:5]

  idx <- 5
  B <- exp(theta[idx + 1])
  s2U <- exp(theta[idx + 2])
  s2V <- exp(theta[idx + 3])
  s2R <- exp(theta[idx + 4])


  a <- gamma + B
  eW <- W - X %*% (gamma * bY + bW)
  eY <- Y - X %*% bY
  Q <- W - Y %*% gamma - X %*% bW
  P <- W - Y %*% a + X %*% (B * bY - bW)
  eY2 <- eY^2
  uYY <- s2U + s2V #u_{yy}
  uYW <- B * s2U + gamma * uYY #u_{yw}

  # Moments. See Appendix 2 of Lewbel et al. (2020), eq. 74-77.
  m1 <- eY * eW - uYW
  m2 <- eY2 - uYY
  m3 <- Q^2 - B^2 * s2U - s2R
  m4 <- Q * P * eY
  m5 <- Q * P * (eY2 - uYY) - 2 * B * s2U * P * eY
  m6 <- eW[, 1] * X
  m7 <- eY[, 1] * X
  M <- cbind(m1, m2, m3, m4, m5, m6, m7)
  M
}

theta <- c(gamma = -1.3, b1 = -0.139, b2 = 0.15, d1 = 1.75, d2 = -0.184,
           beta = log(2), s2U = log(0.02), s2V = log(0.005), s2R = log(0.08))
m <- gmm(g, dt, t0 = theta,
         control=list(fnscale=1e-8), type = "twoStep")
summary(m)


g <- function(tet, x)
{
  m1 <- (tet[1] - x)
  m2 <- (tet[2]^2 - (x - tet[1])^2)
  m3 <- x^3 - tet[1]*(tet[1]^2 + 3*tet[2]^2)
  f <- cbind(m1, m2, m3)
  return(f)
}

gmm(g, x, c(0, 0)))














### Example from vignette
data(Finance)
r <- Finance[1:500, 1:5]
rm <- Finance[1:500, "rm"]
rf <- Finance[1:500, "rf"]
z <- as.matrix(r - rf)
zm <- as.matrix(rm - rf)
res <- gmm(z ~ zm, x = zm)
coef(res)

g5 <- function(tet, x) {
  Y <- x[,1]
  X <- x[, 2:6]
  b1 <- tet[1]
  b2 <- tet[2]
  gmat <- (tet[1] + tet[2] * (1 + Y)) * (1 + X) - 1
  gmat
}
res <- gmm(g5, cbind(rm, r), c(0,0))
specTest(res)


g <- function(theta, x) {
  Y <- x[,1]
  X <- x[,2:3]
  b <- theta[1:2]
  e <- Y - X %*% b
  m <- e[,1] * X
  M <- cbind(m)
  return(M)
}
gmm(g, cbind(rm, 1, rf), t0=c(0,0), method = "BFGS", control=list(fnscale=1e-8))
lm(rm ~ rf, data = data.frame(cbind(rm, 1, rf)))

g <- function(theta, x) {
  m.1 <- x[,"rm"] - theta[1] - theta[2]*x[,"rf"]
  m.z <- (x[,"rm"] - theta[1] - theta[2]*x[,"rf"])*x[,"rf"]
  f <- cbind(m.1, m.z)
  return(f)
}
gmm(g, cbind(rm, rf), t0=c(0,0), method = "BFGS", control=list(fnscale=1e-8))
