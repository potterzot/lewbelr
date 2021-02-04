

### Non exponentiated ----
g <- function(theta, dt) {
  W <- data.matrix(dt[,1])
  Y <- data.matrix(dt[,2])
  gamma <- theta[1]

  idx <- 1
  B <- theta[idx + 1]
  s2U <- theta[idx + 2]
  s2V <- theta[idx + 3]
  s2R <- theta[idx + 4]
  #uWW <- theta[idx + 5]

  uYY <- s2U + s2V
  uYW <- B * s2U + gamma * uYY

  a <- gamma + B
  Q <- W - Y * gamma
  P <- W - Y * a
  YY <- Y * Y
  YW <- Y * W

  # Moments. See Appendix 2 of Lewbel et al. (2020), eq. 74-77.
  m1 <- YW - uYW
  m2 <- YY - uYY
  m3 <- Q * P * Y
  m4 <- Q * P * m2 - 2 * B * s2U * P * Y
  m5 <- Q^2 - B^2 * s2U - s2R
  M <- cbind(m1, m2, m3, m4, m5)
  M
}


res <- lapply(1:100, function(i) {
  dt <- make_simdata(method = 2)

  theta <- c(gamma = 1, beta = 1, tU = 1.72, tV = 1.64, tR = 1.64)
  m <- gmm(g, dt, t0 = theta,
           itermax = 500,
           vcov = "iid",
           method = "L-BFGS-B",
           lower = c(-Inf, rep(1e-3, 4)),
           upper = rep(Inf, 5),
           control=list(
             fnscale=1e-8
           ))
  m
  c(coef(m))
}) %>%
  dplyr::bind_rows()
colMeans(res)



#### Overidentified ---------------------------
g <- function(theta, dt) {
  W <- data.matrix(dt[,1])
  Y <- data.matrix(dt[,2])
  gamma <- theta[1]

  idx <- 1
  B <- exp(theta[idx + 1])
  s2U <- exp(theta[idx + 2])
  s2V <- exp(theta[idx + 3])
  s2R <- exp(theta[idx + 4])
  uWW <- exp(theta[idx + 5])

  uYY <- s2U + s2V
  uYW <- B * s2U + gamma * uYY

  a <- gamma + B
  Q <- W - Y * gamma
  P <- W - Y * a
  YY <- Y * Y
  YW <- Y * W
  WW <- W * W
  Y3 <- Y^3
  Y5a <- Y^5 - 10 * Y3 * uYY

  # Moments. See Appendix 2 of Lewbel et al. (2020), eq. 74-77.
  m1 <- YW - uYW
  m2 <- YY - uYY
  m3 <- Q * P * Y
  m4 <- Q * P * m2 - 2 * B * s2U * P * Y
  m5 <- Q^2 - B^2 * s2U - s2R
  m6 <- WW - uWW
  m7 <- -3 * Y * WW * uYY - 6 * uYW * YY * W - uWW * Y3 +
    Y3 * WW - a^2 * Y5a - (gamma + B) * (-6 * uYY * YY * W - 4 * uYW * Y3 + Y^4 * W - a * Y5a)
  M <- cbind(m1, m2, m3, m4, m5, m6, m7)
  M
}

res <- lapply(1:100, function(i) {
  dt <- make_simdata(method = 2)
  theta <- c(gamma = 1, tB = log(2.7),
             tU = log(1.72), tV = log(1.64), tR = log(1.64), tW = log(9.54))

  g0 <- g(theta, dt)
  g0 <- scale(g0, scale = FALSE)
  class(g0) <- "gmmFct"
  v0 <- crossprod(g0)
  w0 <- corpcor::pseudoinverse(v0)
  m <- gmm(g, dt, t0 = theta,
           itermax = 500,
#           vcov = "TrueFixed",
#           weightsMatrix = w0,
           #method = "BFGS",
           #lower = c(-Inf, rep(1e-3, 4)),
           #upper = rep(Inf, 5),
           control=list(
             fnscale=1e-8
           ))
  c(coef(m)[1], exp(coef(m)[-1]))
}) %>%
  dplyr::bind_rows()
colMeans(res)



dt <- data.frame(
  u = exp(rnorm(N, -0.5, 1)),
  v = runif(N, -2, 2),
  r = runif(N, -2, 2)) %>%
  dplyr::mutate(
    y = u + v,
    w = gamma * y + beta * u + r) %>%
  dplyr::select(w, y, u, v, r) %>%
  dplyr::mutate(
    y = scale(y, scale = FALSE),
    w = scale(w, scale = FALSE)
  )

g <- function(theta, dt) {
  W <- data.matrix(dt[,1])
  Y <- data.matrix(dt[,2])
  gamma <- theta[1]

  idx <- 1
  B <- exp(theta[idx + 1])
  s2U <- exp(theta[idx + 2])
  s2V <- exp(theta[idx + 3])
  s2R <- exp(theta[idx + 4])
  #uWW <- theta[idx + 5]

  uYY <- s2U + s2V
  uYW <- B * s2U + gamma * uYY

  a <- gamma + B
  Q <- W - Y * gamma
  P <- W - Y * a
  YY <- Y * Y
  YW <- Y * W

  # Moments. See Appendix 2 of Lewbel et al. (2020), eq. 74-77.
  m1 <- YW - uYW
  m2 <- YY - uYY
  m3 <- Q * P * Y
  m4 <- Q * P * m2 - 2 * B * s2U * P * Y
  m5 <- Q^2 - B^2 * s2U - s2R
  M <- cbind(m1, m2, m3, m4, m5)
  M
}
theta <- c(gamma = 1, beta = log(1), tU = log(1.72), tV = log(1.33), tR = log(1.33))
m <- gmm(g, dt, t0 = theta,
         itermax = 500,
         vcov = "iid",
#         type = "iterative",
#         prewhite = 0,
         control=list(
           fnscale=1e-4
         ))
m
c(coef(m)[1], exp(coef(m)[-1]))






g0 <- g(theta, dt)
g0 <- scale(g0, scale = FALSE)
v0 <- crossprod(g0)
class(g0) <- "gmmFct"
v0 <- vcovHAC(g0, sandwich = FALSE)
w0 <- solve(v0)
m <- gmm(g, dt, t0 = theta,
         itermax = 500,
         vcov = "TrueFixed",
         weightsMatrix = w0,
         control=list(
           fnscale=1e-8
         ))
m
c(coef(m)[1], exp(coef(m)[-1]))

,
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
