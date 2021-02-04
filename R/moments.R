#' Moment function for heteroskedastic errors.
moments_hetero <- function(nX, nZ, is_simultaneous) {
  function(theta, dt) {
    Y1 <- data.matrix(dt[,1])
    Y2 <- data.matrix(dt[,2])
    X <- data.matrix(dt[, 2 + 1:nX])
    Z <- data.matrix(dt[, 2 + nX + 1:nZ]) # constructed and exogenous

    # b1 is coefficients on X
    # gamma1 is coefficient on Y2
    # b2 is coefficients on X
    # gamma2 is coefficient on Y1 if exists and 0 otherwise
    # u is parameter estimating mean of Z
    gamma1 <- theta[1]
    b1 <- theta[1 + 1:nX]
    gamma2 <- theta[1 + nX + 1]
    gamma2 <- ifelse(!is_simultaneous, 0, theta[1 + nX + 1])
    b2 <- theta[is_simultaneous + 1 + nX + 1 + 1:nX]

    u <- theta[is_simultaneous + 1 + 2*nX + 1 + 1:nZ] # u for constructed Z

    # Moments as in section 3.1 of Lewbel2012
    err1 <- Y1 - X %*% b1 - gamma1 * Y2
    err2 <- Y2 - X %*% b2 - gamma2 * Y1

    m1 <- X * err1[, 1]
    m2 <- X * err2[, 1]
    m3 <- sweep(Z, 2, u)
    m4 <- m3 * err1[, 1] * err2[, 1]
    cbind(m1, m2, m3, m4)
  }
}

moments_nonnormal2 <- function(nX, nZ, overid) {
  function(theta, dt) {
    W <- data.matrix(dt[,1])
    Y <- data.matrix(dt[,2])

    # theta is b1, g,
    n_params <- length(theta)
    gamma <- theta[1]
    idx <- 1
    B <- exp(theta[idx + 1])
    s2U <- exp(theta[idx + 2])
    s2V <- exp(theta[idx + 3])
    s2R <- exp(theta[idx + 4])

    uYY <- s2U + s2V #u_{yy}
    uYW <- B * s2U + gamma * uYY #u_{yw}

    a <- gamma + B
    Q <- W - Y %*% gamma
    P <- W - Y %*% a
    Y2 <- Y^2

    # Moments. See Appendix 2 of Lewbel et al. (2020), eq. 74-77.
    m1 <- Y * W - uYW
    m2 <- Y2 - uYY
    m3 <- Q^2 - (B^2 * s2U) - s2R
    m4 <- Q * P * Y
    m5 <- Q * P * (Y2 - uYY) - 2 * (uYW - gamma * uYY) * P * Y
    M <- cbind(m1, m2, m3, m4, m5)
    M
  }
}



#' Moment function for a triangular system with no instruments.
moments_nonnormal <- function(nX, nZ, overid, exp_variances) {
  function(theta, dt) {
    W <- data.matrix(dt[,1])
    Y <- data.matrix(dt[,2])

    # theta is b1, g,
    n_params <- length(theta)
    gamma <- theta[1]

    if(nX == 0) {
      X <- matrix(rep(0, nrow(W)), ncol = 1)
      bW <- 0
      bY <- 0
    } else {
      X <- data.matrix(dt[, 2 + 1:nX])
      bW <- theta[1 + 1:nX]
      bY <- theta[1 + nX + 1:nX]
    }

    if(nZ == 0) {
      Zd <- 0
    } else {
      d <- theta[1 + 2*nX + 1:nZ]
      Z <- data.matrix(dt[, 2 + nX + 1:nZ])
      Zd <- Z %*% d
    }

    idx <- 1 + 2*nX + nZ
    if(exp_variances) {
      B <- exp(theta[idx + 1])
      s2U <- exp(theta[idx + 2])
      s2V <- exp(theta[idx + 3])
      s2R <- exp(theta[idx + 4])

    } else {
      B <- theta[idx + 1]
      s2U <- theta[idx + 2]
      s2V <- theta[idx + 3]
      s2R <- theta[idx + 4]
    }
    # B <- 2.791
    # s2U <- 0.0248
    # s2V <- 0.00587
    # s2R <- 0.0838

    a <- gamma + B
    eW <- W - X %*% (gamma * bY + bW)
    eY <- Y - X %*% bY - Zd
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
    m5 <- Q * P * m2 - 2 * B * s2U * P * eY
    M <- cbind(m1, m2, m3, m4, m5)
    if(nX > 0) {
      m6 <- eW[, 1] * X
      m7 <- eY[, 1] * X
      M <- cbind(M, m6, m7)
    }
    if(overid) {
      # See Appendix 2 of Lewbel et al. (2020), eq. 78-79
      eW2 <- eW^2
      uWW <- exp(theta[idx + 5])
      eY3 <- eY^3
      eY5a <- eY^5 - 10 * eY3 * uYY
      m8 <- eW2 - uWW
      m9 <- -3 * eY * eW2 * uYY - 6 * uYW * eY2 * eW - uWW * eY3 +
        eY3 * eW2 - a^2 * eY5a -
        (gamma + a) * (-6 * uYY * eY2 * eW - 4 * uYW * eY3 + eY^4 * eW - a * eY5a)
      M <- cbind(M, m8, m9)
    }
    if(nZ > 0) {
      # Eq. 80 from Lewbel et al. (2020) Appendix B
      m10 <- eY[, 1] * Z
      M <- cbind(M, m10)
    }
    M
  }
}
