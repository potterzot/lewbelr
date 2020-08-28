#' Moment function for heteroskedastic errors.
moments_hetero <- function(nX, nZ) {
  function(theta, df) {
    Y1 <- data.matrix(df[,1])
    Y2 <- data.matrix(df[,2])
    X <- data.matrix(df[, 2 + 1:nX]) # Includes intercept
    Z <- data.matrix(df[, 2 + nX + 1:nZ]) # Does not include intercept

    # b1 and b2 include the gamma for each Y
    b1 <- theta[1:nX]
    gamma1 <- theta[nX + 1]
    b2 <- theta[nX + 1 + 1:nX]

    # We store the second gamma last to allow setting to 0
    gamma2 <- ifelse(length(theta) == 2*(1 + nX), theta[length(theta)], 0)

    # Moments as in section 3.1 of Lewbel2012
    err1 <- Y1 - X %*% b1 - gamma1 * Y2
    err2 <- Y2 - X %*% b2 - gamma2 * Y1

    m1 <- err1[, 1] * X
    m2 <- err2[, 1] * X
    m3 <- sweep(Z, 2, colMeans(Z))
    m4 <- m3 * err1[, 1] * err2[, 1]
    cbind(m1, m2, m3, m4)
  }
}

#' Moment function for a triangular system with no instruments.
moments_nonnormal <- function(nX, nZ) {
  function(theta, df) {
    Y1 <- data.matrix(df[,1])
    Y2 <- data.matrix(df[,2])
    X <- data.matrix(df[, 2 + 1:nX]) # Includes intercept
    Z <- data.matrix(df[, 2 + nX + 1:nZ]) # Does not include intercept

    b2 <- theta[1:nX]
    g <- theta[nX + nY1]
    b1 <- theta[nY1 + nX + 1:nX]
    b <- exp(theta[nY1 + 2*nX + 1])
    s2U <- exp(theta[nY1 + 2*nX + 2])
    s2V <- exp(theta[nY1 + 2*nX + 3])
    s2R <- exp(theta[nY1 + 2*nX + 4])

    # Get values for convenience variables
    Y1t <- Y1 - X %*% b1 - Z %*% d
    Y2t <- Y2 - X %*% (g * b1 + b2)
    Q <- Y2 - Y1 %*% g - X %*% b2
    P <- Y2 - Y1 %*% (g + b) + X %*% (b*b1 - b2)

    # moments
    m1 <- Y1t * Y2t - b * s2U - g * (s2U + s2V)
    m2 <- Y1t^2 - s2U - s2V
    m3 <- Q^2 - b^2 * s2U - s2R
    m4 <- Q * P * Y1t
    m5 <- Q * P * (Y1t^2 - s2U - s2V) - 2 * b * s2U * P * Y1t
    m6 <- Q[, 1] * X
    m7 <- Y1t[, 1] * X
    M <- cbind(m1, m2, m3, m4, m5, m6, m7)
    return(M)
  }
}
