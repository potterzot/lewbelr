#' GMM estimator of structural equations under heteroskedastic or non-normal
#' errors.
#'
#' @import gmm
#' @import Formula
#'
#' @param formula a multipart formula like y2 ~ x1 + x2 | y1 | z. See details.
#' @param data a data.frame or object convertable to a matrix.
#' @param init a vector of initial parameters. Estimated by OLS if not provided.
#' @return a gmm object.
gmmlewbel <- function(formula, data,
                      init = NULL,
                      error_type = c("hetero", "nonnormal")) {
  error_type <- match.arg(error_type)
  #flist <- list_formulas(formula, env = parent.frame())
  f <- Formula(formula)
  data <- as.data.frame(data)

  # Get the model matrices
  Y1 <- model.part(f, data = data, lhs = 1)
  Y2 <- model.part(f, data = data, lhs = 2)
  X1 <- model.matrix(f, data = data, rhs = 1)
  X2 <- model.matrix(f, data = data, rhs = 2)

  # X must be the same for all cases
  X <- X1[, setdiff(colnames(X1), names(Y2))]
  if(any(X != X2[, setdiff(colnames(X2), names(Y1))])) {
    stop("The same regressors must be used throughout the system. I.e. X1 == X2.")
  }
  Z <- model.part(f, data = data, rhs = 3)

  df <- cbind(Y1, Y2, X, Z)

  if(error_type == "hetero") {
    moment_fn <- moments_hetero(nX = ncol(X), nZ = ncol(Z))
  } else {
    moment_fn <- moments_nonnormal(nX = ncol(X), nZ = ncol(Z))
  }

  # Initial values should be c(b1, g1, b2, g2)
  if(is.null(init)) {
    # First regression should always include Y2
    f1 <- formula(f, lhs = 1, rhs = 1)
    cr1 <- coef(lm(f1, data = df))[c(colnames(X), names(Y2))]

    # Second may not if it is a triangular system
    f2 <- formula(f, lhs = 2, rhs = 2)
    if(any(grepl(names(Y1), f2))) {
      cr2 <- coef(lm(f2, data = df))[c(colnames(X), names(Y1))]
    } else {
      cr2 <- coef(lm(f2, data = df))[colnames(X)]
    }

    init <- c(cr1, cr2)
  }

  obj <- gmm(moment_fn, df, t0 = init)
  obj
}
