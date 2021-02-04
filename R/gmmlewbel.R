#' GMM estimator of structural equations under heteroskedastic or non-normal
#' errors.
#'
#' @import gmm
#' @import Formula
#'
#' @param formula a multipart formula like y1 | y2 ~ x1 + x2 | z. See Details
#'   for more information.
#' @param data a data.frame or object convertable to a matrix.
#' @param init a vector of initial parameters. Estimated by OLS if not provided.
#' @param overid should additional moments be used? Ignored if error_type is not
#'   "nonnormal".
#' @param error_type Are errors heteroskedastic or non-normal? This determines
#'   the method used.
#' @return a gmm object.
gmmlewbel <- function(formula, data, init = NULL,
                      error_type = c("hetero", "nonnormal"),
                      is_simultaneous = FALSE, use_overid = FALSE,
                      exp_variances = TRUE,
                      ...) {
  error_type <- match.arg(error_type)

  f <- Formula(formula)
  data <- as.data.frame(data)

  # Get the model matrices
  W <- model.part(f, data = data, lhs = 1)
  Y <- model.part(f, data = data, rhs = 1)
  df <- cbind(W, Y)

  w_var <- names(W)
  y_vars <- names(Y)

  # Get X if present
  if(length(f)[[2]] > 1) {
    X <- model.matrix(f, data = data, rhs = 2)
    nX <- ncol(X)
    x_vars <- setdiff(colnames(X), c("(Intercept)"))
    df <- cbind(df, X)
  } else {
    X <- NULL
    nX <- 0
    x_vars <- c()
  }

  # Handle the instruments.
  # Z1 are endogenous instruments from X
  # Z2 are exogenous instruments
  if(length(f)[[2]] > 2) {
    Z <- model.part(f, data = data, rhs = 3)
    z1_vars <- names(Z)[names(Z) %in% x_vars]
    z2_vars <- setdiff(names(Z), z1_vars) # exogenous instruments

    nZ1 <- length(z1_vars)
    nZ2 <- length(z2_vars)

    if(nZ1 > 0) {
      Z1 <- as.matrix(Z[, z1_vars], nrow = nrow(Z)) # Instruments to be constructed
      colnames(Z1) <- z1_vars
    }

    if(nZ2 > 0) {
      Z2 <- as.matrix(Z[, z2_vars], nrow = nrow(Z)) # exogenous instruments
      colnames(Z2) <- z2_vars
    }

  } else {
    Z1 <- NULL
    Z2 <- NULL
    nZ1 <- 0
    nZ2 <- 0
  }

  # Different model types
  if(error_type == "hetero") {

    # Via explicit instrument construction:
    # 1. construct instruments as (Z - u_Z)e, where e is the residual from the
    #    1st stage
    # 2. Run GMM with the constructed instruments

    # Construct the instruments
    f1 <- as.formula(paste0(y_vars, " ~ ",
                            paste(c(x_vars), collapse = " + ")))
    Z1a <- residuals(lm(f1, data = df)) * scale(Z1, scale = FALSE)
    z1a_vars <- paste0("inst_", z1_vars)   # constructed instruments
    colnames(Z1a) <- z1a_vars

    # Formulas for GMM
    f1 <- as.formula(paste0(w_var, " ~ ",
                            paste(c(y_vars, x_vars), collapse = " +")))
    f2 <- as.formula(paste0(" ~ ",
                            paste(c(x_vars, z1a_vars, z2_vars),
                                  collapse = " + ")))

    #obj <- gmm(f1, f2, data = cbind(df, Z1a, Z2), vcov = "iid")

    # Alternative method using specified moments
    moment_fn <- moments_hetero(nX = nX, nZ = nZ1 + nZ2,
                                is_simultaneous)
    obj <- gmm(moment_fn, cbind(df, Z1, Z2), vcov = "iid",
               t0 = rep(1, 1 + 2*nX + nZ1 + nZ2 + is_simultaneous),
               type = type, ...)

    # # estimated. Instruments are not directly instrumented in the first stage
    # if(is.null(init)) {
    #   # Parameter order is: c(b1, g1, b2, g2, u)
    #   # Y1 regression should always include Y2
    #   r1_cols <- c(names(Y1), colnames(X), names(Y2))
    #   f1_str <- paste0(names(Y1), " ~ . - 1")
    #   cr1 <- coef(lm(as.formula(f1_str),
    #                  data = df[, r1_cols]))
    #
    #   # Y2 regression should include X and may include Y1
    #   r2_cols <- c(names(Y2), colnames(X))
    #   if(is_simultaneous) r2_cols <- c(r2_cols, names(Y1))
    #   f2_str <- paste0(names(Y2), " ~ . - 1")
    #   cr2 <- coef(lm(formula(f2_str), data = df[, r2_cols]))
    #
    #   init <- c(cr1, cr2, colMeans(Z))
    # }
    #
    # moment_fn <- moments_hetero(nX = ncol(X), nZ = nZ, is_simultaneous)

  } else if(error_type == "nonnormal") {
    # In the nonnormal method, "external" instruments are part of the Y2
    # equation.
    if(is.null(init)) {
      # Parameter order is c(b1, g, b2, d, b, tU, tV, tR, tY1)

      # Y1 regression should always include Y2
      r1_cols <- c(names(W), names(Y), colnames(X))
      f1_str <- paste0(names(Y), " ~ . - 1")
      cr1 <- coef(lm(as.formula(f1_str),
                     data = df[, r1_cols]))

      # Y2 regression should include X and Z if there are any Z
      r2_cols <- c(names(Y), colnames(X))
      if(nZ2 > 0) r2_cols <- c(r2_cols, colnames(Z2))
      f2_str <- paste0(names(Y), " ~ . - 1")
      cr2 <- coef(lm(formula(f2_str), data = df[, r2_cols]))

      # We also need initial values for b and the taus.
      # To enforce >= 0 requirement, we let all be close to zero and
      # defined as b = exp(tau_b) for example.
      init <- c(cr1, cr2, rep(0.01, 4))
      names(init)[length(init) - 3:0] <- c("t_beta", "t_U", "t_V", "t_R")

      # If overidentified, we also need estimates for tau_Y1
      if(use_overid) init <- c(init, 0.01)
      names(init)[length(init)] <- "t_W"
    }

    #If no covariates, demean the outcomes
    if(nX == 0) df <- scale(df, scale = FALSE)

    moment_fn <- moments_nonnormal(nX = nX, nZ = nZ2,
                                   overid = use_overid,
                                   exp_variances = exp_variances)
    dt <- cbind(as.matrix(df), Z2)
    obj <- gmm(moment_fn, dt, t0 = init, ...)
  }

  obj
}
