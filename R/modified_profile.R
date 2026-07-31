.focus_tangent_coordinates <- function(gradient) {
    gradient <- as.numeric(gradient)
    gradient_norm2 <- sum(gradient^2)
    if (!is.finite(gradient_norm2) ||
        gradient_norm2 <= .Machine$double.eps)
        stop("The focus gradient must be nonzero along the profile.")
    direction <- gradient / gradient_norm2
    Q <- qr.Q(qr(matrix(gradient, ncol = 1L)), complete = TRUE)
    tangent <- Q[, -1L, drop = FALSE]
    list(tangent = tangent,
         transform = cbind(direction, tangent))
}

.focus_log_determinant <- function(x, positive = FALSE) {
    if (nrow(x) == 0L)
        return(0)
    out <- determinant(x, logarithm = TRUE)
    modulus <- as.numeric(out$modulus)
    if (!is.finite(modulus) || (positive && out$sign <= 0))
        stop("Could not compute a valid log determinant.")
    modulus
}


#' Modified profile likelihood and modified signed likelihood-ratio statistic
#'
#' Construct an ordinary profile likelihood for a scalar focus and estimate
#' the nuisance-parameter and information adjustments used for modified
#' profile likelihood and inference based on the modified signed
#' likelihood-ratio statistic \eqn{r^*}.
#'
#' @inheritParams profile_focus
#' @param simulate A function that simulates one dataset at a supplied
#'     parameter vector. It must take the parameter vector as its first
#'     argument.
#' @param nsim Integer. Number of common simulated datasets used to estimate
#'     the sample-space covariance quantities.
#' @param parallelize Logical. If `TRUE`, use
#'     [future.apply::future_lapply()] for the simulation calculations. This
#'     requires the suggested package \pkg{future.apply}.
#'
#' @return
#' An object inheriting from `"profile_focus_list"`. In addition to the
#' columns returned by [profile_focus()], it contains:
#' \describe{
#'   \item{`signed_root`}{The ordinary signed likelihood root
#'     \eqn{r(\psi)}, with sign determined by
#'     \eqn{\hat\psi-\psi}.}
#'   \item{`log_C_inverse`}{The logarithm of the inverse
#'     modified-profile adjustment factor, \eqn{\log C_\psi^{-1}}.}
#'   \item{`NP`}{The nuisance-parameter adjustment
#'     \eqn{NP(\psi)=\log C_\psi^{-1}/r(\psi)}.}
#'   \item{`u_tilde`}{The standardized sample-space derivative
#'     \eqn{\tilde u_\psi} after eliminating the nuisance directions.
#'     Its sign is set equal to that of \eqn{r(\psi)}.}
#'   \item{`log_u_over_r`}{The intermediate quantity
#'     \eqn{\log|\tilde u_\psi/r(\psi)|}.}
#'   \item{`INF`}{The information adjustment
#'     \eqn{INF(\psi)=\log|\tilde u_\psi/r(\psi)|/r(\psi)}.}
#'   \item{`rstar`}{The modified signed likelihood-ratio statistic
#'     \eqn{r^*(\psi)=r(\psi)+NP(\psi)+INF(\psi)}.}
#'   \item{`modified_loglik`}{The modified profile log likelihood
#'     \eqn{\ell_{\mathrm{MP}}(\psi)=
#'     \ell_P(\psi)-\log C_\psi^{-1}}.}
#' }
#'
#' The `"nsim"` attribute records the number of common simulated datasets.
#' The `"expected_information"` attribute contains the simulation estimate
#' of \eqn{\operatorname{Var}_{\hat\theta}\{U(\hat\theta;Y)\}}.
#'
#' @details
#' `likelihood_args` must contain a named `data` element holding the observed
#' dataset. This element is used for the ordinary profile and observed
#' information. It is not passed to `simulate`; the object returned by
#' `simulate` replaces it when evaluating `loglik` and `score` on simulated
#' datasets.
#'
#' The same simulated datasets are used at every point on the profile. The
#' adjustments are undefined at the ordinary maximum in their direct
#' representations because the ordinary signed likelihood root is zero.
#' Quantities that are not evaluated there are stored as `NA`;
#' `log_C_inverse` is set to zero so that `modified_loglik` remains defined.
#' Confidence intervals based on `rstar` are obtained by interpolation away
#' from that point.
#'
#' `u_tilde` and `log_u_over_r` are intermediate quantities and are not
#' parameter estimates or standard errors. They are retained to make the
#' calculation of `INF` reproducible and to help diagnose numerical behavior
#' near the ordinary maximum. Most inferential use will involve
#' `modified_loglik`, `NP`, `INF`, and `rstar` directly.
#'
#' @export
modified_profile_focus <- function(loglik,
                                   score = NULL,
                                   information = NULL,
                                   simulate,
                                   mle,
                                   on = function(theta) theta[1],
                                   on_gradient = NULL,
                                   on_hessian = NULL,
                                   likelihood_args = list(),
                                   nsim = 1000,
                                   parallelize = FALSE,
                                   nleqslv_args = list(),
                                   grid_size = 20,
                                   max_level = 0.995,
                                   approach = c("VM", "focus_grid"),
                                   focus_range = NULL,
                                   auglag_args = list(),
                                   ...) {
    if (!is.function(simulate))
        stop("`simulate` must be a function.")
    if (!is.numeric(nsim) || length(nsim) != 1L || !is.finite(nsim) ||
        nsim < 2L || nsim != as.integer(nsim))
        stop("`nsim` must be an integer of at least 2.")
    if (!is.logical(parallelize) || length(parallelize) != 1L ||
        is.na(parallelize))
        stop("`parallelize` must be `TRUE` or `FALSE`.")
    if (parallelize && !requireNamespace("future.apply", quietly = TRUE))
        stop("Package `future.apply` is required when `parallelize = TRUE`.")
    dots <- list(...)
    model <- .focus_component_functions(loglik,
                                        score,
                                        information,
                                        simulate,
                                        likelihood_args,
                                        require_data = TRUE)
    profile <- do.call(
        profile_focus,
        c(list(loglik = loglik,
               score = score,
               information = information,
               mle = mle,
               on = on,
               on_gradient = on_gradient,
               on_hessian = on_hessian,
               likelihood_args = likelihood_args,
               nleqslv_args = nleqslv_args,
               grid_size = grid_size,
               max_level = max_level,
               approach = approach,
               focus_range = focus_range,
               auglag_args = auglag_args),
          dots)
    )
    theta_hat <- as.numeric(mle)
    theta_profile <- as.matrix(profile$theta)
    p <- length(theta_hat)
    npoints <- nrow(theta_profile)
    on_value <- function(theta)
        do.call(on, c(list(theta), dots))
    on_grad <- if (is.null(on_gradient)) {
        function(theta)
            numDeriv::grad(on_value, theta)
    } else {
        function(theta)
            do.call(on_gradient, c(list(theta), dots))
    }
    on_hess <- if (is.null(on_hessian)) {
        function(theta)
            numDeriv::hessian(on_value, theta)
    } else {
        function(theta)
            do.call(on_hessian, c(list(theta), dots))
    }
    score_at <- if (is.null(score)) {
        function(theta, data)
            numDeriv::grad(model$loglik, theta, data = data)
    } else {
        model$score
    }
    information_at <- if (is.null(information)) {
        if (is.null(score)) {
            function(theta, data)
                -numDeriv::hessian(model$loglik, theta, data = data)
        } else {
            function(theta, data)
                -numDeriv::jacobian(model$score, theta, data = data)
        }
    } else {
        model$information
    }
    coordinates_hat <- .focus_tangent_coordinates(on_grad(theta_hat))
    N_hat <- coordinates_hat$tangent
    T_hat <- coordinates_hat$transform
    J_hat <- information_at(theta_hat, model$data)
    if (!is.numeric(J_hat) ||
        !identical(dim(J_hat), c(p, p)) ||
        any(!is.finite(J_hat)))
        stop("`information` must return a finite numeric matrix with one row ",
             "and one column per parameter.")
    J_hat_c <- crossprod(T_hat, J_hat %*% T_hat)
    J_hat_c_inv <- solve(J_hat_c)
    adjusted_information <- 1 / J_hat_c_inv[1L, 1L]
    if (!is.finite(adjusted_information) || adjusted_information <= 0)
        stop("The adjusted information for the focus must be positive.")
    K_hat <- crossprod(N_hat, J_hat %*% N_hat)
    log_K_hat <- .focus_log_determinant(K_hat, positive = TRUE)
    simulate_one <- function(i) {
        data <- model$simulate(theta_hat)
        score_hat <- score_at(theta_hat, data)
        scores <- t(matrix(
            vapply(
                seq_len(npoints),
                function(j) score_at(theta_profile[j, ], data),
                numeric(p)
            ),
            nrow = p
        ))
        loglik_hat <- model$loglik(theta_hat, data)
        loglik_profile <- vapply(
            seq_len(npoints),
            function(j) model$loglik(theta_profile[j, ], data),
            numeric(1)
        )
        list(score_hat = score_hat,
             scores = scores,
             loglik_hat = loglik_hat,
             loglik_profile = loglik_profile)
    }
    simulations <- if (parallelize) {
        future.apply::future_lapply(seq_len(nsim),
                                    simulate_one,
                                    future.seed = TRUE)
    } else {
        lapply(seq_len(nsim), simulate_one)
    }
    score_hat <- do.call(rbind, lapply(simulations, `[[`, "score_hat"))
    I_hat <- cov(score_hat)
    if (p == 1L)
        I_hat <- matrix(I_hat, 1L, 1L)
    Iinv_J <- solve(I_hat, J_hat)
    loglik_hat <- vapply(simulations, `[[`, numeric(1), "loglik_hat")
    signed_root <- .profile_signed(profile)
    log_C_inverse <- NP <- u_tilde <- log_u_over_r <- INF <- rstar <-
        rep(NA_real_, npoints)
    for (j in seq_len(npoints)) {
        theta_j <- theta_profile[j, ]
        coordinates_j <- .focus_tangent_coordinates(on_grad(theta_j))
        N_j <- coordinates_j$tangent
        T_j <- coordinates_j$transform
        J_j <- information_at(theta_j, model$data)
        if (!is.numeric(J_j) ||
            !identical(dim(J_j), c(p, p)) ||
            any(!is.finite(J_j)))
            stop("Could not evaluate a valid information matrix at profile ",
                 "point ", j, ".")
        if (p == 1L) {
            log_C_inverse[j] <- 0
            NP[j] <- 0
        } else {
            K_j <- crossprod(
                N_j,
                (J_j + profile$lagrange[j] * on_hess(theta_j)) %*% N_j
            )
            scores_j <- do.call(
                rbind,
                lapply(simulations, function(x) x$scores[j, ])
            )
            S_j <- cov(scores_j, score_hat)
            A_j <- S_j %*% Iinv_J
            D_j <- crossprod(N_j, A_j %*% N_hat)
            log_C_inverse[j] <-
                .focus_log_determinant(D_j) -
                0.5 * .focus_log_determinant(K_j, positive = TRUE) -
                0.5 * log_K_hat
            if (signed_root[j] != 0)
                NP[j] <- log_C_inverse[j] / signed_root[j]
        }
        if (signed_root[j] != 0) {
            if (p == 1L) {
                scores_j <- do.call(
                    rbind,
                    lapply(simulations, function(x) x$scores[j, ])
                )
                S_j <- cov(scores_j, score_hat)
                A_j <- S_j %*% Iinv_J
            }
            A_j_c <- crossprod(T_j, A_j %*% T_hat)
            loglik_j <- vapply(
                simulations,
                function(x) x$loglik_profile[j],
                numeric(1)
            )
            q_j <- cov(loglik_hat - loglik_j, score_hat)
            b_j <- drop(q_j %*% Iinv_J %*% T_hat)
            nuisance_adjustment <- if (p == 1L) {
                0
            } else {
                D_j <- A_j_c[-1L, -1L, drop = FALSE]
                drop(b_j[-1L] %*%
                     solve(D_j, A_j_c[-1L, 1L, drop = FALSE]))
            }
            u_raw <- (b_j[1L] - nuisance_adjustment) /
                sqrt(adjusted_information)
            u_tilde[j] <- sign(signed_root[j]) * abs(u_raw)
            if (is.finite(u_tilde[j]) && u_tilde[j] != 0) {
                log_u_over_r[j] <-
                    log(abs(u_tilde[j] / signed_root[j]))
                INF[j] <- log_u_over_r[j] / signed_root[j]
                rstar[j] <- signed_root[j] + NP[j] + INF[j]
            }
        } else {
            log_C_inverse[j] <- 0
        }
    }
    profile$signed_root <- signed_root
    profile$log_C_inverse <- log_C_inverse
    profile$NP <- NP
    profile$u_tilde <- u_tilde
    profile$log_u_over_r <- log_u_over_r
    profile$INF <- INF
    profile$rstar <- rstar
    profile$modified_loglik <- profile$loglik - log_C_inverse
    class(profile) <- c("modified_profile_focus_list", class(profile))
    attr(profile, "nsim") <- nsim
    attr(profile, "expected_information") <- I_hat
    profile
}
