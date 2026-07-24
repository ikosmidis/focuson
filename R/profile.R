#' Profile likelihood confidence intervals
#'
#' Compute a profile likelihood confidence interval for a scalar
#' function of a model parameter vector using the endpoint equations
#' of Venzon and Moolgavkar (1988).
#'
#' @param loglik Function returning the log-likelihood at a supplied
#'     parameter vector. It must take the parameter vector as its
#'     first argument.
#' @param score Optional function returning the score vector at a
#'     supplied parameter vector. It must take the parameter vector as
#'     its first argument.  If `NULL`, the score is computed
#'     numerically from `loglik`.
#' @param information Optional function returning the observed or
#'     expected information matrix at a supplied parameter vector. It
#'     must take the parameter vector as its first argument. If
#'     `NULL`, the information is computed numerically as minus the
#'     Hessian of `loglik`.
#' @param mle Numeric vector. The unrestricted maximum likelihood
#'     estimate of the full model parameter vector.
#' @param on Function specifying the scalar parameter of interest. It
#'     must take the parameter vector as its first argument and return
#'     a numeric scalar.
#' @param on_gradient Optional function returning the gradient of `on`
#'     with respect to the parameter vector. It must take the
#'     parameter vector as its first argument. If `NULL`, the gradient
#'     is computed numerically.
#' @param on_hessian Optional function returning the Hessian matrix of
#'     `on` with respect to the parameter vector. It must take the
#'     parameter vector as its first argument. If supplied together
#'     with `information`, a Jacobian is passed to
#'     [nleqslv::nleqslv()].
#' @param likelihood_args List of additional arguments passed to
#'     `loglik`, `score`, and `information`.
#' @param nleqslv_args List of additional arguments passed to
#'     [nleqslv::nleqslv()]. The `x`, `fn`, and `jac` arguments are
#'     determined internally.
#' @param level Confidence level. Default is `0.95`.
#' @param start Optional named list with elements `"lower"` and
#'     `"upper"`, containing starting values for the corresponding
#'     endpoint equations.  Each element must be a numeric vector of
#'     length `length(mle) + 1`, consisting of the parameter vector
#'     followed by the Lagrange multiplier. This argument is intended
#'     primarily for warm-starting repeated calls using the `"solution"`
#'     attribute returned by a previous call. Users will generally want to
#'     leave it `NULL` (default) so that branch-specific Wald-type starting
#'     values are constructed. Inappropriate or identical starting values for
#'     the two branches may cause both endpoints of the interval to be
#'     identical.
#' @param do_checks Logical. If `TRUE` (default), validate the inputs
#'     and their values at `mle`. Set to `FALSE` only when repeatedly
#'     calling `profile_ci()` with inputs that have already been
#'     checked.
#' @param ... Additional arguments passed to `on`, `on_gradient`, and
#'     `on_hessian`.
#'
#' @return
#' A numeric vector of length 2 with names `"lower"` and `"upper"`. The result
#' has attributes:
#' \describe{
#'   \item{`"type"`}{The string `"profile"`.}
#'   \item{`"max|fvec|"`}{The maximum absolute endpoint equation residual for
#'     each endpoint.}
#'   \item{`"messages"`}{The convergence messages returned by
#'     [nleqslv::nleqslv()] for the lower and upper endpoint solves.}
#'   \item{`"loglik"`}{The log-likelihood evaluated at the lower and upper
#'     endpoint parameter vectors.}
#'   \item{`"solution"`}{A list containing the solutions of the lower and upper
#'     endpoint equations. These can be supplied as `start` in a subsequent
#'     call.}
#'   \item{`"iter"`}{The number of outer iterations used by
#'     [nleqslv::nleqslv()] for each endpoint.}
#' }
#'
#' @details
#'
#' The endpoint equations are
#' \deqn{2\{\ell(\hat\theta) - \ell(\theta)\} - c = 0}
#' and
#' \deqn{\nabla \ell(\theta) - \lambda \nabla g(\theta) = 0,}
#' where \eqn{\ell(\theta)} is `loglik`, \eqn{\hat\theta} is `mle`,
#' \eqn{g} is the function supplied through `on`, and \eqn{c} is the
#' `level` quantile of a chi-squared distribution with one degree of
#' freedom.
#'
#' The function solves these equations twice and returns the corresponding
#' values of `on`. Unless `start` is supplied, Wald-type starting values in
#' opposite directions are used. This is a low-level routine: endpoint convergence
#' diagnostics are returned as attributes, and callers can decide how
#' strictly to enforce them.
#'
#' If `information` and `on_hessian` are both supplied, `profile_ci()`
#' uses them to construct a Jacobian for the endpoint equations and
#' passes it to [nleqslv::nleqslv()]. This Jacobian is exact when
#' `information` is minus the derivative of `score`. If either function
#' is omitted, [nleqslv::nleqslv()] computes its own numerical Jacobian.
#' Supplying `score`, `information`, `on_gradient`, and `on_hessian` can
#' substantially reduce computation.
#'
#' @references
#'
#' Venzon D J, Moolgavkar S H (1988). A method for computing
#' profile-likelihood-based confidence intervals. *Journal of the Royal
#' Statistical Society: Series C (Applied Statistics)*, **37**, 87--94.
#'
#' @examples
#'
#' ## Example from ?profile.glm
#' library("enrichwith")
#'
#' budworm <- data.frame(ldose = rep(0:5, 2),
#'                       numdead = c(1, 4, 9, 13, 18, 20, 0, 2, 6, 10, 12, 16),
#'                       sex = factor(rep(c("M", "F"), c(6, 6)))) |>
#'            transform(numalive = 20 - numdead)
#' bw_fit <- glm(cbind(numalive, numdead) ~ sex * ldose, family = binomial,
#'               data = budworm)
#' aux <- get_auxiliary_functions(bw_fit)
#' ll <- function(theta) {
#'     sum(aux$dmodel(coefficients = theta, log = TRUE))
#' }
#' sc <- function(theta) {
#'     aux$score(coefficients = theta)
#' }
#' info <- function(theta) {
#'     aux$information(coefficients = theta)
#' }
#'
#' parameter <- function(theta, j) theta[j]
#' parameter_gradient <- function(theta, j) {
#'     out <- numeric(length(theta))
#'     out[j] <- 1
#'     out
#' }
#' parameter_hessian <- function(theta, j) {
#'     matrix(0, length(theta), length(theta))
#' }
#'
#' sapply(seq_along(coef(bw_fit)), function(j) {
#'     profile_ci(ll, score = sc, information = info,
#'                mle = coef(bw_fit), level = 0.99,
#'                on = parameter,
#'                on_gradient = parameter_gradient,
#'                on_hessian = parameter_hessian,
#'                j = j)
#' })
#'
#' ## compare with ?profile
#' pr <- profile(bw_fit)
#' pr_ci <- confint(pr, level = 0.99)
#' t(pr_ci)
#'
#' ## Profile likelihood interval for the difference in probabilities
#' ## between the levels of `sex`, for a given `ldose`
#' me <- function(theta, ldose) {
#'    tt <- terms(bw_fit) |> delete.response()
#'    df <- data.frame(ldose = ldose, sex = c("M", "F"))
#'    mm <- model.matrix(tt, data = df)
#'    probs <- drop(plogis(mm %*% theta))
#'    diff(probs)
#' }
#'
#' profile_ci(ll, score = sc, mle = coef(bw_fit),
#'            level = 0.95, on = me, ldose = 0)
#' profile_ci(ll, score = sc, mle = coef(bw_fit),
#'            level = 0.95, on = me, ldose = 5)
#'
#' @export
profile_ci <- function(loglik,
                       score = NULL,
                       information = NULL,
                       mle,
                       on = function(theta) theta[1],
                       on_gradient = NULL,
                       on_hessian = NULL,
                       likelihood_args = list(),
                       nleqslv_args = list(),
                       level = 0.95,
                       ...,
                       start = NULL,
                       do_checks = TRUE) {
    do_checks <- isTRUE(do_checks)
    if (do_checks) {
        required <- list(loglik = loglik, on = on)
        optional <- list(score = score,
                         information = information,
                         on_gradient = on_gradient,
                         on_hessian = on_hessian)
        check_required <- !vapply(required, is.function, logical(1))
        if (any(check_required)) {
            stop("The following arguments must be functions: ",
                 paste(sprintf("`%s`", names(check_required)[check_required]), collapse = ", "), ".")
        }
        check_optional <- !vapply(optional, function(x) is.null(x) || is.function(x), logical(1))
        if (any(check_optional)) {
            stop("The following arguments must be functions or `NULL`: ",
                 paste(sprintf("`%s`", names(optional)[check_optional]), collapse = ", "), ".")
        }
        if (!is.numeric(mle) || length(mle) == 0 || any(!is.finite(mle)))
            stop("`mle` must be a finite numeric vector.")
        if (!is.list(likelihood_args))
            stop("`likelihood_args` must be a list.")
        if (!is.list(nleqslv_args))
            stop("`nleqslv_args` must be a list.")
    }
    ll <- function(theta) do.call(loglik, c(list(theta), likelihood_args))
    if (is.null(on_gradient))
        on_gradient  <- function(theta, ...) numDeriv::grad(on, theta, ...)
    if (is.null(score))
        sc <- function(theta) numDeriv::grad(ll, theta)
    else
        sc <- function(theta) do.call(score, c(list(theta), likelihood_args))
    if (is.null(information))
        info <- function(theta) -numDeriv::hessian(ll, theta)
    else
        info <- function(theta) do.call(information, c(list(theta), likelihood_args))
    use_jacobian <- !is.null(information) && !is.null(on_hessian)
    maxloglik <- ll(mle)
    if (do_checks) {
        if (!is.numeric(maxloglik) || length(maxloglik) != 1 || !is.finite(maxloglik)) {
            stop("`loglik` must return a finite numeric scalar at `mle`.")
        }
        score_mle <- sc(mle)
        if (!is.numeric(score_mle) || length(score_mle) != length(mle) || any(!is.finite(score_mle)))
            stop("`score` must return a finite numeric vector of length `length(mle)` at `mle`.")
        if (max(abs(score_mle)) > 1e-02)
            warning("`mle` is not the maximizer of the log-likelihood;",
                    "the maximum absolute value of the gradient of the log-likelihood at `mle` is ",
                    round(max(abs(score_mle)), 2), ".")
        on_grad <- on_gradient(mle, ...)
        if (!is.numeric(on_grad) || length(on_grad) != length(mle) || any(!is.finite(on_grad)))
            stop("`on_gradient` must return a finite numeric vector of length `length(mle)` at `mle`.")
        if (use_jacobian) {
            on_hess_mle <- on_hessian(mle, ...)
            if (!is.numeric(on_hess_mle) ||
                !identical(dim(on_hess_mle), c(length(mle), length(mle))) ||
                any(!is.finite(on_hess_mle))) {
                stop("`on_hessian` must return a finite numeric matrix with one row and one column per parameter.")
            }
        }
    }
    quant <- qchisq(level, 1)
    npars <- length(mle) + 1
    eq <- function(pars) {
        theta <- pars[1:(npars - 1)]
        c(2 * (maxloglik - ll(theta)) - quant,
          sc(theta) - pars[npars] * on_gradient(theta, ...))
    }
    if (use_jacobian) {
        jacobian <- function(pars) {
            theta <- pars[1:(npars - 1)]
            lambda <- pars[npars]
            on_grad <- on_gradient(theta, ...)
            on_hess <- on_hessian(theta, ...)
            info_theta <- info(theta)
            if (!is.numeric(on_hess) ||
                !identical(dim(on_hess), c(length(mle), length(mle))) ||
                any(!is.finite(on_hess))) {
                stop("`on_hessian` must return a finite numeric matrix with one row and one column per parameter.")
            }
            if (!is.numeric(info_theta) ||
                !identical(dim(info_theta), c(length(mle), length(mle))) ||
                any(!is.finite(info_theta))) {
                stop("`information` must return a finite numeric matrix with one row and one column per parameter.")
            }
            out <- matrix(0, npars, npars)
            out[1, 1:(npars - 1)] <- -2 * sc(theta)
            out[2:npars, 1:(npars - 1)] <- -info_theta - lambda * on_hess
            out[2:npars, npars] <- -on_grad
            out
        }
    } else {
        jacobian <- NULL
    }
    if (is.null(start)) {
        if (!do_checks)
            on_grad <- on_gradient(mle, ...)
        info_mle <- info(mle)
        if (!is.numeric(info_mle) || !identical(dim(info_mle), c(length(mle), length(mle))) ||
            any(!is.finite(info_mle)))
            stop("`information` must return a finite numeric matrix with one row and one column per parameter.")
        iinfo <- try(solve(info_mle), silent = TRUE)
        if (inherits(iinfo, "try-error"))
            stop("Could not invert the information matrix at `mle`.")
        step <- drop(iinfo %*% on_grad)
        var_on <- sum(on_grad * step)
        if (!is.finite(var_on) || var_on <= 0)
            stop("Could not construct profile endpoint starting values; ",
                 "the delta-method variance of `on` at `mle` is non-positive or non-finite.")
        lam <- sqrt(quant / var_on)
        trans <- lam * step
        start <- list(lower = c(mle - trans,  lam),
                      upper = c(mle + trans, -lam))
    }
    endpoints <- lapply(c("lower", "upper"), function(side) {
        args <- c(list(x = start[[side]], fn = eq), nleqslv_args)
        args$jac <- jacobian
        do.call(nleqslv, args)
    })
    ci <- sapply(endpoints, function(end) on(end$x[1:(npars - 1)], ...))
    names(ci) <- c("lower", "upper")
    attr(ci, "max|fvec|") <- c(lower = max(abs(endpoints[[1]]$fvec)),
                               upper = max(abs(endpoints[[2]]$fvec)))
    attr(ci, "messages") <- c(lower = endpoints[[1]]$message,
                              upper = endpoints[[2]]$message)
    attr(ci, "loglik") <- c(lower = ll(endpoints[[1]]$x[1:(npars - 1)]),
                            upper = ll(endpoints[[2]]$x[1:(npars - 1)]))
    attr(ci, "solution") <- list(lower = endpoints[[1]]$x,
                                 upper = endpoints[[2]]$x)
    attr(ci, "iter") <- c(lower = endpoints[[1]]$iter,
                          upper = endpoints[[2]]$iter)
    attr(ci, "type") <- "profile"
    ci
}


#' Profile log-likelihood for a focus object
#'
#' Compute a profile log-likelihood for the scalar parameter defined by a focus
#' object.
#'
#' @param fitted An object of class `"focus_list_glm"`, as returned by
#'     [focus()].
#' @param grid_size Number of profile points computed on each side of the
#'     maximum likelihood estimate. Must be at least 2. The outermost point
#'     provides one grid step of plotting space beyond `max_level`.
#' @param max_level Confidence level in `(0, 1)` whose likelihood roots are
#'     placed one grid step inside the outer boundary of the computed profile.
#' @param nleqslv_args List of additional arguments passed to
#'     [nleqslv::nleqslv()] through [profile_ci()].
#' @param ... Currently unused.
#'
#' @details
#' If the fitted model stored in `fitted` is not an ML fit, it is refitted by
#' maximum likelihood. Profile points are computed outwards from the ML
#' estimate on an equally spaced signed likelihood-root grid. Successive
#' points use secant extrapolation for starting values. If a solve does not
#' attain sufficiently small endpoint-equation residuals, it is retried from
#' the preceding solution and then from Wald-type starting values.
#'
#' The result is likelihood-based and does not use a bias-corrected focus
#' estimate as the centre of the profile. The model-specific likelihood
#' quantities are passed to [profile_focus()] to construct the profile.
#'
#' @return
#' `profile.focus_list_glm()` returns a data frame with columns `psi`,
#' `loglik`, and `signed`, and class `"profile_focus_list"`. The columns contain
#' the focus parameter, profile log-likelihood, and signed likelihood root,
#' respectively. The `"max_loglik"` attribute is the unrestricted maximum
#' log-likelihood, `"mle"` is the focus evaluated at the unrestricted MLE, and
#' `"max_level"` is the supplied maximum confidence level. The matrix-valued
#' column `theta` contains the parameter vectors and `lagrange` contains the
#' corresponding Lagrange multipliers.
#'
#' @examples
#' warp_fit <- glm(breaks ~ wool + tension, family = poisson,
#'                 data = warpbreaks)
#' warp_focus <- focus(warp_fit, correction = "no")
#' warp_profile <- profile(warp_focus, grid_size = 10, max_level = 0.99)
#' plot(warp_profile)
#' plot(warp_profile, signed = TRUE)
#'
#' @seealso [profile_focus()], [profile_ci()], [confint.focus_list()]
#'
#' @export
profile.focus_list_glm <- function(fitted,
                                   grid_size = 20,
                                   max_level = 0.9999,
                                   nleqslv_args = list(),
                                   ...) {
    fit <- fitted$object
    if (!identical(fit$type, "ML")) {
        fit <- update(fit, type = "ML", start = coef(fit, model = "mean"))
    }
    p_mean <- length(coef(fit, model = "mean"))
    theta <- coef(fit, model = "full")
    if (fit$family$family %in% c("poisson", "binomial")) {
        theta <- theta[names(coef(fit, model = "mean"))]
    }
    aux <- enrichwith::get_auxiliary_functions(fit)
    m_inds <- seq_len(p_mean)
    d_ind <- p_mean + 1
    split_theta <- function(theta) {
        if (length(theta) == p_mean) {
            list(coefficients = theta)
        } else {
            list(coefficients = theta[m_inds],
                 dispersion = theta[d_ind])
        }
    }
    loglik <- function(theta)
        sum(do.call(aux$dmodel, c(split_theta(theta), list(log = TRUE))))
    score <- function(theta)
        do.call(aux$score, split_theta(theta))
    information <- function(theta)
        do.call(aux$information, split_theta(theta))
    do.call(profile_focus,
            c(list(loglik = loglik,
                   score = score,
                   information = information,
                   mle = theta,
                   on = fitted$on$on,
                   on_gradient = fitted$on$on_gradient,
                   on_hessian = fitted$on$on_hessian,
                   grid_size = grid_size,
                   max_level = max_level,
                   nleqslv_args = nleqslv_args),
              fitted$dots))
}

#' Profile log-likelihood for a scalar focus parameter
#'
#' Compute the profile log-likelihood for a scalar function of a model
#' parameter vector from supplied likelihood quantities.
#'
#' @inheritParams profile_ci
#' @param grid_size Number of profile points computed on each side of the
#'     maximum likelihood estimate. Must be at least 2. The outermost point
#'     provides one grid step of plotting space beyond `max_level`.
#' @param max_level Confidence level in `(0, 1)` whose likelihood roots are
#'     placed one grid step inside the outer boundary of the computed profile.
#'
#' @details
#' Profile points are computed outwards from the maximum likelihood estimate on
#' an equally spaced signed likelihood-root grid. Each pair of points is
#' obtained with [profile_ci()]. Successive pairs use secant extrapolation from
#' preceding solutions as starting values. If a solve does not attain
#' sufficiently small endpoint-equation residuals, it is retried from the
#' preceding solution and then from the Wald-type starting values constructed
#' by [profile_ci()].
#'
#' The grid extends one step beyond `max_level`, allowing confidence limits at
#' `max_level` to be displayed within the plotting range.
#'
#' @return
#' A data frame with columns `psi`, `loglik`, and `signed`, and class
#' `"profile_focus_list"`. The columns contain the focus parameter, profile
#' log-likelihood, and signed likelihood root, respectively. The
#' `"max_loglik"` attribute is the unrestricted maximum log-likelihood,
#' `"mle"` is the focus evaluated at the unrestricted maximum likelihood
#' estimate, and `"max_level"` is the supplied maximum confidence level. The
#' matrix-valued column `theta` contains the parameter vectors and `lagrange`
#' contains the corresponding Lagrange multipliers.
#'
#' @examples
#' ## Normal model parameterized by theta = (mu, log(sigma)).
#' ## The logarithm ensures that sigma = exp(theta[2]) is positive.
#' set.seed(1)
#' y <- rnorm(20, mean = 2)
#' mle <- c(mu = mean(y), log_sigma = 0.5 * log(mean((y - mean(y))^2)))
#' loglik <- function(theta)
#'     sum(dnorm(y, mean = theta[1], sd = exp(theta[2]), log = TRUE))
#'
#' ## Profile the standardized mean mu / sigma.
#' coef_var <- function(theta)
#'     exp(theta[2]) / theta[1]
#'
#' prof <- profile_focus(mle = mle,
#'                       loglik = loglik,
#'                       on = coef_var)
#' plot(prof)
#' plot(prof, signed = TRUE)
#'
#' @seealso [profile_ci()], [profile.focus_list_glm()],
#'     [plot.profile_focus_list()]
#'
#' @export
profile_focus <- function(loglik,
                           score = NULL,
                           information = NULL,
                           mle,
                           on = function(theta) theta[1],
                           on_gradient = NULL,
                           on_hessian = NULL,
                           likelihood_args = list(),
                           nleqslv_args = list(),
                           grid_size = 20,
                           max_level = 0.9999,
                           ...) {
    dots <- list(...)
    r_target <- qnorm(0.5 + max_level / 2)
    r_step <- r_target / (grid_size - 1)
    r_grid <- seq(r_step, r_target + r_step, by = r_step)
    q_grid <- pchisq(r_grid^2, df = 1)
    on_left <- on_right <- ll <- numeric(length = grid_size)
    p <- length(mle)
    theta_left <- theta_right <-
        matrix(NA_real_, grid_size, p,
               dimnames = list(NULL, names(mle)))
    lagrange_left <- lagrange_right <- numeric(grid_size)
    previous <- current <- NULL
    pro <- function(start) {
        do.call(profile_ci,
                c(list(loglik = loglik,
                       score = score,
                       information = information,
                       mle = mle,
                       on = on,
                       on_gradient = on_gradient,
                       on_hessian = on_hessian,
                       likelihood_args = likelihood_args,
                       level = q_grid[j],
                       nleqslv_args = nleqslv_args,
                       start = start,
                       do_checks = (j == 1)), dots))
    }
    is_converged <- function(object, tolerance = 1e-6) {
        residuals <- attr(object, "max|fvec|")
        all(is.finite(object)) && all(is.finite(residuals)) && max(residuals) <= tolerance
    }
    for (j in seq_len(grid_size)) {
        if (j == 1) {
            start <- NULL
        } else if (j == 2) {
            start <- current
        } else {
            start <- list(lower = 2 * current[["lower"]] - previous[["lower"]],
                          upper = 2 * current[["upper"]] - previous[["upper"]])
        }
        obj <- pro(start)
        ## The secant start is used only from the third point onward.
        if (!is_converged(obj) && j > 2)
            obj <- pro(current)
        ## Retry independently from the Wald starts.
        if (!is_converged(obj) && !is.null(start))
            obj <- pro(NULL)
        ## Fail
        if (!is_converged(obj))
            stop("Could not compute the profile at grid point ", j,
                 " (signed likelihood root = ", format(r_grid[j]), "). ",
                 paste(attr(obj, "messages"), collapse = "; "))
        previous <- current
        current <- attr(obj, "solution")
        theta_left[j, ] <- current$lower[seq_len(p)]
        theta_right[j, ] <- current$upper[seq_len(p)]
        lagrange_left[j] <- current$lower[p + 1L]
        lagrange_right[j] <- current$upper[p + 1L]
        on_left[j] <- obj["lower"]
        on_right[j] <- obj["upper"]
        ll[j] <- attr(obj, "loglik")[1]
    }
    max_loglik <- do.call(loglik, c(list(mle), likelihood_args))
    on_mle <- do.call(on, c(list(mle), dots))
    theta_profile <- rbind(theta_left[rev(seq_len(grid_size)), , drop = FALSE],
                           as.numeric(mle),
                           theta_right)
    colnames(theta_profile) <- names(mle)
    rownames(theta_profile) <- NULL
    out <- data.frame(psi = c(rev(on_left), on_mle, on_right),
                      loglik = c(rev(ll), max_loglik, ll),
                      signed = c(-rev(r_grid), 0, r_grid),
                      theta = I(theta_profile),
                      lagrange = c(rev(lagrange_left), 0, lagrange_right))
    class(out) <- c("profile_focus_list", class(out))
    attr(out, "max_loglik") <- max_loglik
    attr(out, "mle") <- on_mle
    attr(out, "max_level") <- max_level
    out
}


#' Print a focus profile log-likelihood
#'
#' Print the main characteristics of a focus profile log-likelihood.
#'
#' @param x An object of class `"profile_focus_list"`, as returned by
#'     [profile_focus()] or [profile.focus_list_glm()].
#' @param digits Number of significant digits used for printing.
#' @param ... Currently unused.
#'
#' @return `x`, invisibly.
#'
#' @seealso [profile_focus()], [profile.focus_list_glm()],
#'     [plot.profile_focus_list()]
#'
#' @export
print.profile_focus_list <- function(x,
                                     digits = max(3L, getOption("digits") - 2L),
                                     ...) {
    signed_range <- range(x$signed)
    positive <- sort(x$signed[x$signed > 0])
    spacing <- median(diff(c(0, positive)))
    boundary_level <- pchisq(max(abs(x$signed))^2, df = 1)
    format_value <- function(value)
        format(signif(value, digits), trim = TRUE)

    cat("Profile log-likelihood for a scalar focus\n\n")
    cat("Focus at MLE:", format_value(attr(x, "mle")), "\n")
    cat("Maximum log-likelihood:",
        format_value(attr(x, "max_loglik")), "\n")
    cat("Parameter dimension:", ncol(x$theta), "\n")
    cat("Profile points:", nrow(x), "\n")
    cat("Points per branch:",
        sum(x$signed < 0), "left,", sum(x$signed > 0), "right\n")
    cat("Signed-root range:",
        format_value(signed_range[1]), "to",
        format_value(signed_range[2]), "\n")
    cat("Signed-root spacing:", format_value(spacing), "\n")
    cat("Requested maximum level:",
        format_value(attr(x, "max_level")), "\n")
    cat("Boundary confidence level:",
        format_value(boundary_level), "\n")
    invisible(x)
}


#' Plot a focus profile log-likelihood
#'
#' Plot a focus profile log-likelihood either on the log-likelihood scale or the
#' signed likelihood-root scale.
#'
#' @param x An object of class `"profile_focus_list"`, as returned by
#'     [profile_focus()] or [profile.focus_list_glm()].
#' @param level Confidence level used to draw the horizontal cutoff and
#'     vertical confidence limits.
#' @param signed Logical. If `TRUE`, plot the signed likelihood root; otherwise
#'     plot the profile log-likelihood.
#' @param ... Additional graphical arguments passed to [graphics::plot.default()].
#'
#' @details `level` must lie within the range covered by the computed profile.
#'     If it does not, recompute the profile with a larger `max_level`.
#'
#' @return Called for its side effect of drawing a plot.
#'
#' @seealso [profile_focus()], [profile.focus_list_glm()], [profile_ci()]
#'
#' @export
plot.profile_focus_list <- function(x, level = 0.95, signed = FALSE, ...) {
    max_loglik <- attr(x, "max_loglik")
    qua <- qnorm(0.5 + level/2)
    if (qua > max(abs(x$signed))) {
        stop("`level` exceeds the range of the supplied profile; ",
             "recompute the profile with a larger `max_level`.")
    }
    fn <- splinefun(x = x$signed, y = x$psi, method = "monoH.FC")
    r <- seq(min(x$signed), max(x$signed), length.out = 201)
    psi <- fn(r)
    ci <- c(fn(-qua), fn(qua))
    if (signed) {
        plot.default(psi, r, type = "l", xlab = expression(psi), ylab = "Signed likelihood root", ...)
        abline(h = c(-qua, qua), lty = 3, col = "lightgray")
        points(attr(x, "mle"), 0, pch = 21, bg = "lightgray")
    } else {
        plot(psi, max_loglik - r^2/2, type = "l", xlab = expression(psi), ylab = "Log-likelihood", ...)
        cutoff <- max_loglik - qchisq(level, 1) / 2
        abline(h = cutoff, lty = 3, col = "lightgray")
        points(attr(x, "mle"), max_loglik, pch = 21, bg = "lightgray")
    }
    abline(v = ci, lty = 1, col = "lightgray")
}
