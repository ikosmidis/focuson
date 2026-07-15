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
#'     with `information`, an analytic Jacobian is passed to
#'     [nleqslv::nleqslv()].
#' @param likelihood_args List of additional arguments passed to
#'     `loglik`, `score`, and `information`.
#' @param nleqslv_args List of additional arguments passed to
#'     [nleqslv::nleqslv()].
#' @param level Confidence level. Default is `0.95`.
#' @param ... Additional arguments passed to `on` and `on_gradient`.
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
#' The function solves these equations twice, from Wald-type starting
#' values in opposite directions, and returns the corresponding values
#' of `on`. This is a low-level routine: endpoint convergence
#' diagnostics are returned as attributes, and callers can decide how
#' strictly to enforce them.
#'
#' If `information` and `on_hessian` are both supplied, `profile_ci()`
#' passes the analytic Jacobian of the endpoint equations to
#' [nleqslv::nleqslv()]. If either is omitted, [nleqslv::nleqslv()]
#' computes its own numerical Jacobian.  Supplying analytic `score`,
#' `information`, `on_gradient`, and `on_hessian` can substantially
#' reduce computation.
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
#'            level = 0.95, on = me, ldose = 2)
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
                       ...) {
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
    ll <- function(theta) do.call(loglik, c(list(theta), likelihood_args))
    if (is.null(on_gradient))
        on_gradient  <- function(theta, ...) numDeriv::grad(on, theta, ...)
    if (is.null(score)) {
        sc <- function(theta) numDeriv::grad(ll, theta)
    } else {
        sc <- function(theta) do.call(score, c(list(theta), likelihood_args))
    }
    if (is.null(information)) {
        info <- function(theta) -numDeriv::hessian(ll, theta)
    } else {
        info <- function(theta) do.call(information, c(list(theta), likelihood_args))
    }
    use_jacobian <- !is.null(information) && !is.null(on_hessian)
    maxloglik <- ll(mle)
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
    quant <- qchisq(level, 1)
    npars <- length(mle) + 1
    eq <- function(pars) {
        theta <- pars[1:(npars - 1)]
        l <- ll(theta)
        s <- sc(theta)
        c(2 * (maxloglik - l) - quant, s - pars[npars] * on_gradient(theta, ...))
    }
    if (use_jacobian) {
        jacobian <- function(pars) {
            theta <- pars[1:(npars - 1)]
            lambda <- pars[npars]
            s <- sc(theta)
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
            out[1, 1:(npars - 1)] <- -2 * s
            out[2:npars, 1:(npars - 1)] <- -info_theta - lambda * on_hess
            out[2:npars, npars] <- -on_grad
            out
        }
    } else {
        jacobian <- NULL
    }
    ## Starting values
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
    endpoints <- lapply(c(-1, 1), function(s) {
        args <- c(list(x = c(mle + s * trans, - s * lam), fn = eq), nleqslv_args)
        if (!is.null(jacobian)) {
            args$jac <- jacobian
        }
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
    attr(ci, "type") <- "profile"
    ci
}
