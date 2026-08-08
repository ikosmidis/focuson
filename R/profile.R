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
#'     leave it `NULL` (default) so that separate Wald-type starting values are
#'     constructed for the lower and upper endpoints. Inappropriate or
#'     identical starting values may cause both endpoint calculations to
#'     converge to the same solution.
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
#'   \item{`"type"`}{The string `"pl"`, denoting profile-likelihood
#'     inference.}
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
#' Irrespective of `do_checks`, the returned lower and upper endpoints are
#' required to lie below and above, respectively, the value of `on` at `mle`.
#' If they do not, the function cannot confirm that the solutions are the two
#' endpoints surrounding the focus estimate. This can occur when the
#' confidence set is unbounded or disconnected, or when both endpoint
#' calculations converge to the same or otherwise unintended solutions. The
#' function returns an error rather than two potentially misleading endpoints.
#'
#' `profile_ci()` targets the two finite endpoints of a connected confidence
#' interval containing the value of `on` at `mle`. It does not search for
#' additional disconnected components of a confidence set. If the relevant
#' confidence set is unbounded, the function is not intended to return an
#' infinite endpoint, and the endpoint calculation should be regarded as
#' unsuccessful. Because the endpoint equations are stationarity conditions
#' and can have more than one solution, the function also does not provide a
#' global guarantee that the solutions found are the outermost likelihood
#' crossings. Starting values can therefore matter in non-regular problems.
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
    on_mle <- on(mle, ...)
    opposite <- length(on_mle) == 1L && is.finite(on_mle) &&
        isTRUE(ci["lower"] < on_mle) &&
        isTRUE(ci["upper"] > on_mle)
    if (!opposite) {
        stop("Could not identify profile endpoints on opposite sides of ",
             "the focus estimate. The confidence set may be unbounded or ",
             "disconnected, or the endpoint equations may have converged ",
             "to unintended solutions.")
    }
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
    attr(ci, "type") <- "pl"
    ci
}


#' Profile log-likelihood for a focus object
#'
#' Compute a profile log-likelihood for the scalar parameter defined by a focus
#' object.
#'
#' @inheritParams profile_focus
#' @param fitted An object of class `"focus_list_glm"`, as returned by
#'     [focus()].
#' @param ... Currently unused.
#'
#' @details
#' If the fitted model stored in `fitted` is not an ML fit, it is refitted by
#' maximum likelihood. The model-specific likelihood quantities are passed to
#' [profile_focus()] to construct the requested profile.
#'
#' The result is likelihood-based and does not use a bias-corrected focus
#' estimate as the centre of the profile. The signed likelihood root uses the
#' sign of the focus at the unrestricted MLE minus the focus value at the
#' profile point.
#'
#' With `approach = "focus_grid"`, a range that does not contain the focus at
#' the MLE produces a single profile branch. See [profile_focus()] for details.
#'
#' @return
#' `profile.focus_list_glm()` returns a data frame with columns `psi`,
#' `loglik`, `theta`, and `lagrange`, and class `"profile_focus_list"`.
#' The first two columns contain the focus parameter and profile
#' log-likelihood. The matrix-valued column `theta` contains the parameter
#' vectors and `lagrange` contains the corresponding Lagrange multipliers.
#' The `"max_loglik"` attribute is the unrestricted maximum log-likelihood,
#' `"mle"` is the focus evaluated at the unrestricted MLE, and `"approach"`
#' records the profiling approach.
#'
#' @examples
#' warp_fit <- glm(breaks ~ wool + tension, family = poisson,
#'                 data = warpbreaks)
#' warp_focus <- focus(warp_fit, correction = "no")
#' warp_profile <- profile(warp_focus, grid_size = 10, max_level = 0.99)
#' plot(warp_profile)
#' plot(warp_profile, signed = TRUE)
#'
#' \dontrun{
#' ## Exploring the profile log-likelihood
#' data("endometrial", package = "brglm2")
#' endo <- glm(HG ~ NV + PI + EH,
#'             data = endometrial,
#'             family = binomial("logit"))
#'
#' ## Focus on the coefficient of NV, using maximum likelihood
#' (endo2 <- focus(endo, on = function(theta) theta[2], correction = "no"))
#'
#' ## The ML estimate for NV is large in absolute value, and the
#' ## default method for profiling the log-likelihood (with
#' ## `approach = "VM"`) is numerically unstable
#' try(prof <- profile(endo2))
#'
#' ## A closer inspection of the profile log-likelihood reveals that #
#' ## the profile is essentially monotone; the log-likelihood approaches
#' ## its supremum as the NV coefficient tends to infinity.
#' prof <- profile(endo2, approach = "focus_grid", focus_range = c(0, 30))
#' plot(prof)
#'
#' ## which can be verified by
#' library("detectseparation")
#' update(endo, method = "detect_separation")
#'
#' }
#' @seealso [profile_focus()], [profile_ci()], [confint.focus_list()]
#'
#' @export
profile.focus_list_glm <- function(fitted,
                                   grid_size = 20,
                                   max_level = 0.995,
                                   nleqslv_args = list(),
                                   approach = c("VM", "focus_grid"),
                                   focus_range = NULL,
                                   auglag_args = list(),
                                   ...) {
    components <- .profile_glm_components(fitted)
    do.call(profile_focus,
            c(components,
              list(on = fitted$on$on,
                   on_gradient = fitted$on$on_gradient,
                   on_hessian = fitted$on$on_hessian,
                   grid_size = grid_size,
                   max_level = max_level,
                   approach = approach,
                   focus_range = focus_range,
                   nleqslv_args = nleqslv_args,
                   auglag_args = auglag_args),
              fitted$dots))
}

#' Profile log-likelihood for a scalar focus parameter
#'
#' Compute the profile log-likelihood for a scalar function of a model
#' parameter vector from supplied likelihood quantities.
#'
#' @inheritParams profile_ci
#' @param grid_size Number of profile points computed on each side of the
#'     maximum likelihood estimate. Must be at least 2. For the VM approach,
#'     the outermost point provides one grid step of plotting space beyond
#'     `max_level`. For a focus-grid profile whose range does not contain the
#'     focus at the MLE, `grid_size + 1` points are computed across the range.
#' @param max_level Nominal likelihood-ratio level in `(0, 1)` whose
#'     likelihood roots are placed one grid step inside the outer boundary of
#'     the computed profile. Used only when `approach = "VM"`.
#' @param approach Character. `"VM"` uses the endpoint equations of Venzon and
#'     Moolgavkar (1988) on a likelihood-root grid. `"focus_grid"` maximizes the
#'     likelihood at fixed focus values.
#' @param focus_range Optional numeric vector containing the lower and upper
#'     focus values. It is required when `approach = "focus_grid"` and must
#'     contain two distinct finite values. It is ignored otherwise.
#' @param auglag_args Named list of additional arguments passed to
#'     [alabama::auglag()] when `approach = "focus_grid"`. The arguments `par`,
#'     `fn`, `gr`, `hin`, `hin.jac`, `heq`, and `heq.jac` are determined
#'     internally. The argument is ignored when `approach = "VM"`.
#'
#' @details
#' With `approach = "VM"`, profile points are computed outwards from the
#' maximum likelihood estimate on an equally spaced signed likelihood-root
#' grid. Each pair of points is obtained with [profile_ci()]. Successive pairs
#' use secant extrapolation from preceding solutions as starting values.
#'
#' With `approach = "focus_grid"`, if the focus at the MLE lies within
#' `focus_range`, `grid_size` equally spaced constrained fits are computed
#' between it and each available endpoint. If it lies outside `focus_range`,
#' `grid_size + 1` fits are computed across the range, starting at the nearest
#' endpoint and proceeding towards the other endpoint. Successive fits are
#' warm-started from the preceding fit. The `"mle"` and `"max_loglik"`
#' attributes always refer to the supplied unrestricted MLE, which need not be
#' among the returned profile points. Only `loglik`, `mle`, and `on` are
#' required; missing derivatives are evaluated numerically by
#' [alabama::auglag()].
#'
#' `nleqslv_args` applies only to `approach = "VM"`, while `auglag_args`
#' applies only to `approach = "focus_grid"`. The inactive argument is ignored.
#'
#' The VM grid extends one step beyond `max_level`, allowing confidence limits
#' at `max_level` to be displayed within the plotting range.
#' The signed likelihood root is defined as
#' \deqn{\mathop{\rm sign}(\hat\psi-\psi)
#' \left[2\{\ell(\hat\theta)-\ell_p(\psi)\}\right]^{1/2}.}
#'
#' @return
#' A data frame with columns `psi`, `loglik`, `theta`, and `lagrange`, and class
#' `"profile_focus_list"`. The first two columns contain the focus parameter and
#' profile log-likelihood. The matrix-valued column `theta` contains the
#' parameter vectors and `lagrange` contains the corresponding Lagrange
#' multipliers. The `"max_loglik"` attribute is the unrestricted maximum
#' log-likelihood, `"mle"` is the focus evaluated at the unrestricted maximum
#' likelihood estimate, and `"approach"` records the profiling approach. For
#' VM profiles, `"max_level"` is the supplied nominal likelihood-ratio level;
#' for focus-grid profiles, `"focus_range"` is the supplied focus range.
#'
#' @examples
#' ## Normal model parameterized by theta = (mu, log(sigma)).
#' set.seed(1)
#' y <- rnorm(20, mean = 0.5)
#' mle <- c(mu = mean(y), log_sigma = 0.5 * log(mean((y - mean(y))^2)))
#' loglik <- function(theta)
#'     sum(dnorm(y, mean = theta[1], sd = exp(theta[2]), log = TRUE))
#'
#' ## Profile the coefficient of variation sigma / mu.
#' coef_var <- function(theta)
#'     exp(theta[2]) / theta[1]
#'
#' prof <- profile_focus(mle = mle,
#'                       loglik = loglik,
#'                       on = coef_var,
#'                       max_level = 0.99)
#' plot(prof)
#' plot(prof, interpolation = "cubic", level = 0.9, ci = TRUE)
#' plot(prof, signed = TRUE, level = 0.9, ci = TRUE)
#'
#' ## Compare the interpolated limits with directly computed endpoints.
#' confint(prof, level = 0.9, interpolation = "cubic")
#' profile_ci(loglik, on = coef_var, mle = mle, level = 0.9)
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
                          max_level = 0.995,
                          approach = c("VM", "focus_grid"),
                          focus_range = NULL,
                          auglag_args = list(),
                          ...) {
    approach <- match.arg(approach)
    dots <- list(...)
    required <- list(loglik = loglik, on = on)
    optional <- list(score = score,
                     information = information,
                     on_gradient = on_gradient,
                     on_hessian = on_hessian)
    check_required <- !vapply(required, is.function, logical(1))
    if (any(check_required))
        stop("The following arguments must be functions: ",
             paste(sprintf("`%s`", names(required)[check_required]),
                   collapse = ", "), ".")
    check_optional <- !vapply(optional,
                             function(x) is.null(x) || is.function(x),
                             logical(1))
    if (any(check_optional))
        stop("The following arguments must be functions or `NULL`: ",
             paste(sprintf("`%s`", names(optional)[check_optional]),
                   collapse = ", "), ".")
    if (!is.numeric(mle) || length(mle) == 0L || any(!is.finite(mle)))
        stop("`mle` must be a finite numeric vector.")
    if (!is.list(likelihood_args))
        stop("`likelihood_args` must be a list.")
    if (!is.numeric(grid_size) || length(grid_size) != 1L ||
        !is.finite(grid_size) || grid_size < 2L ||
        grid_size != as.integer(grid_size))
        stop("`grid_size` must be an integer of at least 2.")

    ll <- function(theta)
        do.call(loglik, c(list(theta), likelihood_args))
    sc <- if (is.null(score)) NULL else
        function(theta) do.call(score, c(list(theta), likelihood_args))
    on_value <- function(theta)
        do.call(on, c(list(theta), dots))
    on_grad <- if (is.null(on_gradient)) NULL else
        function(theta) do.call(on_gradient, c(list(theta), dots))
    max_loglik <- ll(mle)
    on_mle <- on_value(mle)
    if (!is.numeric(max_loglik) || length(max_loglik) != 1L ||
        !is.finite(max_loglik))
        stop("`loglik` must return a finite numeric scalar at `mle`.")
    if (!is.numeric(on_mle) || length(on_mle) != 1L || !is.finite(on_mle))
        stop("`on` must return a finite numeric scalar at `mle`.")

    p <- length(mle)
    empty_theta <- matrix(numeric(0), nrow = 0L, ncol = p,
                          dimnames = list(NULL, names(mle)))
    on_center <- on_mle
    ll_center <- max_loglik
    theta_center <- matrix(as.numeric(mle), nrow = 1L)
    lagrange_center <- 0
    if (approach == "VM") {
        if (!is.list(nleqslv_args))
            stop("`nleqslv_args` must be a list.")
        if (!is.numeric(max_level) || length(max_level) != 1L ||
            !is.finite(max_level) || max_level <= 0 || max_level >= 1)
            stop("`max_level` must be a number in (0, 1).")

        r_target <- qnorm(0.5 + max_level / 2)
        r_step <- r_target / (grid_size - 1)
        root_grid <- seq(r_step, r_target + r_step, by = r_step)
        q_grid <- pchisq(root_grid^2, df = 1)
        on_left <- on_right <- ll_left <- numeric(grid_size)
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
            all(is.finite(object)) && all(is.finite(residuals)) &&
                max(residuals) <= tolerance
        }
        for (j in seq_len(grid_size)) {
            if (j == 1) {
                start <- NULL
            } else if (j == 2) {
                start <- current
            } else {
                start <- list(
                    lower = 2 * current[["lower"]] - previous[["lower"]],
                    upper = 2 * current[["upper"]] - previous[["upper"]]
                )
            }
            obj <- pro(start)
            if (!is_converged(obj) && j > 2)
                obj <- pro(current)
            if (!is_converged(obj) && !is.null(start))
                obj <- pro(NULL)
            if (!is_converged(obj))
                stop("Could not compute the profile at grid point ", j,
                     " (likelihood-root magnitude = ",
                     format(root_grid[j]), "). ",
                     paste(attr(obj, "messages"), collapse = "; "))
            previous <- current
            current <- attr(obj, "solution")
            theta_left[j, ] <- current$lower[seq_len(p)]
            theta_right[j, ] <- current$upper[seq_len(p)]
            lagrange_left[j] <- current$lower[p + 1L]
            lagrange_right[j] <- current$upper[p + 1L]
            on_left[j] <- obj["lower"]
            on_right[j] <- obj["upper"]
            ll_left[j] <- attr(obj, "loglik")[1]
        }
        ll_right <- ll_left
    } else {
        if (!is.list(auglag_args))
            stop("`auglag_args` must be a list.")
        if (!is.numeric(focus_range) || length(focus_range) != 2L ||
            any(!is.finite(focus_range)))
            stop("`focus_range` must contain two finite numeric values.")
        focus_range <- sort(focus_range)
        if (focus_range[1L] == focus_range[2L])
            stop("`focus_range` must contain two distinct values.")
        control_outer <- list(lam0 = 0,
                              trace = FALSE,
                              kkt2.check = FALSE)
        if ("control.outer" %in% names(auglag_args)) {
            if (!is.null(auglag_args$control.outer))
                control_outer <- modifyList(control_outer,
                                            auglag_args$control.outer)
            auglag_args$control.outer <- NULL
        }
        constraint_tolerance <- if (is.null(control_outer$eps))
            1e-6
        else
            max(1e-6, control_outer$eps)
        equality_scale <- if (is.null(control_outer$e.scale))
            1
        else
            control_outer$e.scale
        solve_at <- function(psi, start) {
            args <- list(par = start,
                         fn = function(theta) -ll(theta),
                         heq = function(theta) on_value(theta) - psi,
                         control.outer = control_outer)
            args <- c(args, auglag_args)
            if (!is.null(sc))
                args$gr <- function(theta)
                    -sc(theta)
            if (!is.null(on_grad))
                args$heq.jac <- function(theta)
                    matrix(on_grad(theta), nrow = 1L)
            fit <- do.call(alabama::auglag, args)
            fit_loglik <- ll(fit$par)
            constraint_error <- abs(on_value(fit$par) - psi)
            if (!isTRUE(fit$convergence == 0L) ||
                any(!is.finite(fit$par)) ||
                !is.finite(fit_loglik) ||
                !is.finite(constraint_error) ||
                constraint_error > constraint_tolerance)
                stop("Could not compute the constrained profile at focus value",
                     format(psi), ".")
            list(psi = psi,
                 loglik = fit_loglik,
                 theta = fit$par,
                 lagrange = - fit$lambda[1L] / equality_scale[1L])
        }
        solve_side <- function(values) {
            out <- vector("list", length(values))
            start <- mle
            for (j in seq_along(values)) {
                out[[j]] <- solve_at(values[j], start)
                start <- out[[j]]$theta
            }
            out
        }
        left <- right <- list()
        if (on_mle < focus_range[1L]) {
            focus_right <- seq(focus_range[1L], focus_range[2L],
                               length.out = grid_size + 1L)
            right <- solve_side(focus_right)
            on_center <- ll_center <- lagrange_center <- numeric(0)
            theta_center <- empty_theta
        } else if (on_mle > focus_range[2L]) {
            focus_left <- seq(focus_range[2L], focus_range[1L],
                              length.out = grid_size + 1L)
            left <- solve_side(focus_left)
            on_center <- ll_center <- lagrange_center <- numeric(0)
            theta_center <- empty_theta
        } else {
            if (on_mle > focus_range[1L]) {
                focus_left <- seq(on_mle, focus_range[1L],
                                  length.out = grid_size + 1L)[-1L]
                left <- solve_side(focus_left)
            }
            if (on_mle < focus_range[2L]) {
                focus_right <- seq(on_mle, focus_range[2L],
                                   length.out = grid_size + 1L)[-1L]
                right <- solve_side(focus_right)
            }
        }
        on_left <- vapply(left, `[[`, numeric(1), "psi")
        on_right <- vapply(right, `[[`, numeric(1), "psi")
        ll_left <- vapply(left, `[[`, numeric(1), "loglik")
        ll_right <- vapply(right, `[[`, numeric(1), "loglik")
        theta_left <- if (length(left))
            do.call(rbind, lapply(left, `[[`, "theta"))
        else
            empty_theta
        theta_right <- if (length(right))
            do.call(rbind, lapply(right, `[[`, "theta"))
        else
            empty_theta
        lagrange_left <- vapply(left, `[[`, numeric(1), "lagrange")
        lagrange_right <- vapply(right, `[[`, numeric(1), "lagrange")
    }

    theta_profile <- rbind(theta_left[rev(seq_len(nrow(theta_left))), ,
                                      drop = FALSE],
                           theta_center,
                           theta_right)
    dimnames(theta_profile) <- list(NULL, names(mle))
    out <- data.frame(psi = c(rev(on_left), on_center, on_right),
                      loglik = c(rev(ll_left), ll_center, ll_right),
                      theta = I(theta_profile),
                      lagrange = c(rev(lagrange_left), lagrange_center,
                                   lagrange_right))
    class(out) <- c("profile_focus_list", class(out))
    attr(out, "max_loglik") <- max_loglik
    attr(out, "theta_mle") <- mle
    attr(out, "mle") <- on_mle
    attr(out, "approach") <- approach
    attr(out, "max_level") <- max_level
    attr(out, "focus_range") <- focus_range
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
#' @details
#' The printed entries have the following meanings:
#'
#' * **Profiling approach:** the method used to construct the profile, either
#'   the Venzon and Moolgavkar approach (`"VM"`) or constrained optimization
#'   over a grid of focus values (`"focus_grid"`).
#' * **Focus at MLE:** the focus evaluated at the supplied unrestricted maximum
#'   likelihood estimate. For a single-branch focus-grid profile, this value
#'   lies outside the supplied focus range and is not a returned profile point.
#' * **Maximum log-likelihood:** the log-likelihood at the supplied
#'   unrestricted maximum likelihood estimate.
#' * **Parameter dimension:** the dimension of the model parameter vector
#'   stored in the matrix-valued `theta` column.
#' * **Profile points:** the number of rows in the computed profile.
#' * **Points by side:** the number of profile points to the left and right of
#'   the focus at the MLE. The MLE itself, when present, is not counted on
#'   either side.
#' * **Requested nominal level:** for a VM profile, the `max_level` used to
#'   construct the likelihood-root grid. The grid extends one step beyond this
#'   level.
#' * **Focus range:** for a focus-grid profile, the supplied lower and upper
#'   focus values.
#' * **Implied level at grid boundary:** the chi-squared reference level
#'   corresponding to the likelihood drop at the outer endpoint of the
#'   computed profile. When both branches are present, the smaller of their
#'   endpoint levels is reported. For a VM profile this is generally larger
#'   than the requested nominal level because of the extra grid step. For a
#'   single-branch profile, it describes that branch only and does not imply
#'   that two-sided confidence limits are available. It is not an estimate of
#'   the finite-sample coverage of an interval.
#'
#' @return `x`, invisibly.
#'
#' @seealso [profile_focus()], [profile.focus_list_glm()],
#'     [plot.profile_focus_list()], [confint.profile_focus_list()]
#'
#' @export
print.profile_focus_list <- function(x,
                                     digits = max(3L, getOption("digits") - 2L),
                                     ...) {
    focus_mle <- attr(x, "mle")
    max_loglik <- attr(x, "max_loglik")
    left <- x$psi < focus_mle
    right <- x$psi > focus_mle
    boundary_drops <- numeric(0)
    if (any(left))
        boundary_drops <- c(boundary_drops,
                            max_loglik - x$loglik[which.min(x$psi)])
    if (any(right))
        boundary_drops <- c(boundary_drops,
                            max_loglik - x$loglik[which.max(x$psi)])
    if (length(boundary_drops)) {
        boundary_lr_level <- pchisq(2 * max(min(boundary_drops), 0),
                                    df = 1)
    } else {
        boundary_lr_level <- NA_real_
    }
    format_value <- function(value)
        format(signif(value, digits), trim = TRUE)
    cat("Profile log-likelihood for a scalar focus\n\n")
    cat("Profiling approach:", attr(x, "approach"), "\n")
    cat("Focus at MLE:", format_value(attr(x, "mle")), "\n")
    cat("Maximum log-likelihood:",
        format_value(attr(x, "max_loglik")), "\n")
    cat("Parameter dimension:", ncol(x$theta), "\n")
    cat("Profile points:", nrow(x), "\n")
    cat("Points by side:",
        sum(x$psi < focus_mle), "left,",
        sum(x$psi > focus_mle), "right\n")
    if (identical(attr(x, "approach"), "VM")) {
        cat("Requested maximum nominal level:",
            format_value(attr(x, "max_level")), "\n")
    } else {
        cat("Focus range:",
            format_value(attr(x, "focus_range")[1L]), "to",
            format_value(attr(x, "focus_range")[2L]), "\n")
    }
    cat("Implied level at grid boundary:",
        if (is.na(boundary_lr_level)) "unavailable" else
            format_value(boundary_lr_level),
        "\n")
    invisible(x)
}

