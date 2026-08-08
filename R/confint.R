#' Confidence intervals for focus objects
#'
#' Compute confidence intervals for objects returned by [focus()].
#'
#' @param object An object of class `"focus_list"`.
#' @param parm Currently unused.
#' @param level Confidence level.
#' @param method Character string specifying the confidence interval
#'     method. One of `"wald"`, `"pl"`, or `"hulc"`. `"pl"` denotes
#'     profile-likelihood inference.
#' @param se_at Character string specifying where the delta-method
#'     standard error is evaluated for `method = "wald"`. `"supplied"`
#'     uses the standard error stored in `object`; `"compatible"`
#'     tries to evaluate the standard error at a reconstructed model
#'     parameter vector compatible with the reported focus estimate.
#' @param V_function Optional function returning the covariance matrix
#'     at a supplied model parameter vector. Required for
#'     `focus_engine_list` objects when `se_at = "compatible"`.
#' @param se_control A list of control parameters passed to
#'     [focus_se()] when `se_at = "compatible"`.
#' @param nleqslv_args A list of control arguments passed to
#'     [nleqslv::nleqslv()] through [profile_ci()] when
#'     `method = "pl"`.
#' @param ... Additional arguments for the confidence interval
#'     method. For `method = "hulc"`, these are passed to [hulc_ci()],
#'     except that the nominal level is determined by `level`. For
#'     `method = "wald"` and `se_at = "compatible"` with
#'     `focus_engine_list` objects, these are passed to `V_function`.
#'
#' @return
#' A numeric vector of length 2 with names `"lower"` and `"upper"`.
#'
#' @details
#' For `method = "wald"`, the interval is computed from the stored point
#' estimate and a delta-method standard error using the normal approximation.
#' By default, the stored standard error is used. If `se_at = "compatible"`,
#' `focus_se()` is used to compute a compatible standard error lazily. If that
#' computation fails, a warning is issued and the stored standard error is used.
#'
#' For `method = "pl"`, [profile_ci()] is used to compute a likelihood
#' profile interval for the scalar parameter defined by the stored `on`
#' function. If the stored fitted object is not an ML fit, it is refitted by
#' maximum likelihood before profiling. This interval is likelihood-based and
#' does not use the bias-corrected focus estimate as the likelihood centre.
#' Profile intervals are currently not available for `focus_engine()` results.
#'
#' For `method = "hulc"`, [hulc_ci()] is applied to data recovered from the
#' stored fitted object using [focus_statistic()] as the statistic evaluated
#' on each partition. The recovered data must be suitable for refitting the
#' stored fitted object with [stats::update()]. HulC intervals are currently
#' not available for `focus_engine()` results.
#'
#' The nominal coverage level is determined by `level`; users should
#' not supply `level` in `...`.
#'
#' @seealso [focus()], [hulc_ci()], [focus_statistic()]
#'
#' @export
confint.focus_list <- function(object,
                               parm,
                               level = 0.95,
                               method = "wald",
                               se_at = "supplied",
                               V_function = NULL,
                               se_control = list(),
                               nleqslv_args = list(),
                               ...) {
    method <- match.arg(method, c("wald", "pl", "hulc"))
    se_at <- match.arg(se_at, c("supplied", "compatible"))
    if (!is.list(se_control)) {
        stop("`se_control` must be a list.")
    }
    if (!is.list(nleqslv_args)) {
        stop("`nleqslv_args` must be a list.")
    }
    alpha <- 1 - level

    if (identical(method, "wald")) {
        se <- object$se
        if (identical(se_at, "compatible")) {
            se_info <- if (inherits(object, "focus_engine_list")) {
                try(focus_se(object,
                             V_function = V_function,
                             control = se_control,
                             ...),
                    silent = TRUE)
            } else {
                try(focus_se(object, control = se_control),
                    silent = TRUE)
            }
            if (inherits(se_info, "try-error")) {
                warning("Could not compute compatible standard error; using supplied standard error instead. ",
                        "Original error: ", conditionMessage(attr(se_info, "condition")),
                        call. = FALSE)
                se_at <- "supplied"
            } else {
                se <- se_info$se
            }
        }
        ci <- unname(object$estimate) + c(-1, 1) * qnorm(1 - alpha / 2) * se
        names(ci) <- c("lower", "upper")
        attr(ci, "level") <- level
        attr(ci, "type") <- "wald"
        attr(ci, "se_at") <- se_at
        if (exists("se_info", inherits = FALSE) && !inherits(se_info, "try-error")) {
            attr(ci, "se_info") <- se_info
        }
        return(ci)
    }

    if (identical(method, "pl")) {
        if (inherits(object, "focus_engine_list")) {
            stop("`method = \"pl\"` is not available for `focus_engine()` results.")
        }
        return(.confint_profile_focus_list_glm(object = object,
                                               level = level,
                                               nleqslv_args = nleqslv_args))
    }

    if (inherits(object, "focus_engine_list")) {
        stop("`method = \"hulc\"` is not available for `focus_engine()` results.")
    }

    correction <- object$correction
    on_funs <- object$on
    on <- on_funs$on
    on_gradient <- on_funs$on_gradient
    on_hessian <- on_funs$on_hessian
    odots <- object$dots
    statistic <- function(data) {
        do.call(focus_statistic,
                c(list(data = data,
                       object = object$object,
                       on = on,
                       correction = correction,
                       on_gradient = on_gradient,
                       on_hessian = on_hessian),
                  odots))
    }
    do.call(hulc_ci,
            c(list(data = .refit_data(object$object),
                   statistic = statistic,
                   level = level),
              list(...)))
}

.refit_data <- function(object) {
    data_call <- object$call$data
    if (!is.null(data_call)) {
        data <- try(eval(data_call, envir = environment(stats::formula(object))),
                    silent = TRUE)
        if (!inherits(data, "try-error")) {
            return(as.data.frame(data))
        }
    }
    data <- try(model.frame(object), silent = TRUE)
    if (!inherits(data, "try-error")) {
        return(as.data.frame(data))
    }
    stop("Could not recover data suitable for refitting `object`; ",
         "`method = \"hulc\"` requires the original model data.")
}

.confint_profile_focus_list_glm <- function(object, level, nleqslv_args) {
    components <- .profile_glm_components(object)
    do.call(profile_ci,
            c(components,
              list(on = object$on$on,
                   on_gradient = object$on$on_gradient,
                   on_hessian = object$on$on_hessian,
                   level = level,
                   nleqslv_args = nleqslv_args),
              object$dots))
}


.profile_signed <- function(object) {
    sign(attr(object, "mle") - object$psi) *
        sqrt(2 * pmax(attr(object, "max_loglik") - object$loglik, 0))
}

.modified_profile_maximum <- function(object) {
    objective <- splinefun(object$psi, object$modified_loglik)
    focus_range <- range(object$psi)
    fit <- optimize(objective,
                    interval = focus_range,
                    maximum = TRUE,
                    tol = sqrt(.Machine$double.eps))
    boundary_loglik <- objective(focus_range)
    tolerance <- sqrt(.Machine$double.eps) *
        max(abs(c(fit$objective, boundary_loglik)), 1)
    if (max(boundary_loglik) >= fit$objective - tolerance)
        stop("The modified profile maximum is at the grid boundary; ",
             "recompute the profile over a wider focus range.")
    c(psi = fit$maximum, loglik = fit$objective)
}

#' Confidence intervals from a focus profile
#'
#' Extract a confidence interval from a computed focus profile by
#' interpolation on the signed likelihood-root scale.
#'
#' @param object An object of class `"profile_focus_list"`, as returned by
#'     [profile_focus()], [modified_profile_focus()], or
#'     [profile.focus_list_glm()].
#' @param parm Currently unused.
#' @param level Nominal confidence level.
#' @param method Character. The likelihood-based method used to construct the
#'     interval. `"pl"` (default) uses the ordinary profile likelihood,
#'     `"mpl"` uses the modified profile likelihood, and `"rstar"` uses the
#'     modified signed likelihood-ratio statistic \eqn{r^*}.
#' @param interpolation Character. Interpolation method used between computed
#'     profile points. `"linear"` (default) uses [stats::approxfun()] and
#'     `"cubic"` uses [stats::splinefun()].
#' @param ... Currently unused.
#'
#' @details
#' This method interpolates the supplied profile and does not perform further
#' likelihood evaluations. The profile must extend beyond the likelihood-root
#' cutoffs for `level` in both directions on the selected root scale.
#' Consequently, a two-sided interval cannot be extracted from a single-branch
#' profile.
#'
#' Methods `"mpl"` and `"rstar"` require an object returned by
#' [modified_profile_focus()]. For `"mpl"`, the modified profile is recentered
#' at the maximum of its cubic-spline interpolant before constructing the
#' signed likelihood root based on the modified profile likelihood. The
#' `interpolation` argument then determines how that root is interpolated to
#' obtain the confidence limits. The modified profile maximum must be in the
#' interior of the supplied grid.
#'
#' In contrast, `confint(focus_object, method = "pl")` uses [profile_ci()]
#' to solve the endpoint equations directly.
#'
#' @return
#' A numeric vector of length 2 with names `"lower"` and `"upper"`. The
#' `"level"`, `"type"`, and `"interpolation"` attributes record the nominal
#' level, interval type (`"pl"`, `"mpl"`, or `"rstar"`), and interpolation
#' method, respectively.
#'
#' @seealso [profile_focus()], [profile.focus_list_glm()],
#'     [modified_profile_focus()], [plot.profile_focus_list()], [profile_ci()]
#'
#' @export
confint.profile_focus_list <- function(object,
                                       parm,
                                       level = 0.95,
                                       method = c("pl", "mpl", "rstar"),
                                       interpolation = c("linear", "cubic"),
                                       ...) {
    if (!is.numeric(level) || length(level) != 1L ||
        !is.finite(level) || level <= 0 || level >= 1)
        stop("`level` must be a number in (0, 1).")
    method <- match.arg(method)
    interpolation <- match.arg(interpolation)
    if (!identical(method, "pl") &&
        !inherits(object, "modified_profile_focus_list"))
        stop("`method = \"", method, "\"` requires an object returned by ",
             "`modified_profile_focus()`.")
    signed_root <- switch(
        method,
        pl = .profile_signed(object),
        mpl = {
            maximum <- .modified_profile_maximum(object)
            sign(maximum["psi"] - object$psi) *
                sqrt(2 * pmax(maximum["loglik"] -
                              object$modified_loglik, 0))
        },
        rstar = object$rstar
    )
    keep <- is.finite(object$psi) & is.finite(signed_root)
    psi <- object$psi[keep]
    signed_root <- signed_root[keep]
    if (length(signed_root) < 2L)
        stop("The supplied profile does not contain enough finite values for ",
             "`method = \"", method, "\"`.")
    qua <- qnorm(0.5 + level / 2)
    if (qua > max(signed_root) || -qua < min(signed_root)) {
        stop("`level` exceeds the range of the supplied profile; ",
             "try recomputing it with a larger `max_level` or a wider ",
             "`focus_range`, or inspect its behavior using ",
             "`approach = \"focus_grid\"`.")
    }
    fn <- if (identical(interpolation, "linear"))
        approxfun(x = signed_root, y = psi)
    else
        splinefun(signed_root, y = psi)
    out <- c(lower = fn(qua), upper = fn(-qua))
    attr(out, "level") <- level
    attr(out, "type") <- method
    attr(out, "interpolation") <- interpolation
    out
}
