#' Confidence intervals for focus objects
#'
#' Compute confidence intervals for objects returned by [focus()].
#'
#' @param object An object of class `"focus_list"`.
#' @param parm Currently unused.
#' @param level Confidence level.
#' @param method Character string specifying the confidence interval
#'     method.  One of `"wald"`, `"profile"`, or `"hulc"`.
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
#'     `method = "profile"`.
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
#' For `method = "profile"`, [profile_ci()] is used to compute a likelihood
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
    method <- match.arg(method, c("wald", "profile", "hulc"))
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

    if (identical(method, "profile")) {
        if (inherits(object, "focus_engine_list")) {
            stop("`method = \"profile\"` is not available for `focus_engine()` results.")
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
    fit <- object$object
    if (!identical(fit$type, "ML")) {
        fit <- update(fit, type = "ML", start = coef(fit, model = "mean"))
    }
    p_mean <- length(coef(fit, model = "mean"))
    theta <- coef(fit, model = "full")
    if (fit$family$family %in% c("poisson", "binomial")) {
        theta <- theta[names(coef(fit, model = "mean"))]
    }
    aux <- enrichwith::get_auxiliary_functions(fit)
    split_theta <- function(theta) {
        if (length(theta) == p_mean) {
            list(coefficients = theta)
        } else {
            list(coefficients = theta[seq_len(p_mean)],
                 dispersion = theta[-seq_len(p_mean)])
        }
    }
    loglik <- function(theta) {
        args <- split_theta(theta)
        sum(do.call(aux$dmodel, c(args, list(log = TRUE))))
    }
    score <- function(theta) {
        do.call(aux$score, split_theta(theta))
    }
    information <- function(theta) {
        do.call(aux$information, split_theta(theta))
    }
    do.call(profile_ci,
            c(list(loglik = loglik,
                   score = score,
                   information = information,
                   mle = theta,
                   on = object$on$on,
                   on_gradient = object$on$on_gradient,
                   on_hessian = object$on$on_hessian,
                   level = level,
                   nleqslv_args = nleqslv_args),
              object$dots))
}
