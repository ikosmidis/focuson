#' Confidence intervals for focus objects
#'
#' Compute confidence intervals for objects returned by [focus()].
#'
#' @param object An object of class `"focus_list"`.
#' @param parm Currently unused.
#' @param level Confidence level.
#' @param method Character string specifying the confidence interval
#'     method.  One of `"wald"` or `"hulc"`.
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
#' For `method = "hulc"`, [hulc_ci()] is applied to the model frame of the
#' stored fitted object using [focus_statistic()] as the statistic evaluated
#' on each partition. This requires that [stats::model.frame()] and
#' [stats::update()] work for the stored fitted object.
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
                               se_at = c("supplied", "compatible"),
                               V_function = NULL,
                               se_control = list(),
                               ...) {
    method <- match.arg(method, c("wald", "hulc"))
    se_at <- match.arg(se_at)
    if (!is.list(se_control)) {
        stop("`se_control` must be a list.")
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

    if (is.null(object$object)) {
        stop("`method = \"hulc\"` requires a `focus()` result with a stored fitted object.")
    }

    correction <- object$correction
    on_funs <- object$on
    on <- on_funs$on
    on_gradient <- on_funs$on_gradient
    on_hessian <- on_funs$on_hessian
    odots <- object$dots
    statistic <- function(data) {
        do.call(
            focus_statistic,
            c(list(data = data,
                   object = object$object,
                   on = on,
                   correction = correction,
                   on_gradient = on_gradient,
                   on_hessian = on_hessian),
              odots
              )
        )
    }
    do.call(hulc_ci,
            c(list(data = model.frame(object$object),
                   statistic = statistic,
                   level = level),
              list(...)))
}
