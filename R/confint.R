#' Confidence intervals for focus objects
#'
#' Compute confidence intervals for objects returned by [focus()].
#'
#' @param object An object of class `"focus_list"`.
#' @param parm Currently unused.
#' @param level Confidence level.
#' @param method Character string specifying the confidence interval method.
#'   One of `"wald"` or `"hulc"`.
#' @param ... Additional arguments for the confidence interval method. For
#'   `method = "hulc"`, these are passed to [hulc_ci()], except that the
#'   miscoverage rate is determined by `level`.
#'
#' @return
#' A numeric vector of length 2 with names `"lower"` and `"upper"`.
#'
#' @details
#' For `method = "wald"`, the interval is computed from the stored point
#' estimate and standard error using the normal approximation.
#'
#' For `method = "hulc"`, [hulc_ci()] is applied to the model frame of the
#' stored fitted object using [focus_statistic()] as the statistic evaluated
#' on each partition. This requires that [stats::model.frame()] and
#' [stats::update()] work for the stored fitted object.
#'
#' The coverage rate is determined by `level`; users should not supply
#' `level` in `...`.
#'
#' @seealso [focus()], [hulc_ci()], [focus_statistic()]
#'
#' @export
confint.focus_list <- function(object,
                               parm,
                               level = 0.95,
                               method = "wald",
                               ...) {
    method <- match.arg(method, c("wald", "hulc"))
    alpha <- 1 - level

    if (identical(method, "wald")) {
        ci <- unname(object$estimate) + c(-1, 1) * qnorm(1 - alpha / 2) * object$se
        names(ci) <- c("lower", "upper")
        attr(ci, "level") <- level
        attr(ci, "type") <- "wald"
        return(ci)
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
