.focus_ci <- function(x) {
    class(x) <- unique(c("focus_ci", class(x)))
    x
}


#' Print a focus confidence interval
#'
#' Print a confidence interval returned by pkg{focuson}, together with a
#' concise summary of the method-specific information used to construct it.
#'
#' @param x An object of class `"focus_ci"`.
#' @param digits Number of significant digits used for printing.
#' @param ... Additional arguments passed to [print()].
#'
#' @details
#' Method-specific diagnostic attributes remain available through
#' [attributes()] but are not all printed. The method reports the location of
#' the standard-error evaluation for Wald intervals, endpoint diagnostics for
#' intervals computed directly by the Venzon--Moolgavkar equations,
#' interpolation for intervals extracted from a computed profile, and the
#' number of partitions and median-bias bound for HulC intervals.
#'
#' @return `x`, invisibly.
#'
#' @export
print.focus_ci <- function(x,
                           digits = max(3L, getOption("digits") - 2L),
                           ...) {
    type <- attr(x, "type")
    level <- attr(x, "level")
    description <- switch(
        type,
        wald = "Wald confidence interval",
        pl = "profile likelihood confidence interval",
        mpl = "modified profile likelihood confidence interval",
        rstar = "confidence interval based on r*",
        hulc = "HulC confidence interval",
        "focus confidence interval"
    )
    level_label <- if (is.numeric(level) && length(level) == 1L &&
                       is.finite(level)) {
        paste0(format(100 * level, digits = digits, trim = TRUE), "% ")
    } else {
        ""
    }
    cat(level_label, description, "\n\n", sep = "")
    interval <- setNames(as.numeric(x), names(x))
    print(interval, digits = digits, ...)

    if (identical(type, "wald")) {
        se_at <- attr(x, "se_at")
        if (!is.null(se_at))
            cat("\nStandard error evaluated at: ", se_at, "\n", sep = "")
    }

    if (type %in% c("pl", "mpl", "rstar")) {
        interpolation <- attr(x, "interpolation")
        if (!is.null(interpolation)) {
            cat("\nInterpolation: ", interpolation, "\n", sep = "")
        } else if (identical(type, "pl")) {
            iter <- attr(x, "iter")
            if (!is.null(iter))
                cat("\nEndpoint iterations: ",
                    paste(paste0(names(iter), " = ", iter), collapse = ", "),
                    "\n", sep = "")
            residuals <- attr(x, "max|fvec|")
            if (!is.null(residuals) && any(is.finite(residuals)))
                cat("Maximum equation residual: ",
                    format(max(residuals[is.finite(residuals)]),
                           digits = digits, trim = TRUE),
                    "\n", sep = "")
        }
    }

    if (identical(type, "hulc")) {
        cat("\nPartitions: ", attr(x, "B"), "\n",
            "Median-bias bound (Delta): ", attr(x, "Delta"), "\n",
            sep = "")
        error <- attr(x, "error")
        if (!is.null(error) && any(!is.na(error)))
            cat("Stored error: ", paste(error[!is.na(error)], collapse = "; "),
                "\n", sep = "")
    }
    invisible(x)
}
