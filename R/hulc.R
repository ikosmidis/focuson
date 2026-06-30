#' Compute the number of HulC partitions
#'
#' Determine the number of data partitions, `B`, required to construct
#' a HulC confidence interval with nominal coverage `level` and median
#' bias parameter `Delta`.
#'
#' The value of `B` is chosen as the smallest integer satisfying the
#' HulC inequality
#'
#' \deqn{(1/2 - \Delta)^B + (1/2 + \Delta)^B \le \alpha.}
#'
#' where \eqn{\alpha} is `1 - level`.
#'
#' This implementation follows the HulC construction of Kuchibhotla et al. (2024)
#'
#'
#' @param level Numeric scalar in `(0, 1)`. Target coverage level.
#' @param Delta Numeric scalar in `[0, 1/2)`. The median bias of the
#'     statistic; see [focus()]. `Delta = 0` corresponds to zero
#'     median bias.
#'
#' @return An integer scalar giving the required number of partitions.
#'
#' @details
#' The search is restricted to a finite interval of candidate values for `B`,
#' then the first value satisfying the HulC inequality is returned.
#'
#' Compared with Kuchibhotla et al. (2024, expression (7)), the upper search
#' bound here uses `ceiling()` to be slightly more conservative.
#'
#' @references
#'
#' Kuchibhotla A K, Balakrishnan S, Wasserman L (2024). The HulC:
#' confidence regions from convex hulls. *Journal of the Royal
#' Statistical Society Series B: Statistical Methodology*, **86**,
#' 586-622. \doi{10.1093/jrsssb/qkad134}.
#'
#' @examples
#' compute_B(level = 0.95, Delta = 0)
#' compute_B(level = 0.90, Delta = 0.05)
#'
#' @export
compute_B <- function(level, Delta){
    stopifnot(length(level) == 1L, is.numeric(level), !is.na(level), level > 0, level < 1)
    stopifnot(length(Delta) == 1L, is.numeric(Delta), !is.na(Delta), Delta >= 0, Delta < 1/2)
    alpha <- 1 - level
    log2 <- log(2)
    loga <- log(alpha)
    Dp <- 1/2 + Delta
    Dm <- 1/2 - Delta
    logDp <- log(Dp)
    log2mloga <- log2 - loga
    ## https://doi.org/10.1093/jrsssb/qkad134 has floor but let's put
    ## ceiling here to be slightly more conservative
    Bl <- max(floor(loga / logDp), floor(log2mloga / log2))
    Bu <- ceiling(- log2mloga / logDp)
    for (B in Bl:Bu)
        if (Dm^B + Dp^B <= alpha) break
    as.integer(B)
}

#' Construct a HulC confidence interval
#'
#' Compute a HulC confidence interval for a user-supplied statistic by
#' partitioning the data into `B` subsets, evaluating the statistic on each
#' subset, and taking the range of the resulting values.
#'
#' The number of partitions is chosen to achieve nominal coverage
#' level `level` under the HulC construction, optionally with
#' randomization as described in the original method.
#'
#' @param data A [`data.frame`] object with observations in rows.
#' @param statistic A function that takes a [`data.frame`] as its first argument
#'   and returns a single numeric value.
#' @param level Numeric scalar in `(0, 1)`. Target coverage level.
#' @param Delta Numeric scalar in `[0, 1/2)`. The median bias of the `statistic`.
#' @param randomize Logical. If `TRUE` (default), randomize between `B` and `B - 1` to
#'   obtain less conservative finite-sample coverage.
#' @param check_statistic Logical. If `TRUE`, first evaluate `statistic` on the
#'   smallest partition to check whether it can be computed on a minimal subset.
#' @param ... Additional arguments passed to `statistic`.
#'
#' @return
#' A numeric vector of length 2 with elements `"lower"` and `"upper"`.
#'
#' The return value has the following attributes:
#' \describe{
#'   \item{`Delta`}{The supplied median bias.}
#'   \item{`level`}{The supplied target coverage level.}
#'   \item{`B`}{The number of partitions used.}
#'   \item{`error`}{`NA` on success, or the captured error message if
#'     evaluating `statistic` failed.}
#'   \item{`type`}{Character string `"hulc"`.}
#' }
#'
#' If the statistic cannot be computed on the smallest partition, or fails on
#' one of the partitions, the function returns `c(NA, NA)`.
#'
#' @details
#' The observations are first randomly permuted, then split into `B` roughly
#' equal-sized partitions. The user-supplied `statistic` is evaluated on each
#' partition, and the confidence interval is defined as the range of these
#' values.
#'
#' When `randomize = TRUE`, the function may reduce the number of
#' partitions from `B` to `B - 1` with a probability chosen to match
#' the randomized HulC construction (see Kuchibhotla et al, 2024,
#' Section 2.1).
#'
#' If `check_statistic = TRUE` and the statistic cannot be evaluated on the
#' smallest partition, the function issues a warning and returns `c(NA, NA)`.
#' If evaluation instead fails later when applying `statistic` across all
#' partitions, the function also issues a warning, returns `c(NA, NA)`, and
#' stores the captured error message in the `"error"` attribute.
#'
#' The function assumes that `statistic` returns a single numeric value for each
#' subset. If `statistic` returns `NA`, throws an error, or is undefined for
#' small subsets, the function returns `c(NA, NA)` for the interval.
#'
#' @references
#'
#' Kuchibhotla A K, Balakrishnan S, Wasserman L (2024). The HulC:
#' confidence regions from convex hulls. *Journal of the Royal
#' Statistical Society Series B: Statistical Methodology*, **86**,
#' 586-622. \doi{10.1093/jrsssb/qkad134}.
#'
#' @examples
#' set.seed(1)
#' x <- data.frame(y = rnorm(100))
#'
#' hulc_ci(
#'   data = x,
#'   statistic = function(d) mean(d$y),
#'   level = 0.95
#' )
#'
#' hulc_ci(
#'   data = x,
#'   statistic = function(d, trim = 0.1) mean(d$y, trim = trim),
#'   level = 0.90,
#'   trim = 0.2
#' )
#'
#' @export
hulc_ci <- function(data,
                    statistic,
                    level = 0.95,
                    Delta = 0,
                    randomize = TRUE,
                    check_statistic = TRUE,
                    ...) {
    nobs <- nrow(data)
    stopifnot(is.data.frame(data), is.function(statistic),
        length(level) == 1L, is.numeric(level), !is.na(level), level > 0, level < 1,
        length(Delta) == 1L, is.numeric(Delta), !is.na(Delta), Delta >= 0, Delta < 1/2,
        length(randomize) == 1L, is.logical(randomize), !is.na(randomize),
        length(check_statistic) == 1L, is.logical(check_statistic), !is.na(check_statistic))
    data <- data[sample(nobs), , drop = FALSE]
    B <- compute_B(level, Delta)
    alpha <- 1 - level
    if (randomize) {
        Dp <- 0.5 + Delta
        Dm <- 0.5 - Delta
        pB <- c(Dm^B, Dp^B)
        prob0 <- sum(pB)
        prob1 <- sum(pB / c(Dm, Dp))
        B <- B - ((alpha - prob0) / (prob1 - prob0) >= runif(1))
    }
    if (B > nobs)
        stop("The required number of partitions (=", B, ") is larger than `nrow(data)` (=", nobs, ").")
    partitions <- ((seq_len(nobs) - 1) * B) %/% nobs + 1
    data <- split(data, partitions)
    ## Check if statistic returns a value or fails for the smallest
    ## number of observations
    test <- FALSE
    error_msg <- NA_character_
    if (check_statistic) {
        min_id <- which.min(vapply(data, nrow, integer(1)))
        small_data <- data[[min_id]]
        stat <- try(statistic(small_data, ...), silent = TRUE)
        got_error <- inherits(stat, "try-error")
        error_msg <- if (got_error) as.character(stat[1]) else NA_character_
        test <- isTRUE(got_error || length(stat) != 1L ||
                       !is.numeric(stat) || is.na(stat) ||
                       !is.finite(stat))
        if (test && !got_error) {
            error_msg <- paste(
                "The statistic must return a single finite numeric value;",
                "got", paste(capture.output(dput(stat)), collapse = "")
            )
        }
    }
    if (test) {
        warning(
            "It has not been possible to evaluate the statistic on the partition ",
            "with the smallest number of observations (=", nrow(small_data), "). ",
            error_msg
        )
        ci <- c(NA_real_, NA_real_)
    } else {
        stats <- try(vapply(data, statistic, numeric(1), ...), silent = TRUE)
        if (inherits(stats, "try-error")) {
            error_msg <- as.character(stats[1])
            ci <- c(NA_real_, NA_real_)
            warning("It has not been possible to evaluate the statistic on all partitions. ", error_msg)
        } else if (anyNA(stats) || any(!is.finite(stats))) {
            error_msg <- "The statistic returned a non-finite value on at least one partition."
            ci <- c(NA_real_, NA_real_)
            warning("It has not been possible to evaluate the statistic on all partitions. ", error_msg)
        } else {
            ci <- range(stats)
        }
    }

    names(ci) <- c("lower", "upper")
    attr(ci, "Delta") <- Delta
    attr(ci, "level") <- level
    attr(ci, "B") <- B
    attr(ci, "error") <- error_msg
    attr(ci, "type") <- "hulc"
    ci
}
