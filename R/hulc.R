#' Compute the number of HulC partitions
#'
#' Determine the number of data partitions, `B`, required to construct a
#' HulC confidence interval with nominal miscoverage level `alpha` and
#' asymmetry parameter `Delta`.
#'
#' The value of `B` is chosen as the smallest integer satisfying the
#' HulC inequality
#'
#' \deqn{(1/2 - \Delta)^B + (1/2 + \Delta)^B \le \alpha.}
#'
#' This implementation follows the HulC construction of Kuchibhotla et al. (2024)
#'
#'
#' @param alpha Numeric scalar in `(0, 1)`. Target miscoverage level.
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
#' Compared with Kuchibholta et al. (2024, expression (7)), the upper search
#' bound here uses `ceiling()` to be slightly more conservative.
#'
#' @references
#' Kuchibhotla A K, Balakrishnan S, Wasserman L (2024). *Journal of
#' the Royal Statistical Society Series B: Statistical
#' Methodology*, **86**, 586-622. \doi{10.1093/jrsssb/qkad134}.
#'
#' @examples
#' compute_B(alpha = 0.05, Delta = 0)
#' compute_B(alpha = 0.1, Delta = 0.05)
#'
#' @export
compute_B <- function(alpha, Delta){
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
#' The number of partitions is chosen to achieve nominal miscoverage level
#' `alpha` under the HulC construction, optionally with randomization as
#' described in the original method.
#'
#' @param data A [`data.frame`] object with observations in rows.
#' @param statistic A function that takes a [`data.frame`] as its first argument
#'   and returns a single numeric value.
#' @param alpha Numeric scalar in `(0, 1)`. Target miscoverage level.
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
#'   \item{`alpha`}{The supplied target miscoverage level.}
#'   \item{`B`}{The number of partitions used.}
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
#' the randomized HulC construction (see Kuchibholta et al, 2025,
#' Section 2.1).
#'
#' The function assumes that `statistic` returns a single numeric value for each
#' subset. If `statistic` returns `NA`, throws an error, or is undefined for
#' small subsets, the function returns `c(NA, NA)` for the interval.
#'
#' @references
#'
#' Kuchibhotla A K, Balakrishnan S, Wasserman L (2024). *Journal of
#' the Royal Statistical Society Series B: Statistical
#' Methodology*, **86**, 586-622. \doi{10.1093/jrsssb/qkad134}.
#'
#' @examples
#' set.seed(1)
#' x <- data.frame(y = rnorm(100))
#'
#' hulc_ci(
#'   data = x,
#'   statistic = function(d) mean(d$y),
#'   alpha = 0.05
#' )
#'
#' hulc_ci(
#'   data = x,
#'   statistic = function(d, trim = 0.1) mean(d$y, trim = trim),
#'   alpha = 0.1,
#'   trim = 0.2
#' )
#'
#' @export
hulc_ci <- function(data,
                    statistic,
                    alpha = 0.05,
                    Delta = 0,
                    randomize = TRUE,
                    check_statistic = TRUE,
                    ...) {
    nobs <- nrow(data)
    data <- data[sample(nobs), , drop = FALSE]
    B <- compute_B(alpha, Delta)
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
    test <- error_msg <- NA
    if (check_statistic) {
        min_id <- which.min(sapply(data, nrow))
        small_data <- data[[min_id]]
        stat <- try(statistic(small_data, ...), silent = TRUE)
        test <- isTRUE(inherits(stat, "try-error") || is.na(stat))
        error_msg <- stat[1]
    }
    if (test) {
        warning("It has not been possible to evaluate the statistic on the partition with the smallest number of observations (=", nrow(small_data), ").")
        ci <- c(NA, NA)
    } else {
        ci <- try(vapply(data, statistic, numeric(1), ...) |> range(), silent = TRUE)
        if (inherits(ci, "try-error")) {
            ci <- c(NA, NA)
            error_msg <- stat[1]
        }
        ## NA's in ci's are handled by `range()` where `na.rm = FALSE`
    }
    names(ci) <- c("lower", "upper")
    attr(ci, "Delta") <- Delta
    attr(ci, "alpha") <- alpha
    attr(ci, "B") <- B
    attr(ci, "error") <- error_msg
    attr(ci, "type") <- "hulc"
    ci
}
