B <- compute_B(alpha = 0.05, Delta = 0)
expect_true(is.integer(B))
expect_length(B, 1)
expect_true(B >= 1)

alpha <- 0.05
Delta <- 0.1
B <- compute_B(alpha, Delta)
Dm <- 0.5 - Delta
Dp <- 0.5 + Delta
expect_true(Dm^B + Dp^B <= alpha)
expect_true(Dm^(B - 1) + Dp^(B - 1) >= alpha)


set.seed(1)
x <- data.frame(y = rnorm(100))
ci <- hulc_ci(
    data = x,
    statistic = function(d) mean(d$y),
    alpha = 0.2,
    Delta = 0.123,
    randomize = FALSE
)
expect_true(is.double(ci))
expect_equal(length(ci), 2)
expect_true(identical(names(ci), c("lower", "upper")))
expect_true(ci["lower"] <= ci["upper"])
expect_identical(attr(ci, "type"), "hulc")


expect_equal(attr(ci, "alpha"), 0.2)
expect_equal(attr(ci, "Delta"), 0.123)
expect_equal(attr(ci, "B"), compute_B(0.2, 0.123))
expect_identical(attr(ci, "type"), "hulc")


set.seed(1)
x <- data.frame(y = rnorm(10))
stat <- function(d) {
    if (nrow(d) < 3) stop("too small")
    mean(d$y)
}

expect_warning(
    ci <- hulc_ci(data = x, statistic = stat, alpha = 0.05,
                  randomize = FALSE, check_statistic = TRUE)
)

expect_true(all(is.na(ci)))
