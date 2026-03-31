set.seed(1)
x <- data.frame(y = rnorm(8))

expect_warning(
    ci <- hulc_ci(
        data = x,
        statistic = function(d) if (nrow(d) == 2) Inf else mean(d$y),
        level = 0.8,
        randomize = FALSE,
        check_statistic = FALSE
    )
)

expect_true(all(is.na(ci)))
expect_true(is.numeric(ci))
expect_identical(
    attr(ci, "error"),
    "The statistic returned a non-finite value on at least one partition."
)


set.seed(1)
x <- data.frame(y = rnorm(20))
ci <- hulc_ci(
    data = x,
    statistic = function(d) mean(d$y),
    level = 0.8,
    randomize = FALSE
)

expect_true(!anyNA(ci))
expect_identical(attr(ci, "error"), NA_character_)


set.seed(1)
x <- data.frame(y = rnorm(10))

expect_warning(
    ci <- hulc_ci(
        data = x,
        statistic = function(d) c(mean(d$y), sd(d$y)),
        level = 0.95,
        randomize = FALSE,
        check_statistic = TRUE
    )
)

expect_true(all(is.na(ci)))
expect_true(is.numeric(ci))
expect_true(grepl(
    "The statistic must return a single finite numeric value;",
    attr(ci, "error"),
    fixed = TRUE
))


set.seed(1)
x <- data.frame(y = rnorm(10))

expect_warning(
    ci <- hulc_ci(
        data = x,
        statistic = function(d) if (nrow(d) < 3) NA_real_ else mean(d$y),
        level = 0.95,
        randomize = FALSE,
        check_statistic = TRUE
    )
)

expect_true(all(is.na(ci)))
expect_true(is.numeric(ci))
expect_true(grepl(
    "The statistic must return a single finite numeric value;",
    attr(ci, "error"),
    fixed = TRUE
))
