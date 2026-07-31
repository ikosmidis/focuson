set.seed(1)
y <- rnorm(20, mean = 0.5)
n <- length(y)
mle <- c(mu = mean(y))

loglik <- function(theta, data, n)
    sum(dnorm(data, mean = theta[1], sd = 1, log = TRUE))
score <- function(theta, data, n)
    c(sum(data - theta[1]))
information <- function(theta, data, n)
    matrix(n, 1, 1)
simulate <- function(theta, n)
    rnorm(n, mean = theta[1])

set.seed(2)
modified <- modified_profile_focus(
    loglik = loglik,
    score = score,
    information = information,
    simulate = simulate,
    mle = mle,
    on_gradient = function(theta) 1,
    on_hessian = function(theta) matrix(0, 1, 1),
    likelihood_args = list(data = y, n = n),
    nsim = 100,
    grid_size = 3,
    max_level = 0.9
)

away_from_mle <- modified$signed_root != 0
expect_true(inherits(modified, "modified_profile_focus_list"))
expect_true(inherits(modified, "profile_focus_list"))
expect_equal(modified$log_C_inverse, rep(0, nrow(modified)))
expect_equal(modified$NP[away_from_mle],
             rep(0, sum(away_from_mle)))
expect_equal(modified$modified_loglik, modified$loglik)
expect_equal(modified$rstar[away_from_mle],
             modified$signed_root[away_from_mle],
             tolerance = 1e-10)
expect_true(is.na(modified$rstar[!away_from_mle]))
expect_identical(attr(modified, "nsim"), 100)
expect_identical(dim(attr(modified, "expected_information")), c(1L, 1L))

printed_modified <- capture.output(print(modified))
expect_true(any(grepl("Profile log-likelihood for a scalar focus",
                      printed_modified, fixed = TRUE)))
expect_true(any(grepl("Modified-profile calculations:",
                      printed_modified, fixed = TRUE)))
expect_true(any(grepl("Simulations: 100",
                      printed_modified, fixed = TRUE)))
expect_true(any(grepl("Usable noncentral points: 6 of 6",
                      printed_modified, fixed = TRUE)))
expect_true(any(grepl("log|u_tilde / r| range:",
                      printed_modified, fixed = TRUE)))

modified_nonfinite <- modified
modified_nonfinite$log_u_over_r[which(away_from_mle)[1L]] <- NA_real_
invisible(capture.output(
    expect_warning(
        print(modified_nonfinite),
        pattern = "non-finite at 1 of 6 noncentral profile points"
    )
))

level <- 0.8
expected_ci <- mle +
    c(lower = -1, upper = 1) * qnorm(0.5 + level / 2) / sqrt(n)
ci_pl <- confint(modified, level = level, method = "pl")
ci_mpl <- confint(modified, level = level, method = "mpl")
ci_mpl_cubic <- confint(modified, level = level, method = "mpl",
                        interpolation = "cubic")
ci_rstar <- confint(modified, level = level, method = "rstar")
expect_equal(ci_pl, expected_ci,
             tolerance = 1e-08, check.attributes = FALSE)
expect_equal(ci_mpl, expected_ci,
             tolerance = 1e-08, check.attributes = FALSE)
expect_equal(ci_mpl_cubic, expected_ci,
             tolerance = 1e-08, check.attributes = FALSE)
expect_equal(ci_rstar, expected_ci,
             tolerance = 1e-08, check.attributes = FALSE)
expect_identical(attr(ci_pl, "type"), "pl")
expect_identical(attr(ci_mpl, "type"), "mpl")
expect_identical(attr(ci_rstar, "type"), "rstar")
modified_plot_file <- tempfile(fileext = ".pdf")
pdf(modified_plot_file)
expect_silent(plot(modified, what = "pl"))
expect_silent(plot(modified, what = "mpl"))
expect_silent(plot(modified, what = "mpl", signed = TRUE))
expect_silent(plot(modified, what = "mpl", level = level, ci = TRUE))
expect_silent(plot(modified, what = "rstar"))
expect_silent(plot(modified, what = "rstar", signed = TRUE))
expect_silent(plot(modified, what = "rstar", signed = TRUE,
                   level = level, ci = TRUE))
dev.off()
unlink(modified_plot_file)

boundary_modified <- modified
boundary_modified$modified_loglik <- seq_len(nrow(boundary_modified))
expect_error(confint(boundary_modified, level = level, method = "mpl"),
             pattern = "grid boundary")


set.seed(3)
y <- rnorm(25, mean = 0.5, sd = 1.2)
n <- length(y)
mle <- c(mu = mean(y),
         log_sigma = 0.5 * log(mean((y - mean(y))^2)))

loglik <- function(theta, data, n)
    sum(dnorm(data, mean = theta[1], sd = exp(theta[2]), log = TRUE))
score <- function(theta, data, n) {
    residuals <- data - theta[1]
    inverse_variance <- exp(-2 * theta[2])
    c(sum(residuals) * inverse_variance,
      -n + sum(residuals^2) * inverse_variance)
}
information <- function(theta, data, n) {
    residuals <- data - theta[1]
    inverse_variance <- exp(-2 * theta[2])
    matrix(c(n * inverse_variance,
             2 * sum(residuals) * inverse_variance,
             2 * sum(residuals) * inverse_variance,
             2 * sum(residuals^2) * inverse_variance),
           2, 2)
}
simulate <- function(theta, n)
    rnorm(n, mean = theta[1], sd = exp(theta[2]))

set.seed(4)
modified_with_nuisance <- modified_profile_focus(
    loglik = loglik,
    score = score,
    information = information,
    simulate = simulate,
    mle = mle,
    on_gradient = function(theta) c(1, 0),
    on_hessian = function(theta) matrix(0, 2, 2),
    likelihood_args = list(data = y, n = n),
    nsim = 100,
    grid_size = 2,
    max_level = 0.8
)

away_from_mle <- modified_with_nuisance$signed_root != 0
expect_true(all(is.finite(
    modified_with_nuisance$log_C_inverse
)))
expect_true(any(abs(
    modified_with_nuisance$log_C_inverse[away_from_mle]
) > 1e-06))
expect_true(all(is.finite(
    modified_with_nuisance$rstar[away_from_mle]
)))
expect_equal(
    modified_with_nuisance$modified_loglik,
    modified_with_nuisance$loglik -
        modified_with_nuisance$log_C_inverse
)
