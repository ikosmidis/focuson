set.seed(1)
y <- rnorm(20, mean = 0.5, sd = 1)
n <- length(y)
theta_hat <- mean(y)
level <- 0.95
cutoff <- qchisq(level, df = 1)

loglik_mean <- function(theta) {
    -0.5 * sum((y - theta[1])^2)
}
score_mean <- function(theta) {
    sum(y - theta[1])
}
information_mean <- function(theta) {
    matrix(n, 1, 1)
}
on_mean <- function(theta) theta[1]
on_mean_gradient <- function(theta) 1
on_mean_hessian <- function(theta) matrix(0, 1, 1)

ci_mean <- profile_ci(loglik = loglik_mean,
                      score = score_mean,
                      mle = theta_hat,
                      on = on_mean,
                      on_gradient = on_mean_gradient,
                      level = level)
expected_mean <- theta_hat + c(lower = -1, upper = 1) * sqrt(cutoff / n)

expect_equal(ci_mean, expected_mean, tolerance = 1e-08, check.attributes = FALSE)
expect_identical(attr(ci_mean, "type"), "profile")
expect_equal(attr(ci_mean, "level"), level)
expect_equal(attr(ci_mean, "cutoff"), cutoff)
expect_true(is.list(attr(ci_mean, "details")))


theta2_hat <- c(mu = mean(y),
                log_sigma = log(sqrt(mean((y - mean(y))^2))))
loglik_normal <- function(theta) {
    mu <- theta[1]
    sigma <- exp(theta[2])
    -n * log(sigma) - 0.5 * sum((y - mu)^2) / sigma^2
}
score_normal <- function(theta) {
    mu <- theta[1]
    sigma <- exp(theta[2])
    rss <- sum((y - mu)^2)
    c(sum(y - mu) / sigma^2,
      -n + rss / sigma^2)
}
information_normal <- function(theta) {
    mu <- theta[1]
    sigma <- exp(theta[2])
    c11 <- n / sigma^2
    c12 <- 2 * sum(y - mu) / sigma^2
    c22 <- 2 * sum((y - mu)^2) / sigma^2
    matrix(c(c11, c12, c12, c22), 2, 2)
}
on_mu <- function(theta) theta[1]
on_mu_gradient <- function(theta) c(1, 0)
on_mu_hessian <- function(theta) matrix(0, 2, 2)

ci_mu <- profile_ci(loglik = loglik_normal,
                    score = score_normal,
                    mle = theta2_hat,
                    on = on_mu,
                    on_gradient = on_mu_gradient,
                    level = level)
rss_hat <- sum((y - theta2_hat[1])^2)
expected_mu <- theta2_hat[1] +
    c(lower = -1, upper = 1) * sqrt(rss_hat * (exp(cutoff / n) - 1) / n)

expect_equal(ci_mu, expected_mu, tolerance = 1e-07, check.attributes = FALSE)
expect_true(max(abs(attr(ci_mu, "details")$lower$residual)) < 1e-08)
expect_true(max(abs(attr(ci_mu, "details")$upper$residual)) < 1e-08)


ci_numeric <- profile_ci(loglik = loglik_normal,
                         on = on_mu,
                         mle = theta2_hat,
                         level = level)

expect_equal(ci_numeric, expected_mu, tolerance = 1e-05, check.attributes = FALSE)
