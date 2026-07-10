library("nleqslv")

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
                      information = information_mean,
                      mle = theta_hat,
                      on = on_mean,
                      on_gradient = on_mean_gradient,
                      level = level)
expected_mean <- theta_hat + c(lower = -1, upper = 1) * sqrt(cutoff / n)

expect_equal(ci_mean, expected_mean, tolerance = 1e-08, check.attributes = FALSE)
expect_identical(attr(ci_mean, "type"), "profile")
expect_true(max(attr(ci_mean, "max|fvec|")) < 1e-08)
expect_true(is.character(attr(ci_mean, "messages")))

expect_error(profile_ci(loglik = loglik_mean,
                        score = function(theta) c(1, 2),
                        mle = theta_hat,
                        on = on_mean,
                        on_gradient = on_mean_gradient,
                        level = level))
expect_error(profile_ci(loglik = loglik_mean,
                        score = score_mean,
                        information = function(theta) matrix(1, 2, 2),
                        mle = theta_hat,
                        on = on_mean,
                        on_gradient = on_mean_gradient,
                        level = level))
expect_error(profile_ci(loglik = loglik_mean,
                        score = score_mean,
                        mle = theta_hat,
                        on = on_mean,
                        on_gradient = function(theta) c(1, 2),
                        level = level))


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
                    information = information_normal,
                    mle = theta2_hat,
                    on = on_mu,
                    on_gradient = on_mu_gradient,
                    level = level)
rss_hat <- sum((y - theta2_hat[1])^2)
expected_mu <- theta2_hat[1] +
    c(lower = -1, upper = 1) * sqrt(rss_hat * (exp(cutoff / n) - 1) / n)

expect_equal(ci_mu, expected_mu, tolerance = 1e-07, check.attributes = FALSE)
expect_true(max(attr(ci_mu, "max|fvec|")) < 1e-08)


ci_numeric <- profile_ci(loglik = loglik_normal,
                         on = on_mu,
                         mle = theta2_hat,
                         level = level)

expect_equal(ci_numeric, expected_mu, tolerance = 1e-05, check.attributes = FALSE)


budworm <- data.frame(ldose = rep(0:5, 2),
                      numdead = c(1, 4, 9, 13, 18, 20, 0, 2, 6, 10, 12, 16),
                      sex = factor(rep(c("M", "F"), c(6, 6))))
budworm <- transform(budworm, numalive = 20 - numdead)
budworm_lg <- glm(cbind(numalive, numdead) ~ sex * ldose,
                  family = binomial, data = budworm)
budworm_profile <- profile(budworm_lg)
budworm_loglik <- function(theta, object) {
    dfun <- enrichwith::get_dmodel_function(object)
    sum(dfun(coefficients = theta, log = TRUE))
}
budworm_score <- function(theta, object) {
    sfun <- enrichwith::get_score_function(object)
    sfun(coefficients = theta)
}
budworm_information <- function(theta, object) {
    ifun <- enrichwith::get_information_function(object)
    ifun(coefficients = theta)
}
second_parameter <- function(theta) theta[2]
second_parameter_gradient <- function(theta) {
    out <- numeric(length(theta))
    out[2] <- 1
    out
}

for (level_budworm in c(0.80, 0.90, 0.95)) {
    ci_budworm <- profile_ci(loglik = budworm_loglik,
                             score = budworm_score,
                             information = budworm_information,
                             mle = coef(budworm_lg),
                             on = second_parameter,
                             on_gradient = second_parameter_gradient,
                             likelihood_args = list(object = budworm_lg),
                             level = level_budworm)
    ci_profile <- confint(budworm_profile, parm = 2, level = level_budworm)
    expect_equal(ci_budworm, ci_profile,
                 tolerance = 1e-04, check.attributes = FALSE)
}
