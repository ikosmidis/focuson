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
                      on_hessian = on_mean_hessian,
                      level = level)
expected_mean <- theta_hat + c(lower = -1, upper = 1) * sqrt(cutoff / n)

expect_equal(ci_mean, expected_mean, tolerance = 1e-08, check.attributes = FALSE)
expect_identical(attr(ci_mean, "type"), "profile")
expect_true(max(attr(ci_mean, "max|fvec|")) < 1e-08)
expect_true(is.character(attr(ci_mean, "messages")))

profile_mean <- profile_focus(mle = theta_hat,
                              loglik = loglik_mean,
                              score = score_mean,
                              information = information_mean,
                              on = on_mean,
                              on_gradient = on_mean_gradient,
                              on_hessian = on_mean_hessian,
                              grid_size = 5,
                              max_level = level)
expect_true(inherits(profile_mean, "profile_focus_list"))
expect_equal(profile_mean$psi,
             theta_hat - profile_mean$signed / sqrt(n),
             tolerance = 1e-08)
expect_identical(names(profile_mean),
                 c("psi", "loglik", "signed", "theta", "lagrange"))
expect_identical(dim(profile_mean$theta),
                 c(nrow(profile_mean), length(theta_hat)))
expect_identical(colnames(profile_mean$theta), names(theta_hat))
expect_equal(as.numeric(profile_mean$theta[, 1]),
             profile_mean$psi,
             tolerance = 1e-08)
expect_equal(as.numeric(profile_mean$theta[6, ]),
             theta_hat,
             tolerance = 1e-08)
expect_equal(profile_mean$lagrange[6], 0)
expect_equal(apply(profile_mean$theta, 1, loglik_mean),
             profile_mean$loglik,
             tolerance = 1e-08)
stationarity_mean <- vapply(
    seq_len(nrow(profile_mean)),
    function(j) {
        score_mean(profile_mean$theta[j, ]) -
            profile_mean$lagrange[j] *
            on_mean_gradient(profile_mean$theta[j, ])
    },
    numeric(1)
)
expect_true(max(abs(stationarity_mean)) < 1e-08)
profile_mean_subset <- profile_mean[c(2, 5, 9), ]
expect_equal(profile_mean_subset$theta[, 1],
             profile_mean$theta[c(2, 5, 9), 1])
expect_equal(profile_mean_subset$lagrange,
             profile_mean$lagrange[c(2, 5, 9)])
printed_profile_mean <- capture.output(
    returned_profile_mean <- print(profile_mean)
)
expect_identical(returned_profile_mean, profile_mean)
expect_true(any(grepl(
    "^Profile log-likelihood for a scalar focus$",
    printed_profile_mean
)))
expect_true(any(grepl("^Parameter dimension: 1\\s*$",
                      printed_profile_mean)))
expect_true(any(grepl("^Profile points: 11\\s*$",
                      printed_profile_mean)))
expect_true(any(grepl("^Points by side: 5 lower, 5 upper\\s*$",
                      printed_profile_mean)))
expect_true(any(grepl("^Requested maximum level: 0.95\\s*$",
                      printed_profile_mean)))
expect_true(any(grepl("^Boundary confidence level:",
                      printed_profile_mean)))

profile_plot_file <- tempfile(fileext = ".pdf")
pdf(profile_plot_file)
expect_silent(plot(profile_mean, interpolation = "linear"))
expect_silent(plot(profile_mean, interpolation = "cubic"))
dev.off()
unlink(profile_plot_file)
expect_error(plot(profile_mean, interpolation = "quadratic"))

loglik_mean_with_args <- function(theta, observations) {
    -0.5 * sum((observations - theta[1])^2)
}
profile_mean_defaults <- profile_focus(
    mle = theta_hat,
    loglik = loglik_mean_with_args,
    likelihood_args = list(observations = y),
    grid_size = 5,
    max_level = level
)
expect_equal(profile_mean_defaults$psi,
             theta_hat - profile_mean_defaults$signed / sqrt(n),
             tolerance = 1e-05)

scaled_mean <- function(theta, scale) scale * theta[1]
profile_mean_scaled <- profile_focus(
    mle = theta_hat,
    loglik = loglik_mean,
    score = score_mean,
    information = information_mean,
    on = scaled_mean,
    grid_size = 5,
    max_level = level,
    scale = 2
)
expect_equal(profile_mean_scaled$psi, 2 * profile_mean$psi,
             tolerance = 1e-08)

warm_level <- 0.99
ci_mean_warm <- profile_ci(loglik = loglik_mean,
                           score = score_mean,
                           information = information_mean,
                           mle = theta_hat,
                           on = on_mean,
                           on_gradient = on_mean_gradient,
                           on_hessian = on_mean_hessian,
                           level = warm_level,
                           start = attr(ci_mean, "solution"),
                           do_checks = FALSE)
expected_mean_warm <- theta_hat + c(lower = -1, upper = 1) *
    sqrt(qchisq(warm_level, 1) / n)
expect_equal(ci_mean_warm, expected_mean_warm, tolerance = 1e-08,
             check.attributes = FALSE)

reversed_start <- rev(attr(ci_mean, "solution"))
ci_mean_reordered <- profile_ci(loglik = loglik_mean,
                                score = score_mean,
                                information = information_mean,
                                mle = theta_hat,
                                on = on_mean,
                                on_gradient = on_mean_gradient,
                                on_hessian = on_mean_hessian,
                                level = level,
                                start = reversed_start)
expect_equal(ci_mean_reordered, ci_mean, tolerance = 1e-08,
             check.attributes = FALSE)

hessian_calls <- 0L
on_mean_hessian_counted <- function(theta) {
    hessian_calls <<- hessian_calls + 1L
    matrix(0, 1, 1)
}
ci_mean_jacobian <- profile_ci(loglik = loglik_mean,
                               score = score_mean,
                               information = information_mean,
                               mle = theta_hat,
                               on = on_mean,
                               on_gradient = on_mean_gradient,
                               on_hessian = on_mean_hessian_counted,
                               level = level)
expect_equal(ci_mean_jacobian, expected_mean, tolerance = 1e-08,
             check.attributes = FALSE)
expect_true(hessian_calls > 0L)

hessian_calls <- 0L
ci_mean_no_jacobian <- profile_ci(loglik = loglik_mean,
                                  score = score_mean,
                                  information = information_mean,
                                  mle = theta_hat,
                                  on = on_mean,
                                  on_gradient = on_mean_gradient,
                                  level = level)
expect_equal(ci_mean_no_jacobian, expected_mean, tolerance = 1e-08,
             check.attributes = FALSE)
expect_equal(hessian_calls, 0L)

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
expect_error(profile_ci(loglik = loglik_mean,
                        score = score_mean,
                        information = information_mean,
                        mle = theta_hat,
                        on = on_mean,
                        on_gradient = on_mean_gradient,
                        on_hessian = function(theta) matrix(1, 2, 2),
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
                    on_hessian = on_mu_hessian,
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

on_coefficient_of_variation <- function(theta)
    exp(theta[2]) / theta[1]
expect_error(
    profile_ci(loglik = loglik_normal,
               mle = theta2_hat,
               on = on_coefficient_of_variation,
               level = 0.998,
               do_checks = FALSE),
    pattern = "Could not identify profile endpoints on opposite sides"
)

on_standardized_mean <- function(theta) theta[1] / exp(theta[2])
on_standardized_mean_gradient <- function(theta) {
    c(exp(-theta[2]), -theta[1] * exp(-theta[2]))
}
on_standardized_mean_hessian <- function(theta) {
    out <- matrix(0, 2, 2)
    out[1, 2] <- out[2, 1] <- -exp(-theta[2])
    out[2, 2] <- theta[1] * exp(-theta[2])
    out
}

ci_standardized_mean <- profile_ci(
    loglik = loglik_normal,
    score = score_normal,
    information = information_normal,
    mle = theta2_hat,
    on = on_standardized_mean,
    on_gradient = on_standardized_mean_gradient,
    on_hessian = on_standardized_mean_hessian,
    level = level
)

## Under the constraint mu / sigma = psi, maximization over sigma has a
## closed-form solution after reparameterizing by t = 1 / sigma.
profile_loglik_standardized_mean <- function(psi) {
    sy <- sum(y)
    sy2 <- sum(y^2)
    t_hat <- (psi * sy + sqrt(psi^2 * sy^2 + 4 * n * sy2)) / (2 * sy2)
    n * log(t_hat) - 0.5 * sum((y * t_hat - psi)^2)
}
psi_hat <- on_standardized_mean(theta2_hat)
lr_equation <- function(psi) {
    2 * (loglik_normal(theta2_hat) - profile_loglik_standardized_mean(psi)) - cutoff
}
find_lr_endpoint <- function(direction) {
    step <- 1
    outer <- psi_hat + direction * step
    while (lr_equation(outer) < 0) {
        step <- 2 * step
        outer <- psi_hat + direction * step
    }
    uniroot(lr_equation, sort(c(psi_hat, outer)), tol = 1e-10)$root
}
expected_standardized_mean <- c(lower = find_lr_endpoint(-1),
                                upper = find_lr_endpoint(1))

expect_equal(ci_standardized_mean, expected_standardized_mean,
             tolerance = 1e-07, check.attributes = FALSE)

profile_standardized_mean <- profile_focus(
    loglik = loglik_normal,
    score = score_normal,
    information = information_normal,
    mle = theta2_hat,
    on = on_standardized_mean,
    on_gradient = on_standardized_mean_gradient,
    on_hessian = on_standardized_mean_hessian,
    grid_size = 5,
    max_level = level
)
expected_profile_loglik <- vapply(
    profile_standardized_mean$psi,
    profile_loglik_standardized_mean,
    numeric(1)
)
expect_equal(profile_standardized_mean$loglik, expected_profile_loglik,
             tolerance = 1e-07)
expect_equal(profile_standardized_mean$psi[6], psi_hat,
             tolerance = 1e-08, check.attributes = FALSE)


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

second_parameter_hessian <- function(theta) {
    matrix(0, length(theta), length(theta))
}
budworm_focus_fit <- glm(cbind(numalive, numdead) ~ sex * ldose,
                         family = binomial, data = budworm,
                         method = brglm2::brglmFit, type = "ML")
budworm_focus <- focus(budworm_focus_fit,
                       on = second_parameter,
                       on_gradient = second_parameter_gradient,
                       on_hessian = second_parameter_hessian,
                       correction = "no")
grid_size <- 10
max_level <- 0.95
focus_profile <- profile(budworm_focus,
                         grid_size = grid_size,
                         max_level = max_level)
r_target <- qnorm(0.5 + max_level / 2)
r_step <- r_target / (grid_size - 1)
r_grid <- seq(r_step, r_target + r_step, by = r_step)

expect_true(inherits(focus_profile, "profile_focus_list"))
expect_equal(nrow(focus_profile), 2 * grid_size + 1)
expect_equal(focus_profile$signed, c(rev(r_grid), 0, -r_grid))
expect_equal(focus_profile$psi[grid_size + 1],
             second_parameter(coef(budworm_focus_fit)), tolerance = 1e-08,
             check.attributes = FALSE)
expect_equal(attr(focus_profile, "max_level"), max_level)
