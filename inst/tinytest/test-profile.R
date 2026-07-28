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
profile_signed <- function(object) {
    sign(attr(object, "mle") - object$psi) *
        sqrt(2 * pmax(attr(object, "max_loglik") - object$loglik, 0))
}

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
             theta_hat - profile_signed(profile_mean) / sqrt(n),
             tolerance = 1e-08)
expect_identical(names(profile_mean),
                 c("psi", "loglik", "theta", "lagrange"))
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

focus_range_mean <- theta_hat + c(-0.5, 0.5)
profile_mean_focus_grid <- profile_focus(
    mle = theta_hat,
    loglik = loglik_mean,
    score = score_mean,
    on = on_mean,
    on_gradient = on_mean_gradient,
    approach = "focus_grid",
    focus_range = focus_range_mean,
    grid_size = 4,
    auglag_args = list(control.outer = list(eps = 1e-08, e.scale = 2))
)
expected_focus_grid <- c(
    seq(focus_range_mean[1], theta_hat, length.out = 5),
    seq(theta_hat, focus_range_mean[2], length.out = 5)[-1]
)
expect_true(inherits(profile_mean_focus_grid, "profile_focus_list"))
expect_identical(attr(profile_mean_focus_grid, "approach"), "focus_grid")
expect_equal(attr(profile_mean_focus_grid, "focus_range"),
             focus_range_mean)
expect_equal(profile_mean_focus_grid$psi, expected_focus_grid,
             tolerance = 1e-06)
expect_equal(profile_mean_focus_grid$psi,
             theta_hat - profile_signed(profile_mean_focus_grid) / sqrt(n),
             tolerance = 1e-06)
expect_equal(profile_mean_focus_grid$loglik,
             vapply(profile_mean_focus_grid$psi, loglik_mean, numeric(1)),
             tolerance = 1e-07)
expect_equal(as.numeric(profile_mean_focus_grid$theta[, 1]),
             profile_mean_focus_grid$psi,
             tolerance = 1e-06)
expected_focus_grid_lagrange <- vapply(
    profile_mean_focus_grid$theta[, 1],
    score_mean,
    numeric(1)
)
expect_equal(profile_mean_focus_grid$lagrange,
             expected_focus_grid_lagrange,
             tolerance = 1e-03)

expect_error(
    profile_focus(mle = theta_hat,
                  loglik = loglik_mean,
                  approach = "focus_grid"),
    pattern = "two finite numeric values"
)
expect_error(
    profile_focus(mle = theta_hat,
                  loglik = loglik_mean,
                  approach = "focus_grid",
                  focus_range = rep(theta_hat, 2)),
    pattern = "two distinct values"
)
profile_mean_right <- profile_focus(
    mle = theta_hat,
    loglik = loglik_mean,
    score = score_mean,
    on = on_mean,
    on_gradient = on_mean_gradient,
    approach = "focus_grid",
    focus_range = theta_hat + c(0.25, 0.5),
    grid_size = 4
)
expect_equal(profile_mean_right$psi,
             seq(theta_hat + 0.25, theta_hat + 0.5, length.out = 5),
             tolerance = 1e-06)
expect_true(all(profile_signed(profile_mean_right) < 0))
expect_equal(attr(profile_mean_right, "mle"), theta_hat)
expect_equal(attr(profile_mean_right, "max_loglik"),
             loglik_mean(theta_hat))
expect_equal(profile_mean_right$loglik,
             vapply(profile_mean_right$psi, loglik_mean, numeric(1)),
             tolerance = 1e-07)
expect_equal(profile_mean_right$lagrange,
             vapply(profile_mean_right$psi, score_mean, numeric(1)),
             tolerance = 1e-03)

profile_mean_left <- profile_focus(
    mle = theta_hat,
    loglik = loglik_mean,
    score = score_mean,
    on = on_mean,
    on_gradient = on_mean_gradient,
    approach = "focus_grid",
    focus_range = theta_hat - c(0.5, 0.25),
    grid_size = 4
)
expect_equal(profile_mean_left$psi,
             seq(theta_hat - 0.5, theta_hat - 0.25, length.out = 5),
             tolerance = 1e-06)
expect_true(all(profile_signed(profile_mean_left) > 0))

profile_mean_endpoint <- profile_focus(
    mle = theta_hat,
    loglik = loglik_mean,
    score = score_mean,
    on = on_mean,
    on_gradient = on_mean_gradient,
    approach = "focus_grid",
    focus_range = theta_hat + c(0, 0.5),
    grid_size = 4
)
expect_equal(profile_mean_endpoint$psi,
             seq(theta_hat, theta_hat + 0.5, length.out = 5),
             tolerance = 1e-06)
expect_equal(profile_mean_endpoint$lagrange[1L], 0)
profile_mean_ignored_arguments <- profile_focus(
    mle = theta_hat,
    loglik = loglik_mean,
    score = score_mean,
    information = information_mean,
    on = on_mean,
    on_gradient = on_mean_gradient,
    on_hessian = on_mean_hessian,
    grid_size = 5,
    max_level = level,
    focus_range = focus_range_mean,
    auglag_args = "ignored"
)
expect_equal(
    profile_mean_ignored_arguments,
    profile_mean,
    tolerance = 1e-08,
    check.attributes = FALSE
)
expect_equal(attr(profile_mean_ignored_arguments, "focus_range"),
             focus_range_mean)

profile_mean_focus_grid_ignored_arguments <- profile_focus(
    mle = theta_hat,
    loglik = loglik_mean,
    score = score_mean,
    on = on_mean,
    on_gradient = on_mean_gradient,
    approach = "focus_grid",
    focus_range = focus_range_mean,
    grid_size = 4,
    max_level = "ignored",
    nleqslv_args = "ignored"
)
expect_equal(
    profile_mean_focus_grid_ignored_arguments$psi,
    profile_mean_focus_grid$psi,
    tolerance = 1e-06
)
expect_identical(attr(profile_mean_focus_grid_ignored_arguments, "max_level"),
                 "ignored")

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
expect_true(any(grepl("^Profiling approach: VM\\s*$",
                      printed_profile_mean)))
expect_true(any(grepl("^Parameter dimension: 1\\s*$",
                      printed_profile_mean)))
expect_true(any(grepl("^Profile points: 11\\s*$",
                      printed_profile_mean)))
expect_true(any(grepl("^Points by side: 5 left, 5 right\\s*$",
                      printed_profile_mean)))
expect_false(any(grepl("^Signed-root range:",
                       printed_profile_mean)))
expect_false(any(grepl("^Median signed-root spacing:",
                       printed_profile_mean)))
expect_true(any(grepl("^Requested maximum nominal level: 0.95\\s*$",
                      printed_profile_mean)))
expect_true(any(grepl("^Implied level at grid boundary:",
                      printed_profile_mean)))

printed_profile_mean_right <- capture.output(print(profile_mean_right,
                                                   digits = 6))
expect_true(any(grepl("^Points by side: 0 left, 5 right\\s*$",
                      printed_profile_mean_right)))
expected_boundary_right <- pchisq(
    2 * (attr(profile_mean_right, "max_loglik") -
         tail(profile_mean_right$loglik, 1)),
    df = 1
)
expect_true(any(grepl(
    paste0("^Implied level at grid boundary: ",
           format(signif(expected_boundary_right, 6), trim = TRUE),
           "\\s*$"),
    printed_profile_mean_right
)))

profile_ci_mean <- confint(profile_mean, level = 0.9)
expected_profile_ci_mean <- theta_hat +
    c(lower = -1, upper = 1) * qnorm(0.95) / sqrt(n)
expect_equal(profile_ci_mean, expected_profile_ci_mean,
             tolerance = 1e-08, check.attributes = FALSE)
expect_equal(attr(profile_ci_mean, "level"), 0.9)
expect_identical(attr(profile_ci_mean, "type"), "profile")
expect_identical(attr(profile_ci_mean, "interpolation"), "linear")
expect_equal(confint(profile_mean, level = 0.9,
                     interpolation = "cubic"),
             expected_profile_ci_mean,
             tolerance = 1e-08, check.attributes = FALSE)
expect_error(confint(profile_mean, level = 0),
             pattern = "number in \\(0, 1\\)")
expect_error(confint(profile_mean, interpolation = "quadratic"))
expect_error(confint(profile_mean, level = 0.999999),
             pattern = "exceeds the range")
expect_error(confint(profile_mean_right),
             pattern = "exceeds the range")

profile_plot_file <- tempfile(fileext = ".pdf")
pdf(profile_plot_file)
expect_silent(plot(profile_mean, interpolation = "linear"))
expect_silent(plot(profile_mean, interpolation = "cubic"))
expect_silent(plot(profile_mean_focus_grid, interpolation = "linear"))
expect_silent(plot(profile_mean_focus_grid, signed = TRUE))
expect_silent(plot(profile_mean, ci = TRUE))
expect_silent(plot(profile_mean, signed = TRUE, ci = TRUE))
expect_silent(plot(profile_mean, level = 0.999999))
expect_silent(plot(profile_mean, signed = TRUE, level = 0.999999))
expect_silent(plot(profile_mean_right))
expect_silent(plot(profile_mean_right, signed = TRUE))
dev.off()
unlink(profile_plot_file)
expect_error(plot(profile_mean, interpolation = "quadratic"))
expect_error(plot(profile_mean, ci = TRUE, level = 0.999999),
             pattern = "exceeds the range")
expect_error(plot(profile_mean_right, ci = TRUE),
             pattern = "exceeds the range")

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
             theta_hat - profile_signed(profile_mean_defaults) / sqrt(n),
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

profile_standardized_mean_focus_grid <- profile_focus(
    loglik = loglik_normal,
    score = score_normal,
    mle = theta2_hat,
    on = on_standardized_mean,
    on_gradient = on_standardized_mean_gradient,
    approach = "focus_grid",
    focus_range = expected_standardized_mean + c(-0.25, 0.25),
    grid_size = 5
)
expect_equal(
    profile_standardized_mean_focus_grid$loglik,
    vapply(profile_standardized_mean_focus_grid$psi,
           profile_loglik_standardized_mean, numeric(1)),
    tolerance = 1e-06
)
expect_true(max(abs(
    apply(profile_standardized_mean_focus_grid$theta, 1,
          on_standardized_mean) -
        profile_standardized_mean_focus_grid$psi
)) < 1e-06)


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
expect_equal(profile_signed(focus_profile),
             c(rev(r_grid), 0, -r_grid),
             tolerance = 1e-06)
expect_equal(focus_profile$psi[grid_size + 1],
             second_parameter(coef(budworm_focus_fit)), tolerance = 1e-08,
             check.attributes = FALSE)
expect_equal(attr(focus_profile, "max_level"), max_level)

budworm_focus_mle <- second_parameter(coef(budworm_focus_fit))
focus_profile_grid <- profile(
    budworm_focus,
    approach = "focus_grid",
    focus_range = budworm_focus_mle + c(-0.5, 0.5),
    grid_size = 3
)
expect_true(inherits(focus_profile_grid, "profile_focus_list"))
expect_identical(attr(focus_profile_grid, "approach"), "focus_grid")
expect_equal(
    focus_profile_grid$psi,
    c(seq(budworm_focus_mle - 0.5, budworm_focus_mle, length.out = 4),
      seq(budworm_focus_mle, budworm_focus_mle + 0.5,
          length.out = 4)[-1]),
    tolerance = 1e-06
)
