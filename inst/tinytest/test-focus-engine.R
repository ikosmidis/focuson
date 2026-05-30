library("brglm2")

data("coalition", package = "brglm2")
coalition_fit <- glm(duration ~ fract + numst2,
                     family = Gamma,
                     data = coalition,
                     method = "brglmFit",
                     type = "ML")

afuns <- enrichwith::get_auxiliary_functions(coalition_fit)
theta <- coef(coalition_fit, model = "full")
V <- vcov(coalition_fit, model = "full")
P <- afuns$Pmat()
Q <- afuns$Qmat()

engine_out <- focus_engine(
    theta = theta,
    components = list(V = V, P = P, Q = Q),
    on = function(theta) exp(theta[1]),
    correction = "mean"
)
expect_identical(engine_out$correction, "mean")
expect_identical(unname(coef(engine_out)), unname(engine_out$estimate))
expect_equal(drop(vcov(engine_out)), engine_out$se^2)
expect_true(inherits(engine_out, "focus_engine_list"))
expect_true(inherits(engine_out, "focus_list"))
expect_true(is.language(engine_out$call))
expect_equal(engine_out$theta, theta, check.attributes = FALSE)
expect_identical(engine_out$components, list(V = V, P = P, Q = Q))
expect_identical(engine_out$estimator, "ML")
expect_true(is.function(engine_out$on$on))
expect_null(engine_out$on$on_gradient)
expect_null(engine_out$on$on_hessian)
expect_identical(engine_out$dots, list())
focus_out <- focus(
    coalition_fit,
    on = function(theta) exp(theta[1]),
    correction = "mean"
)

expect_equal(engine_out$estimate, focus_out$estimate, tolerance = 1e-8, check.attributes = FALSE)
expect_equal(engine_out$se, focus_out$se, tolerance = 1e-8)
expect_equal(
    unname(confint(engine_out)),
    unname(engine_out$estimate) + c(-1, 1) * qnorm(0.975) * engine_out$se,
    check.attributes = FALSE
)
expect_error(confint(engine_out, method = "hulc"))

coalition_mean <- update(coalition_fit, type = "AS_mean")
afuns_mean <- enrichwith::get_auxiliary_functions(coalition_mean)
theta_mean <- coef(coalition_mean, model = "full")
V_mean <- vcov(coalition_mean, model = "full")
P_mean <- afuns_mean$Pmat()
Q_mean <- afuns_mean$Qmat()

engine_out_mean <- focus_engine(
    theta = theta_mean,
    components = list(V = V_mean, P = P_mean, Q = Q_mean),
    on = function(theta) theta[1] - theta[4],
    correction = "mean",
    estimator = "meanBR"
)
engine_out_mean_no_PQ <- focus_engine(
    theta = theta_mean,
    components = list(V = V_mean),
    on = function(theta) theta[1] - theta[4],
    correction = "mean",
    estimator = "meanBR"
)
expect_equal(engine_out_mean$theta, theta_mean, check.attributes = FALSE)
expect_identical(engine_out_mean$components, list(V = V_mean, P = P_mean, Q = Q_mean))
expect_identical(engine_out_mean$estimator, "meanBR")
expect_identical(engine_out_mean_no_PQ$components, list(V = V_mean))
focus_out_mean <- focus(
    coalition_mean,
    on = function(theta) theta[1] - theta[4],
    correction = "mean"
)

expect_equal(engine_out_mean$estimate, focus_out_mean$estimate, tolerance = 1e-8, check.attributes = FALSE)
expect_equal(engine_out_mean$se, focus_out_mean$se, tolerance = 1e-8)
expect_equal(engine_out_mean_no_PQ$estimate, focus_out_mean$estimate, tolerance = 1e-8, check.attributes = FALSE)
expect_equal(engine_out_mean_no_PQ$se, focus_out_mean$se, tolerance = 1e-8)


## Quantile of lm
## Note that glm is in dispersion parameterization so theta = sigma^2
quant <- function(theta, x0, p) {
    theta[1] + theta[2] * x0 + qnorm(p) * sqrt(theta[3])
}
quant_grad <- function(theta, x0, p) {
    c(1, x0, 0.5 * qnorm(p) / sqrt(theta[3]))
}
quant_hess <- function(theta, x0, p) {
    out <- matrix(0, 3, 3)
    out[3, 3] <- - 0.25 * qnorm(p) / theta[3]^(3/2)
    out
}

nobs <- 20
theta_true <- c(1, 2, 1.5^2)
set.seed(123)
x <- seq(0, 1, length.out = nobs)
df <- data.frame(x = x,
                 y = theta_true[1] + theta_true[2] * x + rnorm(nobs, mean = 0, sd = sqrt(theta_true[3])))
temp_mod <- glm(y ~ x, data = df)

c_mod <- do.call(update, list(object = temp_mod, data = df))

ff <- focus(c_mod, on = quant, correction = "median", x0 = 1, p = 0.95)
ff_analytical <- focus(c_mod, on = quant,
                       on_gradient = quant_grad,
                       on_hessian = quant_hess,
                       correction = "median", x0 = 1, p = 0.95)

expect_equal(ff$estimate, ff_analytical$estimate)
expect_equal(ff$se, ff_analytical$se)

## Using focus_engine
Q <- function(theta, x) {
    mat <- matrix(0, 3, 3)
    Q1 <- Q2 <- Q3 <- mat
    Q1[2, 3] <- Q1[3, 2] <- sum(x)
    Q1[1, 3] <- Q1[3, 1] <- length(x)
    Q1 <- - Q1 / theta[3]^2
    Q2[2, 3] <- Q2[3, 2] <- sum(x^2)
    Q2[1, 3] <- Q2[3, 1] <- sum(x)
    Q2 <- - Q2 / theta[3]^2
    Q3[3, 3] <- -length(x) / theta[3]^3
    list(Q1, Q2, Q3)
}
P <- function(theta, x) {
    mat <- matrix(0, 3, 3)
    P1 <- P2 <- P3 <- mat
    P1[2, 3] <- P1[3, 2] <- sum(x)
    P1[1, 3] <- P1[3, 1] <- length(x)
    P1 <- P1 / theta[3]^2
    P2[2, 3] <- P2[3, 2] <- sum(x^2)
    P2[1, 3] <- P2[3, 1] <- sum(x)
    P2 <- P2 / theta[3]^2
    P3[1, 1:2] <- P3[1:2, 1] <- c(length(x), sum(x)) / theta[3]^2
    P3[2, 2] <- sum(x^2) / theta[3]^2
    P3[3, 3] <- length(x) / theta[3]^3
    list(P1, P2, P3)
}
info <- function(theta, x) {
    X <- cbind(1, x)
    xx <- crossprod(X)
    rbind(cbind(xx / theta[3], 0),
          c(0, 0, length(x) / (2 * theta[3]^2)))
}

c_mod_ml <- update(c_mod, method = "brglmFit", type = "ML")
afuns <- enrichwith::get_auxiliary_functions(c_mod_ml)
cml <- coef(c_mod_ml, "full")
ff_engine <- focus_engine(
    cml,
    components = list(V = solve(info(cml, df$x)), P = P(cml, df$x), Q = Q(cml, df$x)),
    on = quant,
    correction = "median",
    x0 = 1,
    p = 0.95
)
expect_identical(ff_engine$correction, "median")
expect_identical(ff_engine$dots, list(x0 = 1, p = 0.95))
expect_equal(ff$estimate, ff_engine$estimate, check.attributes = FALSE)
expect_equal(ff$se, ff_engine$se, check.attributes = FALSE)


## Check afun implementation against Qm
expect_equal(Q(cml, df$x)[[1]], afuns$Qmat()[[1]], check.attributes = FALSE)
expect_equal(Q(cml, df$x)[[2]], afuns$Qmat()[[2]], check.attributes = FALSE)
expect_equal(Q(cml, df$x)[[3]], afuns$Qmat()[[3]], check.attributes = FALSE)

## Check afun implementation against Pm
expect_equal(P(cml, df$x)[[1]], afuns$Pmat()[[1]], check.attributes = FALSE)
expect_equal(P(cml, df$x)[[2]], afuns$Pmat()[[2]], check.attributes = FALSE)
expect_equal(P(cml, df$x)[[3]], afuns$Pmat()[[3]], check.attributes = FALSE)


## estimate_focus_components() diagnostics
norm_loglik <- function(theta, data, ...) {
    sum(dnorm(data, mean = theta[1], sd = 1, log = TRUE))
}

norm_score <- function(theta, data, ...) {
    sum(data - theta[1])
}

norm_info <- function(theta, data, ...) {
    matrix(length(data), 1, 1)
}

sim_norm <- function(theta, n, ...) {
    rnorm(n, mean = theta[1], sd = 1)
}

set.seed(1)
comp_diag <- estimate_focus_components(
    theta = 0,
    loglik = norm_loglik,
    score = norm_score,
    information = norm_info,
    simulate = sim_norm,
    nsim = 6,
    diagnostics = TRUE,
    n = 5
)

expect_true(inherits(comp_diag, "focus_components"))
expect_true("diagnostics" %in% names(comp_diag))
expect_true(all(c("I", "P", "Q") %in% names(comp_diag$diagnostics)))
expect_true(all(c("mcse_max", "mcse_frobenius", "rel_mcse_frobenius",
                  "kappa", "min_eigen") %in% names(comp_diag$diagnostics$I)))
expect_identical(length(comp_diag$diagnostics$P), 1L)
expect_identical(length(comp_diag$diagnostics$Q), 1L)

printed_comp_diag <- capture.output(print(comp_diag))
expect_true(any(grepl("^Monte Carlo component estimates$", printed_comp_diag)))
expect_true(any(grepl("^Parameters: 1\\s*$", printed_comp_diag)))
expect_true(any(grepl("^Simulations: 6\\s*$", printed_comp_diag)))
expect_true(any(grepl("^Diagnostics: yes\\s*$", printed_comp_diag)))
expect_true(any(grepl("^Diagnostics$", printed_comp_diag)))

## estimate_focus_components_fef()
obs_norm <- c(-1, 0, 1, 2, 3)

set.seed(3)
comp_diag_fef <- estimate_focus_components_fef(
    theta = 0,
    data = obs_norm,
    loglik = norm_loglik,
    score = norm_score,
    information = norm_info,
    simulate = sim_norm,
    nsim = 6,
    diagnostics = TRUE,
    n = 5
)

expect_true(inherits(comp_diag_fef, "focus_components"))
expect_true(isTRUE(comp_diag_fef$meta$fef))
expect_equal(drop(comp_diag_fef$V), 1 / length(obs_norm))
expect_equal(drop(comp_diag_fef$Q[[1]]), 0)
printed_comp_diag_fef <- capture.output(print(comp_diag_fef))
expect_true(any(grepl("^Structure: full exponential family shortcut\\s*$", printed_comp_diag_fef)))

## estimate_focus_components_iid()
sim_norm_one <- function(theta, ...) {
    rnorm(1, mean = theta[1], sd = 1)
}

set.seed(2)
comp_diag_iid <- estimate_focus_components_iid(
    theta = 0,
    n = 5,
    loglik = norm_loglik,
    score = norm_score,
    information = norm_info,
    simulate = sim_norm_one,
    nsim = 400,
    diagnostics = TRUE
)

expect_true(inherits(comp_diag_iid, "focus_components"))
expect_true(isTRUE(comp_diag_iid$meta$iid))
expect_identical(comp_diag_iid$meta$n, 5L)
expect_equal(drop(comp_diag_iid$V), 1 / 5, tolerance = 0.15)
printed_comp_diag_iid <- capture.output(print(comp_diag_iid))
expect_true(any(grepl("^Sampling: iid contributions\\s*$", printed_comp_diag_iid)))
expect_true(any(grepl("^Target sample size: 5\\s*$", printed_comp_diag_iid)))
