library("brglm2")

data("endometrial", package = "brglm2")
endo <- glm(HG ~ NV + PI + EH,
            data = endometrial,
            family = binomial("logit"),
            method = "brglmFit")

set.seed(678)
out <- focus(
    endo,
    on = function(theta, randomize = 2) randomize * theta[1],
    ci = "hulc",
    randomize = 2,
    control_ci = ci_control(check_statistic = FALSE)
)

expect_identical(out$ci_type, "hulc")
expect_true(is.numeric(out$estimate))
expect_true(is.numeric(out$confint))
expect_identical(attr(out$confint, "type"), "hulc")


set.seed(678)
out <- focus(
    endo,
    on = function(theta, check_statistic = 2) check_statistic * theta[1],
    correction = "no",
    ci = "hulc",
    check_statistic = 2,
    control_ci = ci_control(check_statistic = FALSE)
)

expect_identical(out$ci_type, "hulc")
expect_true(is.numeric(out$estimate))
expect_identical(attr(out$confint, "type"), "hulc")
expect_equal(
    unname(out$estimate),
    2 * unname(focus(endo, on = function(theta) theta[1], correction = "no")$estimate)
)


for (correction in c("no", "mean", "median")) {
    out <- focus(endo, on = function(theta) pi, correction = correction)

    expect_equal(unname(out$estimate), pi, check.attributes = FALSE)
    expect_equal(out$se, 0)
    expect_equal(unname(out$confint), c(pi, pi), check.attributes = FALSE)

    set.seed(678)
    out_hulc <- focus(
        endo,
        on = function(theta) pi,
        correction = correction,
        ci = "hulc",
        control_ci = ci_control(check_statistic = FALSE)
    )

    expect_equal(unname(out_hulc$estimate), pi, check.attributes = FALSE)
    expect_equal(out_hulc$se, 0)
    expect_equal(unname(out_hulc$confint), c(pi, pi), check.attributes = FALSE)
    expect_identical(attr(out_hulc$confint, "type"), "hulc")
}


expect_error(focus(endo, alpha = 0))
expect_error(focus(endo, alpha = 1))
expect_error(focus(endo, alpha = 2))


f_base <- focus(endo, on = function(theta) theta[1], correction = "median")$estimate
f_scaled <- focus(endo,  on = function(theta) 1e-9 * theta[1],  correction = "median")$estimate

expect_equal(
    unname(f_scaled),
    1e-9 * unname(f_base),
    tolerance = 1e-12,
    check.attributes = FALSE
)


grad_or <- function(theta) {
    out <- rep(0, length(theta))
    out[1] <- exp(theta[1])
    out
}

hessian_or <- function(theta) {
    out <- matrix(0, nrow = length(theta), ncol = length(theta))
    out[1, 1] <- exp(theta[1])
    out
}

out_num <- focus(
    endo,
    on = function(theta) exp(theta[1]),
    correction = "median"
)
out_ana <- focus(
    endo,
    on = function(theta) exp(theta[1]),
    correction = "median",
    on_gradient = grad_or,
    on_hessian = hessian_or
)

expect_equal(out_ana$estimate, out_num$estimate, tolerance = 1e-8, check.attributes = FALSE)
expect_equal(out_ana$se, out_num$se, tolerance = 1e-8)
expect_equal(out_ana$confint, out_num$confint, tolerance = 1e-8, check.attributes = FALSE)

set.seed(678)
out_hulc <- focus(
    endo,
    on = function(theta) exp(theta[1]),
    correction = "median",
    ci = "hulc",
    on_gradient = grad_or,
    on_hessian = hessian_or,
    control_ci = ci_control(check_statistic = FALSE)
)

expect_identical(out_hulc$ci_type, "hulc")
expect_true(is.numeric(out_hulc$estimate))
expect_identical(attr(out_hulc$confint, "type"), "hulc")


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
true_quant <- quant(theta_true, x0 = 1, p = 0.95)

c_mod <- do.call(update, list(object = temp_mod, data = df))
ff <- focus(c_mod, on = quant, correction = "median", x0 = 1, p = 0.95)
ff_analytical <- focus(c_mod, on = quant,
                       on_gradient = quant_grad,
                       on_hessian = quant_hess,
                       correction = "median", x0 = 1, p = 0.95)

expect_equal(ff, ff_analytical)

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
afuns <- get_auxiliary_functions(c_mod_ml)
cml <- coef(c_mod_ml, "full")
ff_engine <- focus_engine(cml, V = solve(info(cml, df$x)), on = quant, correction = "median",
                    P = P(cml, df$x), Q = Q(cml, df$x), x0 = 1, p = 0.95)
expect_equal(ff$estimate, ff_engine$estimate, check.attributes = FALSE)
expect_equal(ff, ff_engine, check.attributes = FALSE)


## Check afun implementation against Qm
expect_equal(Q(cml, df$x)[[1]], afuns$Qmat()[[1]], check.attributes = FALSE)
expect_equal(Q(cml, df$x)[[2]], afuns$Qmat()[[2]], check.attributes = FALSE)
expect_equal(Q(cml, df$x)[[3]], afuns$Qmat()[[3]], check.attributes = FALSE)

## Check afun implementation against Pm
expect_equal(P(cml, df$x)[[1]], afuns$Pmat()[[1]], check.attributes = FALSE)
expect_equal(P(cml, df$x)[[2]], afuns$Pmat()[[2]], check.attributes = FALSE)
expect_equal(P(cml, df$x)[[3]], afuns$Pmat()[[3]], check.attributes = FALSE)
