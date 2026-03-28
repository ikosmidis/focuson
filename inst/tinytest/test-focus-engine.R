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
    V = V,
    on = function(theta) exp(theta[1]),
    correction = "mean",
    P = P,
    Q = Q
)
focus_out <- focus(
    coalition_fit,
    on = function(theta) exp(theta[1]),
    correction = "mean"
)

expect_equal(engine_out$estimate, focus_out$estimate, tolerance = 1e-8, check.attributes = FALSE)
expect_equal(engine_out$se, focus_out$se, tolerance = 1e-8)
expect_equal(engine_out$confint, focus_out$confint, tolerance = 1e-8, check.attributes = FALSE)

coalition_mean <- update(coalition_fit, type = "AS_mean")
afuns_mean <- enrichwith::get_auxiliary_functions(coalition_mean)
theta_mean <- coef(coalition_mean, model = "full")
V_mean <- vcov(coalition_mean, model = "full")
P_mean <- afuns_mean$Pmat()
Q_mean <- afuns_mean$Qmat()

engine_out_mean <- focus_engine(
    theta = theta_mean,
    V = V_mean,
    on = function(theta) theta[1] - theta[4],
    correction = "mean",
    estimator = "meanBR",
    P = P_mean,
    Q = Q_mean
)
engine_out_mean_no_PQ <- focus_engine(
    theta = theta_mean,
    V = V_mean,
    on = function(theta) theta[1] - theta[4],
    correction = "mean",
    estimator = "meanBR"
)
focus_out_mean <- focus(
    coalition_mean,
    on = function(theta) theta[1] - theta[4],
    correction = "mean"
)

expect_equal(engine_out_mean$estimate, focus_out_mean$estimate, tolerance = 1e-8, check.attributes = FALSE)
expect_equal(engine_out_mean$se, focus_out_mean$se, tolerance = 1e-8)
expect_equal(engine_out_mean$confint, focus_out_mean$confint, tolerance = 1e-8, check.attributes = FALSE)
expect_equal(engine_out_mean_no_PQ$estimate, focus_out_mean$estimate, tolerance = 1e-8, check.attributes = FALSE)
expect_equal(engine_out_mean_no_PQ$se, focus_out_mean$se, tolerance = 1e-8)
expect_equal(engine_out_mean_no_PQ$confint, focus_out_mean$confint, tolerance = 1e-8, check.attributes = FALSE)
