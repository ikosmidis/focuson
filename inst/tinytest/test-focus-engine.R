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
