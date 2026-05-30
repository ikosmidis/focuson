library("brglm2")

data("coalition", package = "brglm2")
coalition_fit <- glm(duration ~ fract + numst2,
                     family = Gamma,
                     data = coalition,
                     method = "brglmFit",
                     type = "ML")

focus_square <- focus(coalition_fit,
                      on = function(theta, k) theta[k]^2,
                      correction = "median",
                      k = 2)
theta_square <- focus_se(focus_square)

expect_true(abs(unname(theta_square[2]^2) - unname(coef(focus_square))) < 1e-8)

all_square <- focus_on_all(focus_square)

expect_equal(rownames(all_square), c("estimate", "se"))
expect_equal(colnames(all_square), names(coef(focus_square$object, model = "full")))
expect_equal(length(theta_square), ncol(all_square))

focus_exp <- focus(coalition_fit,
                   on = function(theta) exp(theta[2]),
                   correction = "median")

expect_error(focus_se(focus_exp),
             pattern = "Could not reconstruct the reference parameter vector")
