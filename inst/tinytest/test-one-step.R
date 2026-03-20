## A Gamma example
data("coalition", package = "brglm2")
coalitionML <- glm(duration ~ fract + numst2, family = Gamma, data = coalition,
                   method = "brglmFit", type = "ML")

expect_warning(
fit1s <- update(coalitionML, maxit = 1,  max_step_factor = 1, type = "AS_median",
                start = coef(coalitionML))
)

coefs_medBR <- coef(fit1s, model = "full")
expect_equal(sapply(1:4, function(k) focus(coalitionML, on = function(theta) theta[k])),
             coefs_medBR)

## A binomial logit example
data("lizards", package = "brglm2")
lizardsML <- glm(cbind(grahami, opalinus) ~ height + diameter +
                     light + time, family = binomial(), data = lizards,
                 method = "brglmFit", type = "ML")
expect_warning(
    fit1s <- update(lizardsML, maxit = 1,  max_step_factor = 1, type = "AS_median",
                    start = coef(lizardsML))
)
coefs_medBR <- coef(fit1s)
expect_equal(sapply(1:6, function(k) focus(lizardsML, on = function(theta) theta[k])),
             coefs_medBR)

## A binomial cauchit example
data("lizards", package = "brglm2")
lizardsML <- glm(cbind(grahami, opalinus) ~ height + diameter +
                     light + time, family = binomial(cauchit), data = lizards,
                 method = "brglmFit", type = "ML")
expect_warning(
    fit1s <- update(lizardsML, maxit = 1,  max_step_factor = 1, type = "AS_median",
                    start = coef(lizardsML))
)
coefs_medBR <- coef(fit1s)
expect_equal(sapply(1:6, function(k) focus(lizardsML, on = function(theta) theta[k])),
             coefs_medBR, tolerance = 1e-06)
