data("coalition", package = "brglm2")
c_ml <- glm(duration ~ fract + numst2, family = Gamma, data = coalition,
            method = "brglmFit", type = "ML")
c_mean <- update(c_ml, type = "AS_mean")
c_median <- update(c_ml, type = "AS_median")
c_mixed <- update(c_mean, type = "AS_mixed")
c_cor <- update(c_mean, type = "correction")
c_jef <- update(c_ml, type = "MPL_Jeffreys")

## does not refit
fo_mean <- focus(c_mean, on = function(theta) theta[1] - theta[4], correction = "mean")
fo_cor <- focus(c_cor, on = function(theta) theta[1] - theta[4], correction = "mean")
## refits with AS_mean
fo_mixed <- focus(c_mixed, on = function(theta) theta[1] - theta[4], correction = "mean")
fo_median <- focus(c_median, on = function(theta) theta[1] - theta[4], correction = "mean")
fo_jef <- focus(c_jef, on = function(theta) theta[1] - theta[4], correction = "mean")

## Coeficients from AS_mean
co_mean <- coef(c_mean, model = "full")
co_cor <- coef(c_cor, model = "full")
co_ml <- coef(c_ml, model = "full")

expect_equal(co_mean[1] - co_mean[4], fo_mean, check.attributes = FALSE)
expect_equal(co_mean[1] - co_mean[4], fo_mixed, check.attributes = FALSE)
expect_equal(co_mean[1] - co_mean[4], fo_median, check.attributes = FALSE)
expect_equal(co_cor[1] - co_cor[4], fo_cor, check.attributes = FALSE)

## No correction
fo_ml <- focus(c_ml, on = function(theta) theta[1] - theta[4], correction = "no")
expect_equal(co_ml[1] - co_ml[4], fo_ml, check.attributes = FALSE)

## Check standard errors
fo_mean <- focus(c_mean, on = function(theta) theta[2], correction = "mean")
fo_median <- focus(c_median, on = function(theta) theta[2], correction = "mean")
fo_cor <- focus(c_cor, on = function(theta) theta[2], correction = "mean")
expect_equal(attr(fo_cor, "se"), coef(summary(c_cor))[2, "Std. Error"])
expect_equal(attr(fo_mean, "se"), coef(summary(c_mean))[2, "Std. Error"])
expect_equal(attr(fo_median, "se"), coef(summary(c_mean))[2, "Std. Error"])





