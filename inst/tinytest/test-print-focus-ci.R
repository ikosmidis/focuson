wald <- focuson:::.focus_ci(structure(
    c(lower = -1, upper = 1),
    level = 0.95,
    type = "wald",
    se_at = "compatible"
))
printed_wald <- capture.output(print(wald))
expect_true(any(grepl("95% Wald confidence interval", printed_wald,
                      fixed = TRUE)))
expect_true(any(grepl("Standard error evaluated at: compatible", printed_wald,
                      fixed = TRUE)))

vm <- focuson:::.focus_ci(structure(
    c(lower = -1, upper = 1),
    level = 0.9,
    type = "pl",
    iter = c(lower = 3L, upper = 4L),
    `max|fvec|` = c(lower = 1e-10, upper = 2e-10)
))
printed_vm <- capture.output(print(vm))
expect_true(any(grepl("90% profile likelihood confidence interval", printed_vm,
                      fixed = TRUE)))
expect_true(any(grepl("Endpoint iterations: lower = 3, upper = 4", printed_vm,
                      fixed = TRUE)))
expect_true(any(grepl("Maximum equation residual:", printed_vm,
                      fixed = TRUE)))

mpl <- focuson:::.focus_ci(structure(
    c(lower = -1, upper = 1),
    level = 0.99,
    type = "mpl",
    interpolation = "cubic"
))
printed_mpl <- capture.output(print(mpl))
expect_true(any(grepl("99% modified profile likelihood confidence interval",
                      printed_mpl, fixed = TRUE)))
expect_true(any(grepl("Interpolation: cubic", printed_mpl, fixed = TRUE)))

hulc <- focuson:::.focus_ci(structure(
    c(lower = NA_real_, upper = NA_real_),
    level = 0.95,
    type = "hulc",
    B = 5L,
    Delta = 0,
    error = "statistic failed"
))
printed_hulc <- capture.output(print(hulc))
expect_true(any(grepl("95% HulC confidence interval", printed_hulc,
                      fixed = TRUE)))
expect_true(any(grepl("Partitions: 5", printed_hulc, fixed = TRUE)))
expect_true(any(grepl("Median-bias bound (Delta): 0", printed_hulc,
                      fixed = TRUE)))
expect_true(any(grepl("Stored error: statistic failed", printed_hulc,
                      fixed = TRUE)))
