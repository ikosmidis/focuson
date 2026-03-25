compute_B <- function(alpha, Delta){
    stopifnot(Delta < 1/2)
    log2 <- log(2)
    loga <- log(alpha)
    Dp <- 1/2 + Delta
    Dm <- 1/2 - Delta
    logDp <- log(Dp)
    log2mloga <- log2 - loga
    ## https://doi.org/10.1093/jrsssb/qkad134 has floor for Bu but let's put
    ## ceiling here to be slightly more conservative
    Bl <- max(floor(loga / logDp), floor(log2mloga / log2))
    Bu <- ceiling(- log2mloga / logDp)
    for (B in Bl:Bu) 
        if (Dm^B + Dp^B <= alpha) break
    B
}

hulc_ci <- function(data,
                    statistic,
                    alpha = 0.05,
                    Delta = 0,
                    randomize = TRUE,
                    check_statistic = TRUE,
                    parallel = FALSE, ...) {
    nobs <- nrow(data)
    data <- data[sample(nobs), , drop = FALSE]
    B <- compute_B(alpha, Delta)
    if (randomize) {
        Dp <- 0.5 + Delta
        Dm <- 0.5 - Delta
        ## DpB <- Dp^B
        ## DmB <- Dm^B
        ## B1 <- B - 1
        ## prob0 <- Dm^B + Dp^B
        ## prob1 <- Dm^B1 + Dp^B1
        ## B <- if ((alpha - prob0) / (prob1 - prob0) < runif(1)) B else B1
        pB <- c(Dm^B, Dp^B)
        prob0 <- sum(pB)
        prob1 <- sum(pB / c(Dm, Dp))        
        B <- B - ((alpha - prob0) / (prob1 - prob0) >= runif(1))
    }
    if (B > nobs)
        stop("The required number of partitions (=", B, ") is larger than `nrow(data)` (=", nobs, ").")
    partitions <- ((seq_len(nobs) - 1) * B) %/% nobs + 1
    ## partitions <- cut(seq_len(nobs), breaks = B, labels = FALSE)
    data <- split(data, partitions)
    ## Check if statistic returns a value or fails for the smallest
    ## number of observations
    if (check_statistic) {
        min_id <- which.min(sapply(data, nrow))
        stat <- try(statistic(data[[min_id]], ...))
        test <- isTRUE(inherits(stat, "try-error") || is.na(stat))
    }
    if (test) {
        ci <- c(NA, NA)
    } else {
        if (parallel) {
            plan(multisession)
            ci <- future_sapply(data, statistic, ...) |> range()
        } else {
            ci <- sapply(data, statistic, ...) |> range()
        }
    }
    names(ci) <- c("lower", "upper")
    attr(ci, "Delta") <- Delta
    attr(ci, "alpha") <- alpha
    attr(ci, "B") <- B
    attr(ci, "type") <- "hulc"
    ci
}




