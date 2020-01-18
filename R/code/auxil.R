wgt_mse <- function(n, y, yhat) {
    sum(n * (y - yhat)^2, na.rm = TRUE) / sum(n > 0)
}

calc_point_criterion <- function(pred_occ, pred_pot, n, y, mx, obj_fun = wgt_mse) {
    crit <- matrix(0, ncol(pred_occ), ncol(pred_pot))
    y[is.na(y)] <- 0
    y[y > mx] <- mx
    for(i in seq_len(nrow(crit))) {
        for(j in seq_len(ncol(crit))) {
            tmp <- pred_occ[ , i] * pred_pot[ , j]
            tmp[tmp > mx] <- mx
            crit[i, j] <- obj_fun(n, y, tmp)
        }}
    return(crit)
}

calc_cov_criterion <- function(draws_logpot, sig2, data, n_draw = 250, seed = 1, type_pot = 'arith', scale = 1, size = 0.90) {
    set.seed(seed)

    cov <- length <- loglength <- rep(0, dim(draws_logpot)[[2]])
    N <- nrow(data)

    for(j in seq_len(length(cov))) {
        tmp <- matrix(0, N, n_draw)
        for(k in seq_len(n_draw)) {
            ## variance of potential scales inversely with number of occupied points
            if(type_pot == 'arith') {
                ypot <- rnorm(N, exp(draws_logpot[ , j, k]),
                                        sqrt(sig2[ , j]/(data$points_occ/scale)))
            } else {
                ypot <- rnorm(N, draws_logpot[ , j, k],
                                            sqrt(sig2[ , j]/(data$points_occ/scale)))
                ypot <- exp(ypot)
            }
            tmp[ , k] <- ypot
        }
        qq <- apply(tmp, 1, quantile, c((1-size)/2, 1-(1-size)/2), na.rm = TRUE)
        cov[j] <- mean(data$obs <= qq[2, ] & data$obs >= qq[1, ])
        length[j] <- median(qq[2, ] - qq[1, ])
        loglength[j] <- median(log(qq[2, ]) - log(qq[1, ]))
    }
    names(cov) <- names(length) <- names(loglength) <- names(draws_logpot)
    return(list(cov = cov, length = length, loglength = loglength))
}
