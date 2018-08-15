wgt_mse <- function(n, y, yhat) {
    sum(n * (y - yhat)^2, na.rm = TRUE) / sum(n > 0)
}

calc_cv_criterion <- function(pred_occ, pred_pot, n, y, mx, obj_fun = wgt_mse) {
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
