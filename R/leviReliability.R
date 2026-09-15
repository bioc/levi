# Shared by the observed landscape and every node-label permutation.
# Missing measurements stay undefined until the signal transform, where they
# become neutral. They are recorded separately in result$metadata.
.networkSignals <- function(node_values, edge_index, single_col,
                            signal_mode, logfc_k) {
    edges <- (node_values[edge_index[, 1], , drop = FALSE] +
              node_values[edge_index[, 2], , drop = FALSE]) / 2
    values <- rbind(node_values, edges)
    test <- matrix(values[, 1], ncol = 1)
    control <- matrix(if (single_col) rep(1, nrow(values)) else values[, 2],
                      ncol = 1)
    missing <- !is.finite(test) | !is.finite(control)
    # Do not warn about missing measurements as if they were a zero ratio.
    test[missing] <- if (signal_mode == "ratio") 1 else NA_real_
    control[missing] <- if (signal_mode == "ratio") 1 else NA_real_
    signal <- .computeSignalOut(test, control, signal_mode = signal_mode,
                                logfc_k = logfc_k, single_col = single_col)
    signal[missing] <- 0.5
    # The displayed and tested quantity is always m1. Auxiliary channels are
    # retained for the C++ interface, with finite values on the same scale.
    list(signal = signal, test = signal, control = signal)
}

.adjustLandscapePvalues <- function(pvalues, method) {
    # Both directional families are adjusted together, over occupied cells.
    values <- c(pvalues$over, pvalues$under)
    keep <- is.finite(values)
    values[keep] <- stats::p.adjust(values[keep], method = method)
    n <- length(pvalues$over)
    list(over = matrix(values[seq_len(n)], nrow(pvalues$over)),
         under = matrix(values[n + seq_len(n)], nrow(pvalues$under)))
}

.signalMeaning <- function(mode, single_col) {
    if (mode == "zscore") return("0.5 = mean logFC of network support points")
    if (mode == "ratio" && single_col)
        return("abundance/(abundance + 1); no control comparison")
    "0.5 = no change"
}

.validateScalar <- function(x, name, lower = -Inf, upper = Inf,
                            integer = FALSE) {
    if (!is.numeric(x) || length(x) != 1L || !is.finite(x) ||
        x < lower || x > upper || (integer && x != floor(x)))
        stop("'", name, "' must be a finite ",
             if (integer) "integer" else "number", " in [", lower,
             ", ", upper, "].", call. = FALSE)
}
