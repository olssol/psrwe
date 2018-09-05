#' Cut a sequence of numbers into bins with equal numbers in each bin
#'
#'
#' @export
#'
rweCut <- function(x, y=x, breaks = 5) {
    cuts    <- quantile(x, seq(0, 1,length=breaks+1));
    cuts[1] <- cuts[1] - 0.001;
    rst     <- rep(NA, length(y));
    for (i in 2:length(cuts)) {
        inx      <- which(y > cuts[i-1] & y <= cuts[i]);
        rst[inx] <- i-1;
    }

    rst
}

#' Get unbalance in baseline covariates
#'
#' @param diff If TRUE, get the difference in covariates between groups.
#'     Otherwise, get covariates for each group separately
#'
#' @param var.group Column name in pts that corresponds to treat group
#'
#' @inheritParams simupara
#'
#' @return If diff is TRUE, return a dataset with columns V and Diff. Otherwise,
#'     return a dataset with columns V, Z and Value. In both cases, column V
#'     represents the name of covariates.
#'
#' @export
#'
rweUnbalance <- function(nPat, ..., pts = NULL, covs = NULL, diff = TRUE,
                         var.group = "Z", cov.pattern = "^V[0-9]+$") {
    if (is.null(pts)) {
        pts <- rweSimuTwoArm(nPat, ...);
    }

    ##unbalance
    inx.0 <- which(0 == pts[, var.group]);
    if (is.null(covs)) {
        c.xy  <- colnames(pts);
        c.xy  <- c.xy[grep(cov.pattern, c.xy)];
    } else {
        c.xy <- covs;
    }

    unb   <- NULL;
    for (i in 1:length(c.xy)) {
        x0     <- sample(pts[inx.0,  c.xy[i]], size = nPat, replace = TRUE);
        x1     <- sample(pts[-inx.0, c.xy[i]], size = nPat, replace = TRUE);

        if (diff) {
            x.diff <- x1 - x0;
            unb    <- rbind(unb, data.frame(V=c.xy[i], Diff=x1 - x0));
        } else {
            unb    <- rbind(unb, data.frame(V=c.xy[i], Z=1, Value=x1));
            unb    <- rbind(unb, data.frame(V=c.xy[i], Z=0, Value=x0));
        }
    }

    ## make group column factor if it exists
    if (!diff) {
        unb$Z <- as.factor(unb$Z)
    }

    unb
}


#' Compute distance from F0 to F1
#'
#' @param type type of distances. ovl: overlapping coefficient, kl:
#'     1/(1+Kullback-Leibler divergence)
#' @param n.bins number of bins for KL computation
#' @param epsilon small integer for Dirichlet smoothing
#'
#' @return a vector with the number of samples in group 0, the number of samples
#'     in group 1, and 1/(1+KL divergence) from group 0 to group 1 when type is
#'     kl, or the overlapping coefficient when type is ovl
#' @export
#'
rweDist <- function(sample.F0, sample.F1, n.bins = 10, type = c("kl", "ovl"), epsilon = 10^-6) {

    type     <- match.arg(type);

    smps     <- c(sample.F0, sample.F1);
    n0       <- length(sample.F0);
    n1       <- length(sample.F1);

    if (1 == length(unique(smps))) {
        cut.smps <- rep(1, n0+n1)
        n.bins   <- 1;
        warning("Distributions for computing distances are degenerate.",
                call. = FALSE);
    } else {
        cut.smps <- rweCut(smps, breaks = n.bins);
    }

    rst <- 0;
    for (j in 1:n.bins) {
        n0.j <- length(which(j == cut.smps[1:n0]));
        n1.j <- length(which(j == cut.smps[(n0+1):(n0+n1)]));

        rst  <- rst + switch(type,
                             kl = {ep0  <- (n0.j+epsilon)/(n0 + epsilon * n.bins);
                          ep1  <- (n1.j+epsilon)/(n1 + epsilon * n.bins);
                          ep1 * log(ep1/ep0)},
                          ovl = min(n0.j/n0, n1.j/n1));
    }

    if ("kl" == type)
        rst <- 1/(1+rst);

    c(n0,n1,rst);
}


#' Generate frequency table for factor columns
#'
#' @return a vector with the number of samples in group 0, the number of samples
#'     in group 1, and the KL divergence from group 0 to group 1
#' @export
#'
rweFreqTbl <- function(data, var.groupby, vars = NULL) {

    if (is.null(vars))
        vars <- colnames(data);

    rst <- NULL;
    for (v in vars) {
        if (!is.factor(data[[v]]))
            next;

        cur.freq <- data %>% count_(c(var.groupby, v)) %>%
            group_by_(.dots = var.groupby) %>%
            mutate(Sum = sum(n), Freq = n/sum(n)) %>%
            mutate_if(is.factor, as.character) %>%
            mutate(Cov = v) %>%
            rename_(Value = v);

        rst <- rbind(rst, data.frame(cur.freq));
    }

    rst
}
