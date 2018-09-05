#' Simulate covariates following a mixture of multivariate normal distribution
#'
#' @inheritParams simupara
#'
#' @export
#'
rweSimuCov <- function(nPat, muCov, sdCov, corCov, mix.phi = 1, seed = NULL, cov.breaks = NULL) {

    f.cur <- function(x, i) {
        if (is.array(x)) {
            rst <- x[min(i, nrow(x)),];
        } else {
            rst <- x;
        }

        rst
    }

    stopifnot(is.numeric(mix.phi) | any(mix.phi < 0));
    stopifnot(nPat > 0);

    if (!is.null(seed))
        set.seed(seed);

    n.pts <- rmultinom(1, nPat, mix.phi);
    cov.x <- NULL;
    for (i in 1:length(mix.phi)) {

        if (0 == n.pts[i])
            next;

        cur.mu  <- f.cur(muCov, i);
        cur.sd  <- f.cur(sdCov, i);
        cur.cor <- corCov[min(i, length(corCov))];
        cur.x   <- rmvnorm(n.pts[i],
                           mean=cur.mu,
                           sigma=get.covmat(cur.sd, cur.cor));
        cov.x   <- rbind(cov.x, cur.x);
    }

    colnames(cov.x) <- paste("V", 1:ncol(cov.x), sep="");
    cov.x           <- data.frame(cov.x);
    cov.x           <- get.cov.cat(cov.x, cov.breaks);
}

#' Simulate X multiplied by Beta
#'
#' @inheritParams simupara
#' @param ... Parameters for simulating covariates by function
#'     \code{\link{rweSimuCov}}
#'
#'
#' @export
#'
rweXBeta <- function(..., regCoeff, cov.x = NULL, fmla = NULL) {

    stopifnot(inherits(fmla,"formula") | is.null(fmla));

    if (is.null(cov.x))
        cov.x <- rweSimuCov(...);

    if (is.null(fmla)) {
        fmla <- formula(paste("~",
                              paste(colnames(cov.x), collapse = "+")));
    }

    d.matrix <- model.matrix(fmla, cov.x);
    xbeta    <- get.xbeta(d.matrix, regCoeff);
    xbeta
}


#' Compute standard error of the random error
#'
#' @inheritParams simupara
#' @inheritParams rweGetYSig
#'
#' @return mean of xbeta and standard error or the random error term
#'
#' @export
#'
rweGetYSig <- function(..., nPat=500000, xbeta = NULL, sig2Ratio = 1) {
    if (is.null(xbeta))
        xbeta   <- rweXBeta(nPat=nPat, ...);

    v.xbeta <- var(xbeta);
    ysig    <- sqrt(v.xbeta * sig2Ratio);

    c(mean(xbeta), ysig);
}

#' Get intercept for a binary outcome.
#'
#' The binary outcome may be an outcome or a treatment assignment.
#'
#' @inheritParams simupara
#'
#' @param ... Parameters for simulating covariates by function
#'     \code{\link{rweXBeta}} and \code{\link{rweSimuCov}}
#'
#' @return standard error or the random error term
#'
#' @export
#'
rweGetBinInt <- function(..., regCoeff, nPat=500000, xbeta = NULL, bin.mu = 0.5) {
    ## fill in 0 for intercept temporarily
    if (is.null(xbeta))
        ey <- rweXBeta(nPat, regCoeff = c(0, regCoeff), ...);

    fx <- function(b0) {
        logp <- (b0 + ey) - log(1 + exp(b0+ey));
        m    <- mean(exp(logp));
        abs(m - bin.mu);
    }

    mey <- max(abs(ey));
    rst <- optimize(fx, c(-100-max(ey),100+max(ey)));
    rst$minimum
}

#' Simulate random errors
#'
#' @inheritParams simupara
#'
#' @export
#'
rweSimuError <- function(nPat,
                         error.type = c("normal", "skewed"),
                         ysig   = 1,
                         skew.n = NULL,skew.p = NULL, skew.noise = 0.0001,
                         ...) {

    type <- match.arg(error.type);
    rst <- switch(type,
                  normal = {rnorm(nPat, 0, ysig)},
                  skewed = {
        mu        <- skew.n * (1-skew.p) / skew.p;
        va        <- skew.n * (1-skew.p) / skew.p^2;
        noise.sig <- skew.noise;
        rst       <- rnbinom(nPat, skew.n, skew.p);
        rst       <- rst - mu + rnorm(nPat, 0, noise.sig);
        rst       <- rst/sqrt(va + noise.sig^2)*ysig;
    });

    rst
}

#' Simulate outcomes for a single arm study
#'
#' @inheritParams simupara
#'
#' @export
#'
rweSimuSingleArm <- function(nPat, muCov, sdCov, corCov, regCoeff, mix.phi = 1, cov.breaks = NULL,
                             fmla = NULL,
                             type = c("continuous", "binary"),
                             ysig = NULL, sig2Ratio=1, b0 = NULL, bin.mu = 0.5,
                             ...) {

    type  <- match.arg(type);
    COV.X <- rweSimuCov(nPat       = nPat,
                        muCov      = muCov,
                        sdCov      = sdCov,
                        corCov     = corCov,
                        mix.phi    = mix.phi,
                        cov.breaks = cov.breaks);

    ##simulate Y
    if ("continuous" == type) {
        ## epsilon
        if (is.null(ysig)) {
            ysig <- rweGetYSig(muCov      = muCov,
                               sdCov      = sdCov,
                               corCov     = corCov,
                               mix.phi    = mix.phi,
                               cov.breaks = cov.breaks,
                               regCoeff   = regCoeff,
                               sig2Ratio  = sig2Ratio,
                               fmla       = fmla)[2];
        }

        XBETA   <- rweXBeta(cov.x = COV.X, regCoeff = regCoeff, fmla = fmla);
        EPSILON <- rweSimuError(nPat, ysig = ysig,...);
        Y       <- XBETA + EPSILON;
    } else if ("binary" == type) {
        if (is.null(b0)) {
            stopifnot(!is.null(bin.mu));
            b0  <- rweGetBinInt(bin.mu,
                                muCov      = muCov,
                                sdCov      = sdCov,
                                corCov     = corCov,
                                regCoeff   = regCoeff,
                                mix.phi    = mix.phi,
                                cov.breaks = cov.breaks,
                                fmla       = fmla);
        }
        regCoeff <- c(b0, regCoeff);
        XBETA    <- rweXBeta(cov.x = COV.X, regCoeff = regCoeff, fmla = fmla);
        Y        <- rbinom(nPat, 1, expit(XBETA));
    }

    ##return
    Data           <- cbind(1:nPat, Y, COV.X);
    colnames(Data) <- c("pid", "Y", paste("V", 1:ncol(COV.X), sep=""));
    data.frame(Data);
}

#' Simulate continuous outcomes for a two arm study
#'
#' @inheritParams simupara
#'
#' @export
#'
rweSimuTwoArm <- function(nPat, muCov, sdCov, corCov, trt.effect = 0,
                          regCoeff.y, regCoeff.z=0,
                          mix.phi = 1, cov.breaks = NULL,
                          fmla.y = NULL, fmla.z = NULL, ysig = NULL,
                          b0 = NULL, z1.p = 0.5,
                          sig2Ratio = 2,  ..., do.simu=TRUE) {
    ##treatment assignment
    if (is.null(b0) & !identical(0, regCoeff.z)) {
        b0  <- rweGetBinInt(z1.p,
                            muCov    = muCov,
                            sdCov    = sdCov,
                            corCov   = corCov,
                            regCoeff = regCoeff.z,
                            mix.phi  = mix.phi,
                            fmla     = fmla.z);
    }

    if (is.null(ysig)) {
        ysig <- rweGetYSig(muCov      = muCov,
                           sdCov      = sdCov,
                           corCov     = corCov,
                           mix.phi    = mix.phi,
                           cov.breaks = cov.breaks,
                           regCoeff   = regCoeff.y,
                           sig2Ratio  = sig2Ratio,
                           fmla       = fmla.y)[2];
    }

    simu.data <- NULL;
    if (do.simu) {
        ##covariates
        COV.X   <- rweSimuCov(nPat       = nPat,
                              muCov      = muCov,
                              sdCov      = sdCov,
                              corCov     = corCov,
                              mix.phi    = mix.phi,
                              cov.breaks = cov.breaks);

        if (identical(0, regCoeff.z)) {
            Z  <- rbinom(nPat, 1, z1.p);
        } else {
            xbeta.z <- rweXBeta(cov.x    = COV.X,
                                regCoeff = regCoeff.z,
                                fmla     = fmla.z);
            Z       <- rbinom(nPat, 1, expit(b0 + xbeta.z));
        }

        xbeta.y <- rweXBeta(cov.x    = COV.X,
                            regCoeff = regCoeff.y,
                            fmla     = fmla.y);

        epsilon.y <- rweSimuError(nPat, ysig = ysig, ...);
        Y         <- Z * trt.effect + xbeta.y + epsilon.y;
        simu.data <- cbind(pid=1:nPat, Y=Y, Z=Z, COV.X);
    }

    list(true.effect = trt.effect,
         simu.data   = simu.data,
         b0ysig      = c(b0, ysig));
}

#' Simulate data from an existing dataset
#'
#'
#' @inheritParams simupara
#'
#' @export
#'
rweSimuFromTrial <- function(nPat, trial.data, group = "A", outcome = "Y",
                             with.replacement = TRUE, seed = NULL,
                             permute = TRUE, f.subset = NULL,
                             permute.trteffect = 0, permute.interaction = 0,
                             simu.group = group, simu.outcome = outcome) {
    if (!is.null(seed))
        set.seed(seed);

    if (1 == length(nPat)) {
        ## set the same sample size for the two groups
        nPat <- rep(nPat,2);
    }

    if (!permute) {
        arms   <- unique(trial.data[[group]]);
        rst    <- NULL;
        mean.y <- NULL;
        for (i in 1:length(arms)) {
            cur.d   <- trial.data[arms[i] == trial.data[[group]],];
            cur.n   <- nPat[min(i, length(nPat))];

            stopifnot(with.replacement | nrow(cur.d) > cur.n);

            mean.y   <- c(mean.y, mean(cur.d[[outcome]]));
            cur.smp  <- cur.d[sample(1:nrow(cur.d), cur.n, replace=with.replacement), ];

            cur.smp[[simu.group]]   <- i - 1;
            cur.smp[[simu.outcome]] <- cur.d[[outcome]];

            rst <- rbind(rst, cur.smp);
        }

        trt.effect <- mean.y;
        simu.data  <- rst;
    } else {
        cur.d <- trial.data;
        cur.n <- sum(nPat);

        stopifnot(with.replacement | nrow(cur.d) > cur.n);

        smp.inx <- sample(1:nrow(cur.d), cur.n, replace=with.replacement);
        cur.smp <- cur.d[smp.inx, ];
        grps    <- NULL;
        for (i in 1:length(nPat)) {
            grps <- c(grps, rep(i-1, nPat[i]));
        }
        cur.smp[[simu.group]]   <- grps;
        cur.smp[[simu.outcome]] <- cur.smp[[outcome]];

        ##introduce main effect to the last group
        inx.last <- which(max(grps) == grps);
        cur.smp[inx.last, simu.outcome] <- permute.trteffect + cur.smp[inx.last, simu.outcome];

        ##introduce interaction effect
        if (is.function(f.subset) & permute.interaction != 0) {
            all.subgrp <- f.subset(cur.d);
            stopifnot(all(all.subgrp %in% c(0,1)));

            egbi   <- mean(all.subgrp);
            subgrp <- all.subgrp[smp.inx];
            cur.smp[inx.last, simu.outcome] <- cur.smp[inx.last, simu.outcome] +
                permute.interaction * (subgrp[inx.last] - egbi);
        }

        trt.effect <- permute.trteffect;
        simu.data  <- cur.smp;
    }

    ## randomize the order of simulated patients
    list(true.effect = trt.effect,
         simu.data   = simu.data[sample(1:nrow(simu.data)),]);
}


#' Combining simulation results
#'
#' @param lst.rst List of simulation results. Each element represents a
#'     replication
#' @export
#'
rweSimuCombine <- function(lst.rst, fun = mean, ignore.error = TRUE, ...) {

    if (ignore.error) {
        err.inx <- NULL;
        for (i in 1:length(lst.rst)) {
            if ("try-error" == class(lst.rst[[i]]))
                err.inx <- c(err.inx, i);
        }

        if (!is.null(err.inx))
            lst.rst <- lst.rst[-err.inx];
    }

    nreps       <- length(lst.rst);
    rep1        <- lst.rst[[1]];
    lst.combine <- rep(list(NULL), length(rep1));
    for (i in 1:nreps) {
        for (j in 1:length(rep1)) {
            cur.value        <- lst.rst[[i]][[j]];
            lst.combine[[j]] <- rbind(lst.combine[[j]],
                                      as.vector(as.matrix(cur.value)));
        }
    }

    for (j in 1:length(lst.combine)) {
        cur.rst           <- apply(lst.combine[[j]], 2, fun, ...);
        dim(cur.rst)      <- dim(as.matrix(rep1[[j]]));
        dimnames(cur.rst) <- dimnames(as.matrix(rep1[[j]]));
        lst.combine[[j]]  <- cur.rst;
    }

    names(lst.combine) <- names(rep1);
    lst.combine
}
