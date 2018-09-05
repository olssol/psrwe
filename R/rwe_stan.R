#' Call STAN models
#'
#'
#' @param chains STAN parameter. Number of Markov chainsm
#' @param iter STAN parameter. Number of iterations
#' @param warmup STAN parameter. Number of burnin.
#' @param control STAN parameter. See \code{rstan::stan} for details.
#' @param ... other options to call STAN sampling such as \code{thin},
#'     \code{algorithm}. See \code{rstan::sampling} for details.#'
#'
#'
#' @export
#'
rweSTAN <- function(lst.data, stan.mdl = "powerp",
                    chains = 4, iter = 2000, warmup = 1000,
                    control = list(adapt_delta=0.95), ...) {

    stan.rst <- rstan::sampling(stanmodels[[stan.mdl]],
                                data    = lst.data,
                                chains  = chains,
                                iter    = iter,
                                warmup  = warmup,
                                control = control,
                                ...);

    stan.rst;
}



#' Get Posterior for all stratum
#'
#' @param data class DWITHPS data frame
#' @param type type of outcomes
#' @param A    target number of subjects to be borrowed
#' @param RS   parameters for dirichelet prior
#' @param Fix.RS whether treat RS as fixed or the prior mean of vs
#' @param ...  extra parameters for calling function \code{\link{rweSTAN}}
#'
#' @export
#'
rwePsPowDrawPost <- function(data, type = c("continuous", "binary"), A = 0, RS = NULL, Fix.RS = FALSE,
                             v.outcome = "Y",  ...) {

    stopifnot(v.outcome %in% colnames(data));
    type <- match.arg(type);

    ## prepare stan data
    data   <- data[!is.na(data[["_strata_"]]),];
    S      <- max(data[["_strata_"]]);
    stan.d <- NULL;

    Y1     <- NULL;
    INX1   <- NULL;
    for (i in 1:S) {
        cur.d1 <- data[data[["_strata_"]] == i & data[["_grp_"]] == 1, v.outcome];
        cur.d0 <- data[data[["_strata_"]] == i & data[["_grp_"]] == 0, v.outcome];

        if (0 == length(cur.d1) | 0 == length(cur.d0)) {
            stop(paste("Stratum ", i, " contains no subjects from group 1", sep = ""));
        }

        cur.n1 <- length(cur.d1);
        cur.d  <- c(N0 = length(cur.d0), YBAR0 = mean(cur.d0), SD0   = sd(cur.d0),
                    N1 = cur.n1,         YBAR1 = mean(cur.d1), YSUM1 = sum(cur.d1));

        stan.d <- rbind(stan.d, cur.d);
        Y1     <- c(Y1, cur.d1);
        INX1   <- c(INX1, rep(i, length = cur.n1));
    }

    if (is.null(RS))
        RS <- rep(1/S, S);

    lst.data  <- list(S     = S,
                      A     = A,
                      RS    = RS,
                      FIXVS = as.numeric(Fix.RS),
                      N0    = stan.d[,"N0"],
                      N1    = stan.d[,"N1"],
                      YBAR0 = stan.d[,"YBAR0"],
                      SD0   = stan.d[,"SD0"]);

    ## sampling
    if ("continuous" == type) {
        stan.mdl  <- ifelse(1 == S,  "powerp", "powerps");
        lst.data <- c(lst.data,
                      list(TN1  = length(Y1),
                           Y1   = Y1,
                           INX1 = INX1));
        rst.post  <- rweSTAN(lst.data = lst.data, stan.mdl = stan.mdl, ...);
        rst.theta <- rstan::extract(rst.post, pars = "theta")$theta;
    } else {
        lst.data <- c(lst.data,
                      list(YBAR1 = as.numeric(stan.d[,"YBAR1"]),
                           YSUM1 = as.numeric(stan.d[,"YSUM1"])));

        if (1 < S) {
            rst.post  <- rweSTAN(lst.data = lst.data, stan.mdl = "powerpsbinary", ...);
            rst.theta <- rstan::extract(rst.post, pars = "theta")$theta;
        } else {
            rst.post  <- NULL;
            rst.theta <- with(lst.data, rbeta(2000, YSUM1+A*YBAR0+1, N1-YSUM1+A*(1-YBAR0)+1));
        }
    }

    ## return
    list(post.theta = rst.theta,
         stan.rst   = rst.post);
}

#' Summary Posterior theta
#'
#'
#' @param post.theta posterior samples from STAN
#' @param true.theta true value of theta
#' @param quants     quantiles
#' @param weights    number of subjects in each strata as the weight
#'
#' @export
#'
rweSummaryPost <- function(post.theta, true.theta, quants = c(0.025, 0.975), weights=1) {

    ##cur.post      <- apply(post.theta, 2, function(x) {mean(weights*x)});
    cur.post        <- post.theta;
    cur.m           <- mean(cur.post);
    cur.var         <- var(cur.post);
    cur.ci          <- quantile(cur.post, quants);
    ##post.var        <- apply(post.theta, 1, var);
    ##names(post.var) <- paste("var", 1:length(post.var), sep = "");
    range.ci        <- range(cur.ci);

    rst <- c(thetahat = cur.m,
             thetavar = cur.var,
             bias     = cur.m - true.theta,
             mse      = (cur.m - true.theta)^2,
             width    = range.ci[2] - range.ci[1],
             cover    = true.theta >= range.ci[1] & true.theta <= range.ci[2],
             lb       = as.numeric(cur.ci[1]),
             ub       = as.numeric(cur.ci[2]));
    rst

}

#' Get Posterior for all stratum
#'
#'
#' @param data.ps class DWITHPS data frame
#' @param A    power of the power prior for each strata
#' @param ...  extra parameters for calling function \code{\link{rweSTAN}}
#'
#' @export
#'
rwePowerDrawPost <- function(data.ps, v.outcome = "Y", As=0, type = c("continuous", "binary"), ...) {
    stopifnot("RWE_DWITHPS" %in% class(data.ps));

    data <- data.ps$data;
    stopifnot(v.outcome %in% colnames(data));
    type     <- match.arg(type);
    stan.mdl <- switch(type,
                       continuous = "powerp",
                       binary     = "powerpbinary");

    data      <- data[!is.na(data[["_strata_"]]),];
    nstrata   <- max(data[["_strata_"]]);
    rst.theta <- NULL;
    for (i in 1:nstrata) {
        cur.d1 <- data[data[["_strata_"]] == i & data[["_grp_"]] == 1, v.outcome];
        cur.d0 <- data[data[["_strata_"]] == i & data[["_grp_"]] == 0, v.outcome];

        if (0 == length(cur.d1)) {
            warning("Stratum contains no subjects from group 1");
            next;
        }

        Y1 <- cur.d1;
        N1 <- length(cur.d1);
        N0 <- length(cur.d0);
        if (0 == N0) {
            ##dummy value
            Y0 <- 0;
        } else {
            Y0 <- cur.d0;
        }

        cur.A     <- As[min(i, length(As))];
        cur.post  <- rweSTAN(lst.data = list(A=cur.A, N0=N0, N1=N1, Y0=Y0, Y1=Y1),
                             stan.mdl = stan.mdl, ...);
        cur.theta <- rstan::extract(cur.post, pars = "theta")$theta;
        rst.theta <- rbind(rst.theta, cur.theta);
    }

    rst.theta
}
