#' Get propensity scores
#'
#' @param ... parameters to get propensity scores
#'
#' @export
#'
rwePS <- function(data, ps.fml = NULL,
                  v.grp = "group",
                  v.covs = "V1",
                  d1.grp = 1,
                  nstrata = 5, ...,
                  caliper = 0) {

    dnames <- colnames(data);
    stopifnot(v.grp %in% dnames);

    ## generate formula
    if (is.null(ps.fml))
        ps.fml <- as.formula(paste(v.grp, "~", paste(v.covs, collapse="+"), sep=""));

    all.ps  <- get.ps(data, ps.fml = ps.fml, ...);
    D1.ps   <- all.ps[which(d1.grp == data[[v.grp]])];

    ## add caliper width
    if (caliper > 0) {
        lgt.d1.ps  <- log(D1.ps/(1-D1.ps));
        sd.lgt     <- sd(lgt.d1.ps[!is.infinite(lgt.d1.ps)]);
        D1.ps[which.min(D1.ps)] <- expit(min(lgt.d1.ps) - sd.lgt * caliper);
        D1.ps[which.max(D1.ps)] <- expit(max(lgt.d1.ps) + sd.lgt * caliper);
    }

    ## stratification
    strata  <- rweCut(D1.ps, all.ps, breaks = nstrata);
    grp     <- rep(1, nrow(data));
    grp[which(data[[v.grp]] != d1.grp)] <- 0;

    data[["_ps_"]]     <- all.ps;
    data[["_strata_"]] <- strata;
    data[["_grp_"]]    <- grp;


    rst <- list(data    = data,
                ps.fml  = ps.fml,
                nstrata = nstrata);

    class(rst) <- get.rwe.class("DWITHPS");

    rst
}

#' Get number of subjects and the distances of PS distributions for each PS
#' strata
#'
#'
#' @export
#'
rwePSDist <- function(data.withps, n.bins = 10, type = c("ovl", "kl"), ...) {
    stopifnot(inherits(data.withps,
                       what = get.rwe.class("DWITHPS")));

    type <- match.arg(type);

    dataps   <- data.withps$data;
    nstrata  <- data.withps$nstrata;
    rst     <- NULL;
    for (i in 1:nstrata) {
        ps0 <- dataps[which(i == dataps[["_strata_"]] &
                                 0 == dataps[["_grp_"]]),
                           "_ps_"];
        ps1 <- dataps[which(i == dataps[["_strata_"]] &
                                 1 == dataps[["_grp_"]]),
                           "_ps_"];

        if (0 == length(ps0) | 0 == length(ps1))
            warning("No samples in strata");

        if (any(is.na(c(ps0, ps1))))
            warning("NA found in propensity scores in a strata");

        cur.dist <- rweDist(ps0, ps1, n.bins = n.bins, type = type, ...);
        rst      <- rbind(rst, c(i, cur.dist));
    }

    ##overall
    ps0        <- dataps[which(0 == dataps[["_grp_"]]), "_ps_"];
    ps1        <- dataps[which(1 == dataps[["_grp_"]]), "_ps_"];
    all.dist   <- rweDist(ps0, ps1, n.bins = nstrata*n.bins, type = type, ...);
    rst        <- rbind(rst, c(0, all.dist));


    colnames(rst) <- c("Strata", "N0", "N1", "Dist");
    rst           <- data.frame(rst);
    class(rst)    <- append(get.rwe.class("PSDIST"), class(rst));

    rst
}

#' Get the actual power term in the power prior
#'
#' @param psdist      A RWE_PSDIST type object
#' @param a           power term
#' @param overall.ess ratio of overall added number of patients to N1
#' @param adjust.size whether adjust for sizes in group 0 and 1 in the power term
#' @param adjust.dist whether adjust for distance in ps scores in the power term
#'
#' @export
#'
rweGetPowerA <- function(psdist, a = NULL, overall.ess = 0.3, adjust.size = TRUE, adjust.dist = TRUE) {

    stopifnot(inherits(psdist, what = get.rwe.class("PSDIST")));

    ## compute a
    if (is.null(a)) {
        stopifnot(1 == adjust.size);
        stopifnot(overall.ess >= 0);

        if (1 == adjust.dist) {
            a <- overall.ess / mean(psdist$Dist);
        } else {
            a <- overall.ess;
        }
    }


    ## compute as, power term for each strata
    rst <- rep(a, nrow(psdist));
    if (1 == adjust.size) {
        rst <- rst * psdist$N1 / psdist$N0;
    }

    if (1 == adjust.dist) {
        rst <- rst * psdist$Dist;
    }

    list(a  = a,
         as = rst);
}
