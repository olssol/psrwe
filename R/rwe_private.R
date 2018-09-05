##----------------------------------------------------------------------------
##                SETTING
##----------------------------------------------------------------------------
get.rwe.class <- function(c.str = c("DWITHPS", "PSDIST")) {
    c.str <- match.arg(c.str);
    switch(c.str,
           DWITHPS = "RWE_DWITHPS",
           PSDIST  = "RWE_PSDIST");
}


##----------------------------------------------------------------------------
##                GENERAL FUNCTIONS
##----------------------------------------------------------------------------

make.global <- function(alist, dest.env='.GlobalEnv') {
    for (i in 1:length(alist)) {
        assign(names(alist[i]), alist[[i]], dest.env );
    }
}


expit <- function(x) {
    ex <- exp(x);
    ex/(1+ex);
}

get.xbeta <- function(covX, regCoeff) {

    if (length(regCoeff) > 0 &
        length(regCoeff) != ncol(covX))
        warning("Number of coefficients does not match with the design matrix.");

    apply(covX, 1, function(x) {sum(x * regCoeff)});
}


get.covmat <- function(StDevCovar, corrCovar) {
    n.x      <- length(StDevCovar);
    Vars     <- StDevCovar*StDevCovar;
    CovarMat <- matrix(NA, n.x, n.x);
    for (i in 1:n.x) {
        CovarMat[i,i] <- Vars[i];
        for (j in i:n.x) {
            if (j == i) {
                CovarMat[i,i] <- Vars[i];
                next;
            }
            CovarMat[i, j] <- corrCovar*StDevCovar[i]*StDevCovar[j];
            CovarMat[j, i] <- CovarMat[i, j];
        }
    }

    CovarMat
}

## cut covariates into categories
get.cov.cat <- function(covX, breaks = NULL) {
    f.cut <- function(x, bs) {
        if (is.null(bs))
            return(x);

        bs  <- sort(unique(c(-Inf, bs, Inf)));
        rst <- as.numeric(cut(x, breaks = bs)) - 1;
        factor(rst);
    }

    if (is.null(breaks))
        return(covX);

    if (is.numeric(breaks)) {
        rst <- apply(covX, 1, f.cut, breaks);
        rst <- t(rst);
    } else if (is.list(breaks)) {
        rst <- covX;
        for (i in 1:min(ncol(covX), length(breaks))) {
            rst[,i] <- f.cut(covX[,i], breaks[[i]]);
        }
    }

    rst
}

##----------------------------------------------------------------------------
##                PROPENSITY SCORES
##----------------------------------------------------------------------------
##compute propensity scores
get.ps <- function(dta, ps.fml, type = c("randomforest", "logistic"), ntree = 5000,
                   ..., grp = NULL, ps.cov = NULL) {

    type <- match.arg(type);

    ## generate formula
    if (is.null(ps.fml))
        ps.fml <- as.formula(paste(grp, "~", paste(ps.cov, collapse="+"), sep=""));

    ## identify grp if passed from formula
    grp <- all.vars(ps.fml)[1];

    ## fit model
    switch(type,
           logistic = {glm.fit <- glm(ps.fml, family=binomial, data=dta, ...);
                       est.ps <- glm.fit$fitted;},
           randomforest = {dta[[grp]] <- as.factor(dta[[grp]]);
                           rf.fit     <- randomForest(ps.fml, data = dta, ntree = ntree, ...);
                           est.ps     <- predict(rf.fit, type = "prob")[,2];
                          });

    est.ps
}


##----------------------------------------------------------------------------
##                PRESENTATION
##----------------------------------------------------------------------------
## plot density of propensity score
plotRwePs <- function(data.withps, overall.inc = TRUE, add.text = TRUE, facet.scales = "free_y", ...) {
    stopifnot(inherits(data.withps,
                       what = get.rwe.class("DWITHPS")));

    pskl        <- rwePSDist(data.withps, ...);
    nstrata     <- data.withps$nstrata;
    dtaps       <- data.withps$data;

    pskl$Strata <- as.factor(c(paste("Stratum ", 1:nstrata, sep = ""), "Overall"));
    xlim        <- range(dtaps[which(!is.na(dtaps[["_strata_"]])),"_ps_"], na.rm = TRUE);

    all.data <- NULL;
    for (i in 1:nstrata) {
        cur.sub  <- dtaps[which(i == dtaps[["_strata_"]]),];
        cur.data <- data.frame(Strata = paste("Stratum ", i, sep = ""),
                               Ps     = cur.sub[["_ps_"]],
                               Group  = cur.sub[["_grp_"]]);
        all.data <- rbind(all.data, cur.data);
    }

    if (overall.inc) {
        cur.data <-  data.frame(Strata = "Overall",
                                Ps     = dtaps[["_ps_"]],
                                Group  = dtaps[["_grp_"]]);

        all.data <- rbind(all.data, cur.data);
    } else {
        pskl <- pskl %>% filter(Strata != "Overall");
    }

    all.data$Group <- as.factor(all.data$Group);
    rst <- ggplot(data = all.data, aes(x = Ps)) +
        geom_density(alpha = 0.2,
                       aes(group = Group,
                           color = Group,
                           fill  = Group,
                           linetype = Group),
                     trim  = TRUE,
                     na.rm = TRUE) +
        labs(x = "Propensity Score", y = "Density") +
        scale_y_continuous(breaks = NULL) +
        scale_x_continuous(limits = xlim) +
        theme_bw() +
        theme(strip.background = element_blank(),
              panel.grid = element_blank(),
              panel.border = element_rect(colour = "black"),
              panel.spacing = unit(0, "lines")) +
        facet_grid(Strata ~ ., scales = facet.scales);

    if (add.text) {
        rst <- rst +
            geom_text(x = Inf, y = Inf, hjust = 1, vjust = 1,
                      aes(label = paste('n0=', N0,
                                        ", n1=", N1,
                                        ", OVL=", format(Dist, digits = 3),
                                        sep = "")),
                      data = pskl, size = 4);
    }
    rst
}

plot.balance.fac <- function(dtaps, v, overall.inc = TRUE) {
    cur.d <- rweFreqTbl(dtaps, var.groupby = c("Strata", "Group"), vars = v);
    cur.d <- cur.d %>% dplyr::filter(!is.na(Strata));
    cur.d$Strata <- paste("Stratum ", cur.d$Strata, sep = "")
    if (overall.inc) {
        cur.overall <- rweFreqTbl(dtaps, var.groupby = "Group", vars = v);
        cur.overall$Strata <- "Overall";
        cur.d <- rbind(cur.d, cur.overall);
    }
    cur.d$Group <- as.factor(cur.d$Group);
    cur.d$Value <- as.factor(cur.d$Value);

    rst <- ggplot(data = cur.d, aes(x = Value, y = Freq)) +
        geom_bar(alpha = 0.4,
                 stat = "identity",
                 position = "dodge",
                 aes(group = Group,
                     linetype = Group,
                     fill  = Group)) +
        scale_y_continuous(breaks = NULL, limits = c(0,1)) +
        labs(x = "", y = "") +
        facet_grid(Strata ~ .);
    rst
}

plot.balance.cont <- function(dtaps, v, nstrata, overall.inc = TRUE, facet.scales = "free_y") {
    cur.d <- NULL;
    for (i in 1:nstrata) {
        cur.sub      <- dtaps[which(i == dtaps[["_strata_"]]),];
        cur.v        <- data.frame(Cov    = v,
                                   Value  = cur.sub[[v]],
                                   Group  = cur.sub[["_grp_"]]);
        cur.v$Strata <- paste("Stratum ", i, sep = "");
        cur.d        <- rbind(cur.d, cur.v);
    }

    if (overall.inc) {
        cur.sub      <- dtaps;
        cur.v        <- data.frame(Cov    = v,
                                   Value  = cur.sub[[v]],
                                   Group  = cur.sub[["_grp_"]]);
        cur.v$Strata <- paste("Overall");
        cur.d        <- rbind(cur.d, cur.v);
    }
    cur.d$Group <- as.factor(cur.d$Group);

    rst <- ggplot(data = cur.d, aes(x = Value)) +
        geom_density(alpha = 0.2,
                     aes(group = Group,
                         color = Group,
                         fill  = Group,
                         linetype = Group),
                     na.rm = TRUE) +
        scale_y_continuous(breaks = NULL) +
        labs(x = "", y = "") +
        facet_grid(Strata ~ ., scales = facet.scales);
    rst
}


plotRweBalance <- function(data.withps, overall.inc = TRUE, v.cov = NULL,
                           facet.scales = "free_y", label.cov = v.cov, legend.width = 0.08,
                           ...) {

    if (is.null(v.cov))
        v.cov <- all.vars(data.withps$ps.fml)[-1];

    if (is.null(label.cov))
        label.cov <- v.cov;

    nstrata      <- data.withps$nstrata;
    dtaps        <- data.withps$data;
    dtaps$Strata <- dtaps[["_strata_"]];
    dtaps$Group  <- dtaps[["_grp_"]];

    rst <- list();
    for (v in v.cov) {
        if (is.factor(dtaps[[v]])) {
            cur.p <- plot.balance.fac(dtaps, v, overall.inc = overall.inc);
        } else {
            cur.p <- plot.balance.cont(dtaps, v, nstrata = nstrata,
                                       overall.inc = overall.inc, facet.scales = facet.scales);
        }
        cur.p <- cur.p +
            labs(title = label.cov[v == v.cov]) +
            theme_bw() +
            theme(strip.background = element_blank(),
                  strip.placement  = "right",
                  strip.text       = element_blank(),
                  panel.grid       = element_blank(),
                  panel.border     = element_blank(),
                  panel.spacing    = unit(0, "lines"),
                  plot.title = element_text(hjust = 0.5),
                  legend.position  = "none",
                  plot.margin      = unit(c(1,0,1,-0.5), "lines"));

        rst[[v]] <- cur.p;
    }

    rst[[length(rst)]]     <- rst[[length(rst)]] + theme(strip.text = element_text(size=8),
                                                         legend.position = "right");
    rst$nrow               <- 1;
    rst$rel_widths         <- c(rep(1, length(v.cov)-1), 1+legend.width*length(v.cov));
    do.call(plot_grid, rst);
}
