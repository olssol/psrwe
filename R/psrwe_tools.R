#' @title RWE related class names
#'
#' @noRd
get_rwe_class <- function(c.str) {
  switch(c.str,
         DWITHPS   = "RWE_PS_DTA",
         DPSMATCH  = "RWE_PS_DTA_MAT",
         PSDIST    = "RWE_PS_BOR",
         PSRST     = "RWE_PS_RST",
         CLRST     = "RWE_CL_RST",
         PPRST     = "RWE_POWERPRST",
         ANARST    = "RWE_PS_RST",
         ANARSTPLT = "RWE_PS_RST_PLOT",
         ovl       = "overlapping area",
         ksd       = "Kullback-Leibler distance",
         std       = "standardized difference in means",
         abd       = "absolute difference in means",
         ley       = "Levy distance",
         mhb       = "Mahalanobis distance",
         omkss     = "one minus Kolmogorov-Smirnov statistic",
         n_current = "the number of current control subjects",
         distance  = "the distance in PS distributions",
         inverse_distance = "one minus the inverse of the distance in PS distributions"
         )
}


#' @title Compute propensity scores
#'
#' @noRd
get_ps <- function(dta,
                   ps_method = c("logistic", "randomforest"),
                   ps_fml,
                   ntree = 5000,
                   grp = NULL,
                   ps_cov = NULL,
                   ...) {

  type <- match.arg(ps_method)

  # Generate formula
  if (is.null(ps_fml))
    ps_fml <- as.formula(paste(grp, "~",
                               paste(ps_cov, collapse = "+"),
                               sep = ""))

  # Identify `grp` if passed from formula
  grp <- all.vars(ps_fml)[1]

  # Fit model
  switch(type,
         logistic = {
           glm_fit <- glm(ps_fml, family = binomial, data = dta, ...)
           est_ps  <- glm_fit$fitted
         },
         randomforest = {
           dta[[grp]] <- as.factor(dta[[grp]])
           rf_fit     <- randomForest(ps_fml, data = dta,
                                      ntree = ntree, ...)
           est_ps     <- predict(rf_fit, type = "prob")[, 2]
         })

  return(est_ps)
}


#' @title Get number of subjects borrowed
#'
#' @param total_borrow integer. Target number of subjects to be borrowed.
#' @param ns0 vector. Number of subjects in historical data (control) arm for
#'   each stratum.
#' @param rs vector. Similarity measure; for example, overlapping coefficient
#'   for each stratum.
#' @param m.lambda character. Method to split \code{total_borrow}, which can be
#'   based on distance (\code{m.lambda = "dist"}) or inverse distance
#'   (\code{m.lambda = "inverse"}).
#'
#' @return Data frame with proportion for borrowing
#'
#' @noRd
#'
get_aborrow <- function(total_borrow, ns0, ns1, rs,
                        m_lambda = c("distance",
                                     "inverse_distance",
                                     "n_current")) {

    m_lambda   <- match.arg(m_lambda)
    proportion <- switch(m_lambda,
                         distance = {
                             proportion <- rs / sum(rs)
                         },
                         inverse_distance = {
                             mrs        <- 1 / (1 - rs)
                             proportion <- mrs / sum(mrs)
                         },
                         n_current = {
                             proportion <- ns1 / sum(ns1)
                         })

    borrow     <- apply(cbind(ns0, total_borrow * proportion),
                        1, min)

    alpha      <- borrow / ns0

    data.frame(Proportion = proportion,
               N_Borrow   = borrow,
               Alpha      = alpha)
}


#' @title Get data from current stratum
#'
#' @noRd
#'
get_cur_d <- function(data, i, v_covs, grp_cur = 1, arm_ctl = 0) {
    f_d <- function(x) {
        if (0 == nrow(x)) {
            rst <- NULL
        } else {
            rst <- x[, v_covs, drop = TRUE]
        }

        rst
    }

    cur_d1 <- data[data[["_strata_"]] == i       &
                   data[["_grp_"]]    == grp_cur &
                   data[["_arm_"]]    == arm_ctl,
                   v_covs,
                   drop = FALSE]

    cur_d0 <- data[data[["_strata_"]] == i       &
                   data[["_grp_"]]    != grp_cur &
                   data[["_arm_"]]    == arm_ctl,
                   v_covs,
                   drop = FALSE]

    cur_d1t <- data[data[["_strata_"]] == i       &
                    data[["_grp_"]]    == grp_cur &
                    data[["_arm_"]]    != arm_ctl,
                    v_covs,
                    drop = FALSE]

    ## check data
    if (0 == nrow(cur_d1) | 0 == nrow(cur_d0)) {
        warning(paste(i, "contains no subjects from group 1 or 0"))
    }

    out <- list(cur_d1  = f_d(cur_d1),
                cur_d0  = f_d(cur_d0),
                cur_d1t = f_d(cur_d1t))

    return(out)
}

#' @title Get observed outcome
#'
#' @noRd
#'
get_observed <- function(data, v_covs) {
    data <- data %>%
        mutate(Stratum = `_strata_`,
               Group  = `_grp_`,
               Arm    = `_arm_`,
               Y      = (!!as.name(v_covs))) %>%
        filter(!is.na(Stratum))

    rst1 <- data %>%
        group_by(Group, Arm, Stratum) %>%
        summarize(N    = n(),
                  Mean = mean(Y),
                  SD   = sd(Y))

    rst2 <- data %>%
        mutate(Stratum = "Overall") %>%
        group_by(Group, Arm, Stratum) %>%
        summarize(N    = n(),
                  Mean = mean(Y),
                  SD   = sd(Y))

    rbind(rst1, rst2)
}

#' @title Generate frequency table for factor columns
#'
#' @return A vector with the number of samples in group 0, the number of samples
#'   in group 1, and the distance measure from group 0 to group 1.
#'
#' @noRd
get_freq_tbl <- function(data, var.groupby, vars = NULL) {

  if (is.null(vars))
    vars <- colnames(data)

  rst <- NULL
  for (v in vars) {
    if (!is.factor(data[[v]]))
      next

    cur_freq <- data %>%
        count(.dots = c(var.groupby, v)) %>%
        group_by(!!as.name(var.groupby)) %>%
        mutate(Sum  = sum(.data$n),
               Freq = .data$n / sum(.data$n)) %>%
        ungroup() %>%
        mutate_if(is.factor, as.character) %>%
        mutate(Cov = v) %>%
        rename(Value = v)

    rst <- rbind(rst,
                 data.frame(cur_freq))
  }

  return(rst)
}


#' @title Plot density distribution of propensity score
#'
#' @noRd
plot_ps <- function(data.withps,
                    overall.inc = TRUE,
                    add.text = TRUE,
                    facet.scales = "free_y",
                    ...) {

    N0 <- N1 <- Dist <- Ps <- Group <- NULL

    ## stopifnot(inherits(data.withps,
    ##                    what = get_rwe_class("DWITHPS")))

    rst_sum <- summary(data.withps, ...)
    nstrata <- data.withps$nstrata
    dtaps   <- data.withps$data %>%
        filter(!is.na(`_strata_`))

    xlim   <- range(dtaps[which(!is.na(dtaps[["_strata_"]])), "_ps_"],
                    na.rm = TRUE)
    strata <- levels(dtaps[["_strata_"]])

    all.data <- NULL
    for (i in strata) {
        cur.sub  <- dtaps[which(i == dtaps[["_strata_"]]), ]
        cur.data <- data.frame(Strata = i,
                               Ps     = cur.sub[["_ps_"]],
                               Group  = cur.sub[["_grp_"]])
        all.data <- rbind(all.data, cur.data)
    }
    ps_kl <- data.frame(Strata = rst_sum$Summary$Stratum,
                        N0     = rst_sum$Summary$N_RWD,
                        N1     = rst_sum$Summary$N_Current)

    if (overall.inc) {
        cur.data <-  data.frame(Strata = "Overall",
                                Ps     = dtaps[["_ps_"]],
                                Group  = dtaps[["_grp_"]])

        all.data <- rbind(all.data, cur.data)
        ps_kl    <- rbind(ps_kl,
                          data.frame(Strata = "Overall",
                                     N0     = rst_sum$Overall$N_RWD,
                                     N1     = rst_sum$Overall$N_Current))
    }

    all.data$Group <- as.factor(all.data$Group)

    rst <- ggplot(data = all.data, aes(x = Ps)) +
        geom_density(alpha = 0.2,
                     aes(group = Group,
                         fill  = Group,
                         linetype = Group),
                     trim  = TRUE,
                     na.rm = TRUE) +
        labs(x = "Propensity Score", y = "Density") +
        scale_y_continuous(breaks = NULL) +
        scale_x_continuous(limits = xlim) +
        scale_fill_manual(values = c("gray20", "gray80")) +
        theme_bw() +
        theme(strip.background = element_blank(),
              panel.grid = element_blank(),
              panel.border = element_rect(colour = "black"),
              panel.spacing = unit(0, "lines")) +
        facet_grid(Strata ~ ., scales = facet.scales)

    if (add.text) {
        rst <- rst +
            geom_text(x = Inf, y = Inf, hjust = 1, vjust = 1,
                      aes(label = paste('N0 =', N0, ", N1 =", N1,
                                        sep = "")),
                      data = ps_kl, size = 4)
    }
    return(rst)
}


#' @title Plot balance distributions for factors
#'
#' @noRd
plot_balance_fac <- function(dtaps, v,
                             overall.inc = TRUE) {

    cur.d <- get_freq_tbl(dtaps,
                          var.groupby = c("Strata", "Group"),
                          vars = v)
    cur.d <- cur.d %>%
        dplyr::filter(!is.na(Strata))

    if (overall.inc) {
        cur.overall <- get_freq_tbl(dtaps,
                                    var.groupby = "Group",
                                    vars = v)
        cur.overall$Strata <- "Overall"
        cur.d <- rbind(cur.d, cur.overall)
    }

    cur.d$Group <- as.factor(cur.d$Group)
    cur.d$Value <- as.factor(cur.d$Value)

    rst <- ggplot(data = cur.d, aes(x = .data$Value, y = .data$Freq)) +
        geom_bar(alpha = 0.4,
                 stat = "identity",
                 position = "dodge",
                 color = "black",
                 aes(group = .data$Group,
                     fill  = .data$Group)) +
        scale_fill_manual(values = c("gray20", "gray80")) +
        scale_y_continuous(breaks = NULL, limits = c(0, 1)) +
        labs(x = "", y = "") +
        facet_grid(Strata ~ .)

    return(rst)
}


#' @title Plot balance distributions for factors
#'
#' @noRd
plot_balance_cont <- function(dtaps, v, strata,
                              overall.inc = TRUE,
                              facet.scales = "free_y") {

  Value <- Group <- NULL
  cur.d <- NULL
  for (i in strata) {
    cur.sub      <- dtaps[which(i == dtaps[["_strata_"]]), ]
    cur.v        <- data.frame(Cov    = v,
                               Value  = cur.sub[[v]],
                               Group  = cur.sub[["_grp_"]])
    cur.v$Strata <- i
    cur.d        <- rbind(cur.d, cur.v)
  }

  if (overall.inc) {
    cur.sub      <- dtaps
    cur.v        <- data.frame(Cov    = v,
                               Value  = cur.sub[[v]],
                               Group  = cur.sub[["_grp_"]])
    cur.v$Strata <- paste("Overall")
    cur.d        <- rbind(cur.d, cur.v)
  }
  cur.d$Group <- as.factor(cur.d$Group)

  rst <- ggplot(data = cur.d, aes(x = Value)) +
    geom_density(alpha = 0.2,
                 aes(group = Group,
                     fill  = Group,
                     linetype = Group),
                 na.rm = TRUE) +
    scale_y_continuous(breaks = NULL) +
    scale_fill_manual(values = c("gray20", "white")) +
    labs(x = "", y = "") +
    facet_grid(Strata ~ ., scales = facet.scales)

  return(rst)
}


#' @title Plot the balance of baseline variables
#'
#' @noRd
plot_balance <- function(data.withps,
                         overall.inc = TRUE,
                         v.cov = NULL,
                         facet.scales = "free_y",
                         label.cov = v.cov,
                         legend.width = 0.08,
                         ...) {

    if (is.null(v.cov)) {
        v.cov <- all.vars(data.withps$ps_fml)[-1]
    }

    ## remove stratification covs
    if (!is.null(data.withps$strata_covs)) {
        s_cov_inx <- which(v.cov %in% data.withps$strata_covs)
        v.cov     <- v.cov[-s_cov_inx]
    }

    if (is.null(label.cov)) {
        label.cov <- v.cov
    }


    nstrata      <- data.withps$nstrata
    dtaps        <- data.withps$data
    dtaps$Strata <- dtaps[["_strata_"]]
    dtaps$Group  <- dtaps[["_grp_"]]
    strata       <- levels(dtaps[["_strata_"]])

    rst <- list()
    for (v in v.cov) {
        if (is.factor(dtaps[[v]])) {
            cur.p <- plot_balance_fac(dtaps, v, overall.inc = overall.inc)
        } else {
            cur.p <- plot_balance_cont(dtaps, v, strata = strata,
                                       overall.inc = overall.inc,
                                       facet.scales = facet.scales)
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
                  plot.title       = element_text(hjust = 0.5),
                  legend.position  = "none",
                  plot.margin      = unit(c(1, 0, 1, -0.5), "lines"))

        rst[[v]] <- cur.p
    }

    rst[[length(rst)]] <- rst[[length(rst)]] +
        theme(strip.text = element_text(size = 8),
              legend.position = "right")

    rst$nrow       <- 1
    rst$rel_widths <- c(rep(1, length(v.cov) - 1),
                        1 + legend.width * length(v.cov))
    do.call(plot_grid, rst)
}

#' @title Get strata by covariates
#'
#' @noRd
#'
get_strata <- function(data, strata_covs  = NULL) {

    if (is.null(strata_covs)) {
        data[["_strata_"]] <- "Stratum All"
        nstrata            <- 1
    } else {
        strata <- data[, strata_covs, drop = F] %>%
            distinct()
        nstrata <- nrow(strata)

        s_label <- NULL
        for (i in 1:nstrata) {
            cur_label <- paste(strata_covs, "=", sep = "")
            cur_label <- paste(cur_label, strata[i, ], sep = "")
            cur_label <- paste(cur_label, collapse = ",")
            s_label   <- c(s_label, cur_label)
        }
        strata$label  <- factor(s_label)

        suppressMessages(
            data <- data %>%
                left_join(strata) %>%
                mutate(`_strata_` = label) %>%
                select(-label)
        )
    }

    return(list(data    = data,
                nstrata = nstrata))

}

#' @title Print information
#'
#'
#' @noRd
#'
cat_ps_dta <- function(x, rst_sum) {

    ss <- paste("This is a",
                ifelse(x$is_rct, "randomized", "sing-arm"),
                "study. A total of", as.character(rst_sum$N["RWD"]),
                "RWD subjects and", rst_sum$N["Current"],
                "current study subjects are used to",
                "estimate propensity",
                "scores by", x$ps_method, "model.",
                "A total of", rst_sum$N["Trimmed"],
                "RWD subjects are trimmed",
                "and excluded from the final analysis.",
                "The following covariates are adjusted in the propensity",
                "score model:",
                paste(all.vars(x$ps_fml[-1]), collapse = ", "),
                sep = " ")
    ss <- paste(ss, ".", sep = "")

    cat_paste(ss)
}

#' Break string into words and cat by fixed width
#'
#' @noRd
cat_paste <- function(ss, fill = 60) {
    do.call(cat, c(strsplit(ss, " "), fill = fill))
}

#' Jackknife variance
#'
#'
#'
#' @noRd
#'
get_jk_sd <- function(overall_mean, jk_all) {
    var_theta <- (length(jk_all) - 1) / length(jk_all)
    var_theta <- var_theta * sum((jk_all - overall_mean)^2)
    sqrt(var_theta)
}

#' Get estimates for composite likelihood and survival
#'
#'
#'
#' @noRd
#'
get_ps_cl_km <- function(dta_psbor,
                         v_outcome  = NULL,
                         v_event    = NULL,
                         v_time     = NULL,
                         f_stratum  = get_cl_stratum,
                         ...) {

    ## prepare data
    is_rct  <- dta_psbor$is_rct
    data    <- dta_psbor$data
    data    <- data[!is.na(data[["_strata_"]]), ]

    strata  <- levels(data[["_strata_"]])
    nstrata <- length(strata)
    borrow  <- dta_psbor$Borrow$N_Borrow

    ## estimate
    ctl_theta <- NULL
    trt_theta <- NULL

    for (i in seq_len(nstrata)) {
        cur_01  <- get_cur_d(data,
                             strata[i],
                             c(v_outcome, v_time, v_event))

        cur_d1  <- cur_01$cur_d1
        cur_d0  <- cur_01$cur_d0
        cur_d1t <- cur_01$cur_d1t

        ## control with borrowing
        cur_ctl   <- f_stratum(cur_d1, cur_d0, n_borrow = borrow[i], ...)
        ctl_theta <- rbind(ctl_theta, cur_ctl)
        if (is_rct) {
            cur_trt   <- f_stratum(cur_d1t, ...)
            trt_theta <- rbind(trt_theta, cur_trt)
        }
    }

    ## summary
    rst_trt    <- NULL
    rst_effect <- NULL
    if (is_rct) {
        rst_trt    <- get_overall_est(trt_theta[, 1],
                                      trt_theta[, 2],
                                      dta_psbor$Borrow$N_Cur_TRT)
        rst_effect <- get_overall_est(trt_theta[, 1] - ctl_theta[, 1],
                                      sqrt(trt_theta[, 2] + ctl_theta[, 2]),
                                      dta_psbor$Borrow$N_Current)
        n_ctl      <- dta_psbor$Borrow$N_Cur_CTL
    } else {
        n_ctl      <- dta_psbor$Borrow$N_Current
    }
    rst_ctl <- get_overall_est(ctl_theta[, 1], ctl_theta[, 2], n_ctl)

    ## return
    rst <-  list(Control   = rst_ctl,
                 Treatment = rst_trt,
                 Effect    = rst_effect,
                 Borrow    = dta_psbor$Borrow,
                 Total_borrow = dta_psbor$Total_borrow,
                 is_rct       = is_rct)
}

#' Summarize overall theta
#'
#'
#' @noRd
#'
get_overall_est <- function(theta, sds, weights) {
    ws         <- weights / sum(weights)
    overall    <- sum(theta * ws)
    sd_overall <- sqrt(sum(ws^2 * sds^2))

    list(Stratum_Estimate = cbind(Mean = theta,
                                  SD   = sds),
         Overall_Estimate = c(Mean = overall,
                              SD   = sd_overall))
}
