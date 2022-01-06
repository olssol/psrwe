#' @title RWE related class names
#'
#' @noRd
get_rwe_class <- function(c_str) {
  switch(c_str,
         DWITHPS   = "PSRWE_DTA",
         DPSMATCH  = "PSRWE_DTA_MAT",
         PSDIST    = "PSRWE_BOR",
         PSRST     = "PSRWE_RST",
         CLRST     = "RWE_CL_RST",
         PPRST     = "RWE_POWERPRST",
         ANARST    = "PSRWE_RST",
         OUTANA    = "PSRWE_RST_OUTANA",
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
#' @param total_borrow integer. Target number of subjects to be borroqwed.
#' @param ns0 vector. Number of subjects in historical data (control) arm for
#'   each stratum.
#' @param rs vector. Similarity measure; for example, overlapping coefficient
#'   for each stratum.
#' @param m_lambda character. Method to split \code{total_borrow}, which can be
#'   based on distance (\code{m_lambda = "dist"}) or inverse distance
#'   (\code{m_lambda = "inverse"}).
#'
#' @return Data frame with proportion for borrowing
#'
#' @noRd
#'
get_aborrow <- function(total_borrow, ns0, ns1, rs,
                        m_lambda = c("distance",
                                     "inverse_distance",
                                     "n_current"),
                        ...) {

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
        summarize(N      = n(),
                  Mean   = mean(Y),
                  StdErr = sd(Y))

    rst2 <- data %>%
        mutate(Stratum = "Overall") %>%
        group_by(Group, Arm, Stratum) %>%
        summarize(N      = n(),
                  Mean   = mean(Y),
                  StdErr = sd(Y))

    rst <- data.frame(rbind(rst1, rst2))
    rst
}

#' @title Generate frequency table for factor columns
#'
#' @return A vector with the number of samples in group 0, the number of samples
#'   in group 1, and the distance measure from group 0 to group 1.
#'
#' @noRd
get_freq_tbl <- function(data, var_groupby, vars = NULL) {

  if (is.null(vars))
    vars <- colnames(data)

  rst <- NULL
  for (v in vars) {
    if (!is.factor(data[[v]]))
      next

    cur_freq <- data %>%
        count(.dots = c(var_groupby, v)) %>%
        group_by(!!as.name(var_groupby)) %>%
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
plot_ps <- function(data_withps,
                    overall_inc = TRUE,
                    add_text = TRUE,
                    facet_scales = "free_y",
                    ...) {

    N0 <- N1 <- Dist <- Ps <- Group <- NULL

    ## stopifnot(inherits(data_withps,
    ##                    what = get_rwe_class("DWITHPS")))

    rst_sum <- summary(data_withps, ...)
    nstrata <- data_withps$nstrata
    dtaps   <- data_withps$data %>%
        filter(!is.na(`_strata_`))

    xlim   <- range(dtaps[which(!is.na(dtaps[["_strata_"]])), "_ps_"],
                    na.rm = TRUE)
    strata <- levels(dtaps[["_strata_"]])

    all_data <- NULL
    for (i in strata) {
        cur_sub  <- dtaps[which(i == dtaps[["_strata_"]]), ]
        cur_data <- data.frame(Strata = i,
                               Ps     = cur_sub[["_ps_"]],
                               Group  = cur_sub[["_grp_"]])
        all_data <- rbind(all_data, cur_data)
    }
    ps_kl <- data.frame(Strata = rst_sum$Summary$Stratum,
                        N0     = rst_sum$Summary$N_RWD,
                        N1     = rst_sum$Summary$N_Current)

    if (overall_inc) {
        cur_data <-  data.frame(Strata = "Overall",
                                Ps     = dtaps[["_ps_"]],
                                Group  = dtaps[["_grp_"]])

        all_data <- rbind(all_data, cur_data)
        ps_kl    <- rbind(ps_kl,
                          data.frame(Strata = "Overall",
                                     N0     = rst_sum$Overall$N_RWD,
                                     N1     = rst_sum$Overall$N_Current))
    }

    all_data$Group <- as.factor(all_data$Group)

    rst <- ggplot(data = all_data, aes(x = Ps)) +
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
        facet_grid(Strata ~ ., scales = facet_scales)

    if (add_text) {
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
                             overall_inc = TRUE) {

    cur_d <- get_freq_tbl(dtaps,
                          var_groupby = c("Strata", "Group"),
                          vars = v)
    cur_d <- cur_d %>%
        dplyr::filter(!is.na(Strata))

    if (overall_inc) {
        cur_overall <- get_freq_tbl(dtaps,
                                    var_groupby = "Group",
                                    vars = v)
        cur_overall$Strata <- "Overall"
        cur_d <- rbind(cur_d, cur_overall)
    }

    cur_d$Group <- as.factor(cur_d$Group)
    cur_d$Value <- as.factor(cur_d$Value)

    rst <- ggplot(data = cur_d, aes(x = .data$Value, y = .data$Freq)) +
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
                              overall_inc = TRUE,
                              facet_scales = "free_y") {

  Value <- Group <- NULL
  cur_d <- NULL
  for (i in strata) {
    cur_sub      <- dtaps[which(i == dtaps[["_strata_"]]), ]
    cur_v        <- data.frame(Cov    = v,
                               Value  = cur_sub[[v]],
                               Group  = cur_sub[["_grp_"]])
    cur_v$Strata <- i
    cur_d        <- rbind(cur_d, cur_v)
  }

  if (overall_inc) {
    cur_sub      <- dtaps
    cur_v        <- data.frame(Cov    = v,
                               Value  = cur_sub[[v]],
                               Group  = cur_sub[["_grp_"]])
    cur_v$Strata <- paste("Overall")
    cur_d        <- rbind(cur_d, cur_v)
  }
  cur_d$Group <- as.factor(cur_d$Group)

  rst <- ggplot(data = cur_d, aes(x = Value)) +
    geom_density(alpha = 0.2,
                 aes(group = Group,
                     fill  = Group,
                     linetype = Group),
                 na.rm = TRUE) +
    scale_y_continuous(breaks = NULL) +
    scale_fill_manual(values = c("gray20", "white")) +
    labs(x = "", y = "") +
    facet_grid(Strata ~ ., scales = facet_scales)

  return(rst)
}


#' @title Plot the balance of baseline variables
#'
#' @noRd
plot_balance <- function(data_withps,
                         overall_inc = TRUE,
                         v_cov = NULL,
                         facet_scales = "free_y",
                         label_cov = v_cov,
                         legend_width = 0.08,
                         ...) {

    if (is.null(v_cov)) {
        v_cov <- all.vars(data_withps$ps_fml)[-1]
    }

    ## remove stratification covs
    if (!is.null(data_withps$strata_covs)) {
        s_cov_inx <- which(v_cov %in% data_withps$strata_covs)
        v_cov     <- v_cov[-s_cov_inx]
    }

    if (is.null(label_cov)) {
        label_cov <- v_cov
    }

    nstrata      <- data_withps$nstrata
    dtaps        <- data_withps$data
    dtaps$Strata <- dtaps[["_strata_"]]
    dtaps$Group  <- dtaps[["_grp_"]]
    strata       <- levels(dtaps[["_strata_"]])

    rst <- list()
    for (v in v_cov) {
        if (is.factor(dtaps[[v]])) {
            cur_p <- plot_balance_fac(dtaps, v, overall_inc = overall_inc)
        } else {
            cur_p <- plot_balance_cont(dtaps, v, strata = strata,
                                       overall_inc = overall_inc,
                                       facet_scales = facet_scales)
        }

        cur_p <- cur_p +
            labs(title = label_cov[v == v_cov]) +
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

        rst[[v]] <- cur_p
    }

    rst[[length(rst)]] <- rst[[length(rst)]] +
        theme(strip.text = element_text(size = 8),
              legend.position = "right")

    rst$nrow       <- 1
    rst$rel_widths <- c(rep(1, length(v_cov) - 1),
                        1 + legend_width * length(v_cov))
    do.call(plot_grid, rst)
}


#' @title Plot the standardized (absolute) mean differences of baseline
#'     variables
#'
#' @noRd
plot_astd <- function(data_withps,
                      metric = c("std", "astd"),
                      avg_only = FALSE,
                      ...) {

    v_cov <- all.vars(data_withps$ps_fml)[-1]

    ## check arguments
    d_metric <- match.arg(metric)
    if (d_metric == "astd") {
        xlab       <- "Standardized Absolute Mean Difference"
        xintercept <- c(0.2, 0.4)
    } else {
        xlab       <- "Standardized Mean Difference"
        xintercept <- c(-0.4, -0.2, 0.2, 0.4)
    }

    ## remove stratification covs
    if (!is.null(data_withps$strata_covs)) {
        s_cov_inx <- which(v_cov %in% data_withps$strata_covs)
        v_cov     <- v_cov[-s_cov_inx]
    }

    ## prepare data
    dtaps        <- data_withps$data
    dtaps$Strata <- dtaps[["_strata_"]]
    dtaps$Group  <- dtaps[["_grp_"]]
    strata       <- levels(dtaps[["_strata_"]])

    dta_asd <- data.frame()
    for (v in v_cov) {
        ## original data without any trimming
        cov0    <- as.numeric(dtaps[[v]][dtaps$Group == 0])
        cov1    <- as.numeric(dtaps[[v]][dtaps$Group == 1])
        std_all <- get_distance(cov0, cov1, metric = d_metric)
        dta_asd <- rbind(dta_asd,
                         data.frame(v_cov = v,
                                    Group = "Before Stratification",
                                    asd   = std_all))

        ## PS stratified data with trimming
        std_ws <- NULL
        for (s in strata) {
            cov0 <- as.numeric(dtaps[[v]][dtaps$Group == 0 &
                                          dtaps$Strata == s &
                                          !is.na(dtaps$Strata)])
            cov1 <- as.numeric(dtaps[[v]][dtaps$Group == 1 &
                                          dtaps$Strata == s &
                                          !is.na(dtaps$Strata)])
            std_s <- get_distance(cov0, cov1, metric = d_metric)

            if (!avg_only) {
                dta_asd <- rbind(dta_asd,
                                 data.frame(v_cov = v,
                                            Group = s,
                                            asd   = std_s))
	    }

            std_ws <- c(std_ws, std_s)
        }

        if (avg_only) {
            dta_asd <- rbind(dta_asd,
                             data.frame(v_cov = v,
                                        Group = "After Stratification",
                                        asd   = mean(std_ws)))

            dta_asd$Group <- factor(dta_asd$Group,
                                    levels = c("Before Stratification",
                                               "After Stratification"))
        }
    }

    ## plot
    rst <- ggplot(dta_asd,
                  aes(x = asd, y = v_cov, shape = Group, color = Group)) +
           geom_point() +
           geom_vline(xintercept = 0, linetype = 2) +
           geom_vline(xintercept = xintercept, linetype = 3) +
           labs(x = xlab, y = "Variables") +
           theme_bw() +
           theme(strip.background = element_blank(),
                 panel.border = element_rect(colour = "black"),
                 panel.spacing = unit(0, "lines"))

    return(rst)
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
                         v_outcome     = NULL,
                         v_event       = NULL,
                         v_time        = NULL,
                         f_stratum     = get_cl_stratum,
                         f_overall_est = get_overall_est,
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
        rst_trt    <- f_overall_est(trt_theta, dta_psbor$Borrow$N_Cur_TRT)
        rst_effect <- f_overall_est(trt_theta, dta_psbor$Borrow$N_Current,
                                    ctl_theta)
        n_ctl      <- dta_psbor$Borrow$N_Cur_CTL
    } else {
        n_ctl      <- dta_psbor$Borrow$N_Current
    }
    rst_ctl <- f_overall_est(ctl_theta, n_ctl)

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
get_overall_est <- function(ts1, weights, ts2 = NULL) {

    if (is.null(ts2)) {
        theta0 <- ts1[, 1]
        sds0   <- ts1[, 2]
    } else {
        theta0 <- ts1[, 1] - ts2[, 1]
        sds0   <- sqrt(ts1[, 2] + ts2[, 2])
    }

    ws         <- weights / sum(weights)
    nstrata    <- length(ws)
    theta      <- matrix(theta0, ncol = nstrata)
    sds        <- matrix(sds0,   ncol = nstrata)

    overall    <- as.vector(theta %*% ws)
    sd_overall <- as.vector(sqrt(sds^2 %*% ws^2))

    ## stratum est
    s_est <- data.frame(Mean   = theta0,
                        StdErr = sds0)
    o_est <- data.frame(Mean   = overall,
                        StdErr = sd_overall)

    if (ncol(ts1) > 2) {
        prept <- matrix(ts1[, 3], ncol = nstrata)
        s_est <- cbind(s_est, T = prept[, 1],
                       Stratum = rep(1:nstrata, each = nrow(theta)))
        o_est <- cbind(o_est, T = prept[, 1])
    }

    list(Stratum_Estimate = s_est,
         Overall_Estimate = o_est)
}


#'  Plot density for power prior results
#'
#'
#' @noRd
#'
plot_pp_rst <- function(x,
                        add_stratum = FALSE,
                        split_rct_arm = FALSE,
                        ...) {
    ## prepare data
    if (x$is_rct){
      label_Arm <- "Control"
    } else {
      label_Arm <- "Single"
    }

    label_Type <- c("Arm Specific", "Arm Specific", "Treatment Effect")
    if (x$is_rct && split_rct_arm) {
      label_Type <- c("Arm Control", "Arm Treatment", "Treatment Effect")
    }

    rst <- data.frame(Type    = label_Type[1],
                      Arm     = label_Arm,
                      Stratum = "Overall",
                      theta   = x$Control$Overall_Samples)

    if (x$is_rct) {
        rst <- rbind(rst,
                     data.frame(Type    = label_Type[2],
                                Arm     = "Treatment",
                                Stratum = "Overall",
                                theta   = x$Treatment$Overall_Samples),
                     data.frame(Type    = label_Type[3],
                                Arm     = "Effect",
                                Stratum = "Overall",
                                theta   = x$Effect$Overall_Samples))
    }

    ## add stratum
    if (add_stratum) {
        rst <- rbind(rst,
                     data.frame(Type    = label_Type[1],
                                Arm     = label_Arm,
                                Stratum = x$Borrow$Stratum,
                                theta   = as.vector(x$Control$Stratum_Samples)))

        if (x$is_rct) {
            rst <- rbind(rst,
                         data.frame(Type    = label_Type[2],
                                    Arm     = "Treatment",
                                    Stratum = x$Borrow$Stratum,
                                    theta   = as.vector(x$Treatment$Stratum_Samples)),
                         data.frame(Type    = label_Type[3],
                                    Arm     = "Effect",
                                    Stratum = x$Borrow$Stratum,
                                    theta   = as.vector(x$Effect$Stratum_Samples))
	                 )
        }
    }

    ## plot
    rst$Arm_by_Stratum <- interaction(rst$Arm, rst$Stratum)
    if (x$is_rct) {
        rst_plt <- ggplot(data = rst, aes(x = theta)) +
            stat_density(aes(group = Arm_by_Stratum,
			     color = Arm,
                             linetype = Stratum),
                         position  = "identity",
                         geom      = "line",
                         adjust    = 1.2,
                         trim      = TRUE) +
            facet_wrap(~ Type, scales = "free")
    } else {
        rst_plt <- ggplot(data = rst, aes(x = theta)) +
            stat_density(aes(group    = Arm_by_Stratum,
                             color = Stratum,
                             linetype = Stratum),
                         position  = "identity",
                         geom      = "line",
                         adjust    = 1.2)
    }
    rst_plt <- rst_plt +
        theme_bw() +
        labs(x = expression(theta), y = "Density")

    rst_plt
}

#' @title Plot KM at all time points
#'
#'
#' @noRd
#'
plot_km_rst <- function(x,
                        xlab = "Time",
                        ylab = "Survival Probability",
                        add_ci = TRUE,
                        add_stratum = FALSE,
                        ...) {

    ## prepare data
    if (x$is_rct){
      label_Arm <- "Control"
    } else {
      label_Arm <- "Single"
    }

    rst <- data.frame(Arm = paste(label_Arm, "Overall", sep = " "),
                      x$Control$Overall_Estimate)

    if (x$is_rct) {
        rst <- rbind(rst,
                     cbind(Arm = "Treatment Overall",
                           x$Treatment$Overall_Estimate))
    }

    ## add stratum
    if (add_stratum) {
        name_v <- c("Mean", "StdErr", "T")
        name_s <- x$Borrow$Stratum[x$Control$Stratum_Estimate$Stratum]
        rst <- rbind(rst,
                     cbind(Arm = paste(label_Arm, name_s, sep = " "),
                           x$Control$Stratum_Estimate[, name_v]))

        if (x$is_rct) {
            name_s <- x$Borrow$Stratum[x$Treatment$Stratum_Estimate$Stratum]
            rst <- rbind(rst,
                         cbind(Arm = paste("Treatment", name_s, sep = " "),
                               x$Treatment$Stratum_Estimate[, name_v]))
        }
    }

    ## CI
    if (add_ci) {
      ci  <- get_ci_km(rst$Mean, rst$StdErr, ...)
      rst <- cbind(rst, Lower = ci$Lower, Upper = ci$Upper)
    }

    ## check arguments
    args <- list(...)
    if ("xlim" %in% names(args)) {
        xlim <- args[['xlim']]
    } else {
        xlim <- range(rst$T)
    }

    if ("ylim" %in% names(args)) {
        ylim <- args[['ylim']]
    } else {
        ylim <- c(min(0, rst$Mean), max(1, rst$Mean))

        if (add_ci) {
            ylim <- c(min(0, rst$Lower), max(1, rst$Upper))
        }
    }

    ## plot
    # lt_a <- rep(2, length(rst$Arm))
    # lt_a[grep(".* Overall$", rst$Arm)] <- 1
    rst_plt <- ggplot(data = rst) +
        geom_step(aes(x = T, y = Mean, col = Arm, linetype = Arm)) +
        scale_y_continuous(limits = ylim) +
        scale_x_continuous(limits = xlim) +
        labs(x = xlab, y = ylab) +
        theme_bw()

    if (add_ci) {
      rst_plt <- rst_plt +
          geom_step(aes(x = T, y = Lower, col = Arm), linetype = 3) +
          geom_step(aes(x = T, y = Upper, col = Arm), linetype = 3)
    }

    rst_plt
}

## ------------------------------------------------------------
##
##                MATCHING METHODS
##
## ------------------------------------------------------------

#' @title optmatch method
#'
#' @noRd
#'
get_match_optm <- function(data, ratio, caliper, ...) {
    ## prepare data
    dta_sub <- data.frame(gid = data[["_grp_"]],
                          psv = data[["_ps_"]],
                          sid = data[["_strata_"]])

    ## build distance matrix by stratum and within caliper distance
    mat_dm <- optmatch::match_on(gid ~ psv + strata(sid), data = dta_sub,
                                 method = "euclidean")
    mat_dm <- mat_dm + optmatch::caliper(mat_dm, width = caliper)

    ## optmatch
    pm <- optmatch::pairmatch(mat_dm, data = dta_sub, controls = ratio)

    ## match
    id_matched <- !is.na(pm)
    to_match   <- data %>%
        dplyr::filter(1 == `_grp_` & 0 == `_arm_`)

    data[["_matchn_"]]   <- NA
    data[["_matchid_"]]  <- NA

    for (i in seq_len(nrow(to_match))) {
        cur_id    <- to_match[i, "_id_"]
        cur_match <- data[id_matched &
                          pm == pm[cur_id] &
                          data$"_id_" != cur_id, ]

        cur_matchn <- nrow(cur_match)

        ## update
        data[cur_id, "_matchn_"] <- cur_matchn
        if (cur_matchn > 0) {
            cur_matchid <- cur_match[1:cur_matchn, "_id_"]
            data[cur_matchid, "_matchid_"] <- cur_id
        }
    }

    data[which(0 == data[["_grp_"]] &
               is.na(data[["_matchid_"]])),
         "_strata_"] <- NA

    return(data)
}


#' @title Nearest neighbor without replacement matching method by CG
#'
#' @noRd
#'
get_match_nnwor <- function(data, ratio, caliper, ...) {
    ## match
    to_match <- data %>%
        dplyr::filter(1 == `_grp_` & 0 == `_arm_`)

    data[["_matchn_"]]   <- NA
    data[["_matchid_"]]  <- NA

    ## random order nearest neighbor match
    to_match_id <- sample(nrow(to_match))
    for (i in to_match_id) {
        cur_id    <- to_match[i, "_id_"]
        cur_stra  <- to_match[i, "_strata_"]
        cur_ps    <- to_match[i, "_ps_"]

        cur_match <- data %>%
            filter(0        == `_grp_`    &
                   cur_stra == `_strata_` &
                   is.na(`_matchid_`)) %>%
            mutate(dif_ps = abs(`_ps_` - cur_ps)) %>%
            filter(dif_ps <= caliper) %>%
            arrange(dif_ps)

        cur_matchn <- min(nrow(cur_match), ratio)

        ## update
        data[cur_id, "_matchn_"] <- cur_matchn
        if (cur_matchn > 0) {
            cur_matchid <- cur_match[1:cur_matchn, "_id_"]
            data[cur_matchid, "_matchid_"] <- cur_id
        }
    }

    data[which(0 == data[["_grp_"]] &
               is.na(data[["_matchid_"]])),
         "_strata_"] <- NA

    return(data)
}
