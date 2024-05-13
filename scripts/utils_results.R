table_adf <- function(data, start, vars_names, title = "Teste Dickey-Fuller Aumentado",
  table_names = c("Variáveis", "Estatística", "Lag", "P-valor"), ...) {
  result_adf <- data %>%
    filter(date >= start) %>%
    select(-date) %>%
    map_dfr(~ tseries::adf.test(na.omit(.x))[c("statistic", "parameter", "p.value")])
  
  result_adf %>%
    mutate(
      across(-parameter, ~ formatC(.x, digits = 2, format = "f")),
      var = vars_names, .before = 1,
      across(everything(), unname)
    ) %>%
    set_names(table_names) %>%
    stargazer(summary = FALSE, title = title, label = "tb:testadf", ...)
}


ccf_ba <- function(data_list, labels, ...) {
  data_ccf <- imap_dfr(data_list, function(df, name) {
    ci <- qnorm((1 - 0.95) / 2) / sqrt(nrow(df))
    ccf_values <- ccf(df$cases, df$deaths, plot = FALSE, ...)
    
    tibble(
      name = as.factor(name),
      lag = ccf_values$lag[,,1],
      value = ccf_values$acf[,,1],
      ci = ci
    )
  })
  
  ggplot(data_ccf, aes(x = lag, y = value)) +
    geom_vline(xintercept = 0, color = "darkgray") +
    geom_hline(yintercept = 0) +
    geom_segment(aes(xend = lag, yend = 0)) + 
    geom_ribbon(aes(ymin = -ci, ymax = ci), linetype = 2, fill = NA, color = pal[1]) +
    facet_wrap(vars(name), ncol = 3, labeller = labeller(
      name = set_names(labels, names(data_list))
    )) +
    labs(y = "Correlação", x = "Lag (casos)")
}


granger_ba = function(data, p, test, labels = c("Antes", "Depois", "Depois+"), ...){
  pvalue.ast = function(x) {
    c("", "*", "**", "***")[pmap_dbl(tibble(x), ~ sum(. < c(0.01, 0.05, 0.10)) + 1)]
  }
  
  mods <- map(data, ~ VAR(select(.x, -one_of("vaccines")), p))
  mods_vcov <- map(mods, ~ vcovHC(.x, "HC1")) #heteroscedasticity correction
  
  data_granger <- map2(mods, mods_vcov, function(mod, vcov) {
    causality(mod, cause = "cases", vcov. = vcov)[[test]] %>% 
      `[`(c("statistic", "parameter", "p.value")) %>%
      unlist()
  }) %>%
    reduce(rbind) %>%
    as_tibble() %>%
    set_names(if (test == "Granger") {
      c("Estatística", "GL 1", "GL 2", "P-valor")
    } else {
      c("Estatística", "GL", "P-valor")
    })
  
  data_granger %>%
    mutate(
      Estatística = round(`Estatística`, 2),
      across(contains("GL"), round),
      `P-valor` = paste0(round(`P-valor`, 4), pvalue.ast(`P-valor`))
    ) %>%
    mutate(Variável = labels, .before = 1) %>%
    stargazer(
      summary = FALSE, title = "Causalidade de Granger",
      label = "tb:grangerba", notes = "GL: graus de liberdade",
      rownames = FALSE, ...
    )
}


ggvar_select <- function(mod_select, lag.max, remove = "FPE(n)") {
  mod_select %>%
    imap_dfr(function(sel, mod) {
      crit <- as_tibble(apply(t(sel$criteria), 2, \(x) x/x[1]))
      tibble(lag = 1:lag.max, mod = mod, select(crit, -one_of(remove)))
    }) %>%
    pivot_longer(-c(lag, mod)) %>%
    mutate(mod = factor(mod, unique(mod))) %>%
    ggplot(aes(lag, value, color = name)) +
    geom_line() +
    facet_wrap(vars(mod), scales = "free_y", labeller = labeller(mod = namings$ba_mods)) +
    scale_x_continuous(breaks = scales::pretty_breaks())
}


get_yse <- function(object, n.ahead, n_pred, ci = 0.95) {
  Zy <- as.matrix(object$datamat[, 1:(object$K * (object$p + 1))])
  yse <- matrix(NA, nrow = n.ahead, ncol = object$K)
  sig.y <- vars:::.fecov(x = object, n.ahead = n.ahead)
  for (i in 1:n.ahead) {
    yse[i, ] <- sqrt(diag(sig.y[, , i]))
  }
  yse <- -1 * qnorm((1 - ci)/2) * yse
  colnames(yse) <- names(object$varresult)
  
  l <- nrow(yse)
  rbind(yse[1:(l - 1),], yse[rep(l, n_pred - l - object$p + 1),])
}


create_cf_1 <- function(mod, nnn = NULL, ...) {
  exogen <<- create_exogen(data, include, week, slice = n_vac:n) %>%
    `[`(names(.) %in% colnames(mod$datamat)) #super assign to solve {vars} bug
  
  result <- custom_predict(mod, n.ahead = n - n_vac + 1, dumvar = exogen, ..., nnn = nnn)
  
  rmv <- grep("vaccines", names(result$fcst))
  if (length(rmv) == 1) result$fcst[[rmv]] <- NULL
  result
}

create_cf_2 <- function(mod, ...) {
  exogen <<- create_exogen(data, include, week, slice = (n_vac - mod$p):n) %>%
    `[`(names(.) %in% colnames(mod$datamat)) #super assign to solve {vars} bug
  data <- data %>% slice((n_vac - mod$p):n)
  
  pred <- map_dfc(1:mod$p, function(x) { #create data for "prediction"
    to <- nrow(data) + x - 1
    map_dfc(select(data, c(cases, deaths)), ~ c(.x, numeric(x))[x:to]) %>%
      set_names(paste0(names(.), ".l", x))
  }) %>%
    tibble(
      const = if (type %in% c("const", "both")) 1,
      trend = if (type %in% c("trend", "both")) mod$p:(nrow(.) + mod$p - 1),
      exogen
    ) %>%
    as.matrix()
  
  coefs <- mod$varresult$deaths$coefficients %>% .[!grepl("vaccines", names(.))]
  cols <- grepl("deaths", colnames(pred))
  
  for (t in 1:nrow(pred)) {
    values <- pred[t,] %>% .[names(.) %in% names(coefs)]
    pred[t, cols] <- c(sum(coefs * values), pred[t, cols][-mod$p])
  }
  
  yse <- get_yse(mod, n.ahead, n_pred = nrow(pred), ...) %>% .[,"deaths"]
  
  result <- list(
    fcst = list(
      deaths = as.matrix(tibble(
        fcst = pred[-(1:mod$p), "deaths.l1"],
        lower = fcst - yse, upper = fcst + yse, CI = yse
      ))
    ),
    model = mod
  )
  
  class(result) <- "varprd"
  result
}


custom_predict <- function(object, ..., n.ahead = 10, ci = 0.95, dumvar = NULL, nnn) {
  K <- object$K
  p <- object$p
  obs <- object$obs
  type <- object$type
  data.all <- object$datamat
  ynames <- colnames(object$y)
  n.ahead <- as.integer(n.ahead)
  Z <- object$datamat[, -c(1:K)]
  B <- Bcoef(object)
  if (type == "const") {
    Zdet <- matrix(rep(1, n.ahead), nrow = n.ahead, ncol = 1)
    colnames(Zdet) <- "const"
  }
  else if (type == "trend") {
    trdstart <- nrow(Z) + 1 + p
    Zdet <- matrix(seq(trdstart, length = n.ahead), nrow = n.ahead, 
                   ncol = 1)
    colnames(Zdet) <- "trend"
  }
  else if (type == "both") {
    trdstart <- nrow(Z) + 1 + p
    Zdet <- matrix(c(rep(1, n.ahead), seq(trdstart, length = n.ahead)), 
                   nrow = n.ahead, ncol = 2)
    colnames(Zdet) <- c("const", "trend")
  }
  else if (type == "none") {
    Zdet <- NULL
  }
  if (!is.null(eval(object$call$season))) {
    season <- eval(object$call$season)
    seas.names <- paste("sd", 1:(season - 1), sep = "")
    cycle <- tail(data.all[, seas.names], season)
    seasonal <- as.matrix(cycle, nrow = season, ncol = season - 
                            1)
    if (nrow(seasonal) >= n.ahead) {
      seasonal <- as.matrix(cycle[1:n.ahead, ], nrow = n.ahead, 
                            ncol = season - 1)
    }
    else {
      while (nrow(seasonal) < n.ahead) {
        seasonal <- rbind(seasonal, cycle)
      }
      seasonal <- seasonal[1:n.ahead, ]
    }
    rownames(seasonal) <- seq(nrow(data.all) + 1, length = n.ahead)
    if (!is.null(Zdet)) {
      Zdet <- as.matrix(cbind(Zdet, seasonal))
    }
    else {
      Zdet <- as.matrix(seasonal)
    }
  }
  if (!is.null(eval(object$call$exogen))) {
    if (is.null(dumvar)) {
      stop("\nNo matrix for dumvar supplied, but object varest contains exogenous variables.\n")
    }
    if (!all(colnames(dumvar) %in% colnames(data.all))) {
      stop("\nColumn names of dumvar do not coincide with exogen.\n")
    }
    if (!identical(nrow(dumvar), n.ahead)) {
      stop("\nRow number of dumvar is unequal to n.ahead.\n")
    }
    if (!is.null(Zdet)) {
      Zdet <- as.matrix(cbind(Zdet, dumvar))
    }
    else {
      Zdet <- as.matrix(dumvar)
    }
  }
  Zy <- as.matrix(object$datamat[, 1:(K * (p + 1))])
  yse <- matrix(NA, nrow = n.ahead, ncol = K)
  sig.y <- vars:::.fecov(x = object, n.ahead = n.ahead)
  for (i in 1:n.ahead) {
    yse[i, ] <- sqrt(diag(sig.y[, , i]))
  }
  yse <- -1 * qnorm((1 - ci)/2) * yse
  colnames(yse) <- paste(ci, "of", ynames)
  forecast <- matrix(NA, ncol = K, nrow = n.ahead)
  lasty <- c(Zy[nnn %||% nrow(Zy), ])
  for (i in 1:n.ahead) {
    lasty <- lasty[1:(K * p)]
    Z <- c(lasty, Zdet[i, ])
    forecast[i, ] <- B %*% Z
    temp <- forecast[i, ]
    lasty <- c(temp, lasty)
  }
  colnames(forecast) <- paste(ynames, ".fcst", sep = "")
  lower <- forecast - yse
  colnames(lower) <- paste(ynames, ".lower", sep = "")
  upper <- forecast + yse
  colnames(upper) <- paste(ynames, ".upper", sep = "")
  forecasts <- list()
  for (i in 1:K) {
    forecasts[[i]] <- cbind(forecast[, i], lower[, i], upper[, 
                                                             i], yse[, i])
    colnames(forecasts[[i]]) <- c("fcst", "lower", "upper", 
                                  "CI")
  }
  names(forecasts) <- ynames
  result <- list(fcst = forecasts, endog = object$y, model = object, 
                 exo.fcst = dumvar)
  class(result) <- "varprd"
  return(result)
}


custom_svars_fun <- function(svars_fun, x, ...) {
  keep <- if (x$type == "both") 1:2 else 1
  x$A_hat <- x$A_hat[, c(keep, (ncol(x$A_hat) - x$K*x$p + 1):(ncol(x$A_hat)))]
  svars_fun(x, ...)
}


custom_boot <- function(x, scales = "free_y", lowerq = 0.16, upperq = 0.84, 
                        percentile = "standard", selection = NULL, cumulative = NULL, 
                        ..., base) {
  probs <- NULL
  V1 <- NULL
  value <- NULL
  cutt <- function(x, selind) {
    x$irf <- x$irf[selind]
    return(x)
  }
  agg <- function(x, aggind) {
    x$irf[, aggind] <- cumsum(x$irf[aggind])
    return(x)
  }
  kk <- ncol(x$true$irf)
  k <- sqrt(kk - 1)
  if (!is.null(cumulative)) {
    aggind <- list()
    temp <- (k * cumulative) + 1
    for (i in 1:length(temp)) {
      aggind[[i]] <- temp[i] - (k - c(1:k))
    }
    aggind <- unlist(aggind)
    x$true <- agg(x$true, aggind)
    x$bootstrap <- lapply(x$bootstrap, agg, aggind)
  }
  if (!is.null(selection)) {
    selind <- list()
    temp <- (k * selection[[1]]) + 1
    for (i in 1:length(temp)) {
      selind[[i]] <- temp[i] - (k - selection[[2]])
    }
    selind <- c(1, unlist(selind))
    x$true <- cutt(x$true, selind)
    x$bootstrap <- lapply(x$bootstrap, cutt, selind)
    kk <- ncol(x$true$irf)
  }
  n.ahead <- nrow(x$true$irf)
  bootstrap <- x$bootstrap
  nboot <- length(bootstrap)
  rest <- x$rest_mat
  n.probs <- length(lowerq)
  if (length(lowerq) != length(upperq)) {
    stop("Vectors 'lowerq' and 'upperq' must be of same length!")
  }
  intervals <- array(0, c(n.ahead, kk, nboot))
  for (i in 1:nboot) {
    intervals[, , i] <- as.matrix(bootstrap[[i]]$irf)
  }
  lower <- array(0, dim = c(n.ahead, kk, n.probs))
  upper <- array(0, dim = c(n.ahead, kk, n.probs))
  if (percentile == "standard" | percentile == "hall") {
    for (i in 1:n.ahead) {
      for (j in 1:kk) {
        lower[i, j, ] <- quantile(intervals[i, j, ], 
                                  probs = lowerq)
        upper[i, j, ] <- quantile(intervals[i, j, ], 
                                  probs = upperq)
      }
    }
    if (percentile == "hall") {
      for (p in 1:n.probs) {
        lower[, -1, p] <- as.matrix(2 * x$true$irf)[, 
                                                    -1] - lower[, -1, p]
        upper[, -1, p] <- as.matrix(2 * x$true$irf)[, 
                                                    -1] - upper[, -1, p]
      }
    }
  }
  else if (percentile == "bonferroni") {
    rest <- matrix(t(rest), nrow = 1)
    rest[is.na(rest)] <- 1
    rest <- c(1, rest)
    for (i in 1:n.ahead) {
      for (j in 1:kk) {
        if (rest[j] == 0) {
          lower[i, j, ] <- quantile(intervals[i, j, 
          ], probs = (lowerq/n.ahead))
          upper[i, j, ] <- quantile(intervals[i, j, 
          ], probs = 1 + (((upperq - 1)/n.ahead)))
        }
        else {
          lower[i, j, ] <- quantile(intervals[i, j, 
          ], probs = (lowerq/(n.ahead + 1)))
          upper[i, j, ] <- quantile(intervals[i, j, 
          ], probs = 1 + (((upperq - 1)/(n.ahead + 
                                           1))))
        }
      }
    }
  }
  else {
    stop("Invalid choice of percentile; choose between standard, hall and bonferroni")
  }
  alp <- 0.7 * (1 + log(n.probs, 10))/n.probs
  irf <- reshape2::melt(x$true$irf, id = "V1")
  cbs <- data.frame(V1 = rep(irf$V1, times = n.probs),
                    variable = rep(irf$variable, times = n.probs),
                    probs = rep(1:n.probs, each = (kk - 1) * n.ahead),
                    lower = c(lower[, -1, ]), upper = c(upper[, -1, ]))
  ggplot() +
    geom_ribbon(data = cbs, aes(x = V1, ymin = lower, ymax = upper, group = probs),
                alpha = alp, fill = NA, linetype = 2, color = pal[1]) + 
    geom_line(data = irf, aes(x = V1, y = value)) +
    geom_hline(yintercept = 0, color = "black") +
    facet_wrap(~variable, scales = scales, labeller = label_parsed) +
    xlab("Horizon") + ylab("Response") + 
    theme_bw()
}


custom_irfs <- function(x, base, scales = "free_y", selection = NULL, cumulative = NULL, 
                        ...) {
  cutt <- function(x, selind) {
    x$irf <- x$irf[selind]
    return(x)
  }
  agg <- function(x, aggind) {
    x$irf[, aggind] <- cumsum(x$irf[aggind])
    return(x)
  }
  kk <- ncol(x$irf)
  k <- sqrt(kk - 1)
  if (!is.null(cumulative)) {
    aggind <- list()
    temp <- (k * cumulative) + 1
    for (i in 1:length(temp)) {
      aggind[[i]] <- temp[i] - (k - c(1:k))
    }
    aggind <- unlist(aggind)
    x <- agg(x, aggind)
  }
  if (!is.null(selection)) {
    selind <- list()
    temp <- (k * selection[[1]]) + 1
    for (i in 1:length(temp)) {
      selind[[i]] <- temp[i] - (k - selection[[2]])
    }
    selind <- c(1, unlist(selind))
    x <- cutt(x, selind)
  }
  impulse <- reshape2::melt(x$irf, id = "V1")
  ggplot(impulse, aes_(x = ~V1, y = ~value)) + geom_line() + 
    geom_hline(yintercept = 0, color = "black") + facet_wrap(~variable, 
                                                             scales = scales, labeller = label_parsed) + xlab("Horizon") + 
    ylab("Response") + theme_bw()
}
