get_basic_data <- function(data_raw, trans){
  data_basic <- data_raw %>%
    transmute(
      date = date,
      cases_s = new_cases_smoothed,
      deaths_s = new_deaths_smoothed,
      vaccines_s = new_vaccinations_smoothed,
      cases_p = new_cases_smoothed_per_million,
      deaths_p = new_deaths_smoothed_per_million,
      vaccines_p = new_vaccinations_smoothed_per_million
    ) %>%
    filter(
      !if_all(-date, is.na) & date < as.Date("2023-03-23") & cases_s + deaths_s > 0
    ) %>%
    na_approx() %>%
    mutate(across(ends_with("_s"), ~ ifelse(.x == 0, 0, log(.x)), .names = "{str_remove(.col, '_s')}_l"))
  
  switch(trans,
         "level" = select(data_basic, !matches("_p$|_l$")),
         "prop"  = select(data_basic, !matches("_s$|_l$")),
         "log"   = select(data_basic, !matches("_p$|_s$"))
  )
}


na_approx <- function(data) {
  numeric_cols <- sapply(data, is.numeric)
  data[numeric_cols] <- map_dfc(data[numeric_cols], function(col) {
    index <- (which(!is.na(col))[1]):length(col)
    col[index] <- zoo::na.approx(col[index], na.rm = FALSE)
    col
  })
  data
}


create_exogen <- function(data, include = "none", week, breaks = vars_breaks, slice = NULL) {
  exogen <- list()
  
  if (include == "none") {
    return(NULL)
  }
  
  if (grepl("season", include)) {
    exogen$season <- data %>%
      transmute(
        season_high = format(date, "%m") %in% c("12", "01", "02"),
        season_mid = format(date, "%m") %in% c("10", "11", "03", "04")
      ) %>%
      mutate(across(everything(), as.numeric))
  }
  
  if (grepl("intramonth", include)) {
    exogen$intramonth <- as.numeric(format(data$date, "%U")) %% 4 %>%
      fastDummies::dummy_cols(
        remove_first_dummy = TRUE,
        remove_selected_columns = TRUE,
        ignore_na = TRUE
      ) %>%
      set_names(c("week1", "week2", "week3"))
  }
  
  if (grepl("waves", include)) {
    x <- findInterval(data$date, breaks)
    x <- ifelse(x %in% c(1, 3, 5, 7, 9), x, 0)
    
    inds <- lapply(setdiff(unique(x), 0), \(u) x == u)
    types <- names(breaks) %>% .[. != "none"]
    
    mat <- NULL
    for (i in seq_along(inds)) {
      ind <- inds[[i]]
      temp <- numeric(length(x))
      if (types[i] == "const") {
        temp[ind] <- rep(1, sum(ind))
        add <- tibble("wave{i}.const" := temp)
      } 
      if (types[i] %in% c("lin", "sqr")) {
        temp[ind] <- 1:sum(ind)
        add <- tibble("wave{i}.lin" := temp)
      } 
      if (types[i] == "sqr") {
        add <- tibble("wave{i}.lin" := temp, "wave{i}.sqr" := temp^2)
      } 
      mat <- bind_cols(mat, add)
    }
    exogen$waves <- mat
  }
  
  if (week && grepl("intraweek", include)) {
    x <- as.numeric(data$date %>% format("%w"))
    exogen$weekdays <- ifelse(x %in% c(6, 0),
      "weekend", ifelse(x %in% c(1, 2), "weekstart", "weekmid")
    ) %>%
      fastDummies::dummy_cols(remove_selected_columns = TRUE) %>%
      select(-weekmid)
  }
  if (length(exogen) == 0) {
    NULL
  } else {
    bind_cols(exogen) %>%
      slice(slice %||% 1:nrow(.))
  }
}


select_var <- function(data, lag.max, slice = 1:nrow(data), remove_match = "$^", include, ...) {
  data <- data %>% slice(slice)
  exogen <- create_exogen(data, include, week)
  
  data %>%
    select(-matches(remove_match)) %>%
    VARselect(lag.max, exogen = exogen, ...)
}


create_var <- function(data, p, slice = 1:nrow(data), remove_match = "$^", include, ...) {
  data <- data %>% slice(slice)
  exogen <- create_exogen(data, include, week)
  
  data %>%
    select(-matches(remove_match)) %>%
    VAR(p = p, exogen = exogen, ...)
}


custom_svars_fun <- function(svars_fun, x, ...) {
  keep <- if (x$type == "both") 1:2 else 1
  x$A_hat <- x$A_hat[, c(keep, (ncol(x$A_hat) - x$K*x$p + 1):(ncol(x$A_hat)))]
  svars_fun(x, ...)
}

name <- function(x, obj = NULL, value = TRUE) {
  names <- if (is.null(obj)) {
    vars_order
  } else if (is.matrix(obj)) {
    colnames(obj)
  } else {
    names(obj)
  }
  grep(x, names, value = value)
} 

custom_predict <- function (object, ..., n.ahead = 10, ci = 0.95, dumvar = NULL, nnn) {
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
