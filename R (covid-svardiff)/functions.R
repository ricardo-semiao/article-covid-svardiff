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


create_exogen <- function(data, include = "none", week) {
  exogen <- list()
  breaks_waves <- as.Date(c(
    start1 = "2021-02-21", end1 = "2021-06-22",
    start2 = "2022-01-03", end2 = "2022-03-12"
  ))
  
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
    x <- findInterval(data$date, breaks_waves)
    l <- length(unique(x))
    if (l > 1) {
      exogen$waves <- x %>%
        fastDummies::dummy_cols(remove_first_dummy = TRUE, remove_selected_columns = TRUE) %>%
        set_names(paste0("wave", 2:l))
    }
  }
  
  if (week && grepl("intraweek", include)) {
    x <- as.numeric(data$date %>% format("%w"))
    exogen$weekdays <- ifelse(x %in% c(6, 0),
      "weekend", ifelse(x %in% c(1, 2), "weekstart", "weekmid")
    ) %>%
      fastDummies::dummy_cols(remove_selected_columns = TRUE) %>%
      select(-weekmid)
  }
  
  bind_cols(exogen)
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
