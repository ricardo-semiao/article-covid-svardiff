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


create_exogen <- function(data, include = "none") {
  exogen <- list()
  
  if (include == "none") {
    return(NULL)
  }
  
  if (grepl("season", include)) {
    exogen$season <- data %>%
      transmute(
        season_high = format(date, "%m") %in% c("12", "01", "02"),
        season_mid = format(date, "%m") %in% c("10", "11", "03", "04")
      )
  }
  
  if (grepl("intramonth", include)) {
    exogen$intramonth <- as.numeric(format(data$date, "%U")) %% 4 %>%
      fastDummies::dummy_cols(remove_first_dummy = TRUE, remove_selected_columns = TRUE) %>%
      set_names(c("week1", "week2", "week3"))
  }
  
  bind_cols(exogen)
}

create_var <- function(data, p, slice = 1:nrow(data), remove_match = "$^", include) {
  data <- data %>% slice(slice)
  exogen <- create_exogen(data, include)
  
  data %>%
    select(-matches(remove_match)) %>%
    VAR(p = p, exogen = exogen)
}


#custom_svars_irf <- function(x, ..., n.ahead = 20) 
#{
#  if (!(is(x, "svars"))) {
#    stop("\nPlease provide an object of class 'svars'.\n")
#  }
#  if (x$type == "const") {
#    A_hat <- x$A_hat[, -1]
#  }
#  else if (x$type == "trend") {
#    A_hat <- x$A_hat[, -1]
#  }
#  else if (x$type == "both") {
#    A_hat <- x$A_hat[, -c(1, 2)]
#  }
#  else {
#    A_hat <- x$A_hat
#  }
#  A_hat <- A_hat[, (ncol(A_hat) - x$K*x$p + 1):(ncol(A_hat))]
#  B_hat <- x$B
#  IR <- array(unlist(svars:::IRF(A_hat, B_hat, n.ahead)), c(x$K, x$K, 
#                                                    n.ahead))
#  impulse <- matrix(0, ncol = dim(IR)[2]^2 + 1, nrow = dim(IR)[3])
#  colnames(impulse) <- rep("V1", ncol(impulse))
#  cc <- 1
#  impulse[, 1] <- seq(0, dim(IR)[3] - 1)
#  for (i in 1:dim(IR)[2]) {
#    for (j in 1:dim(IR)[2]) {
#      cc <- cc + 1
#      impulse[, cc] <- IR[i, j, ]
#      colnames(impulse)[cc] <- paste("epsilon[", colnames(x$y)[j], 
#                                     "]", "%->%", colnames(x$y)[i])
#    }
#  }
#  impulse <- list(irf = as.data.frame(impulse))
#  class(impulse) <- "svarirf"
#  return(impulse)
#}
