get_basic_data <- function(data_raw, trans, n_end = max(data_raw$date)){
  data_raw %>%
    transmute(
      date = date,
      cases_level = new_cases_smoothed,
      deaths_level = new_deaths_smoothed,
      vaccines_level = new_vaccinations_smoothed,
      cases_prop = new_cases_smoothed_per_million,
      deaths_prop = new_deaths_smoothed_per_million,
      vaccines_prop = new_vaccinations_smoothed_per_million
    ) %>%
    mutate(across(ends_with("_level"),
      ~ ifelse(.x == 0, 0, log(.x)), .names = "{str_remove(.col, '_level')}_log"
    )) %>%
    select(c(date, matches(paste0("_", trans, "$")))) %>%
    rename_with(.cols = -date, ~ str_remove(.x, "_.+$")) %>%
    filter(!if_all(-date, is.na) & date <= n_end & cases + deaths > 0) %>%
    na_approx()
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
  breaks <- unlist(breaks)
  
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


select_var <- function(data, lag.max, slice = 1:nrow(data), remove = matches("$^"), include, ...) {
  data <- data %>% slice(slice)
  exogen <- create_exogen(data, include, week)
  select(data, -c(date, one_of(remove))) %>% VARselect(lag.max, exogen = exogen, ...)
}


create_var <- function(data, p, slice = 1:nrow(data), remove = matches("$^"), include, ...) {
  data <- data %>% slice(slice)
  exogen <- create_exogen(data, include, week)
  select(data, -c(date, one_of(remove))) %>% VAR(p = p, exogen = exogen, ...)
}


sci <- function(x) format(x, digits = 1, scientific = TRUE)


custom_ggsave = function(filename, w, h, ...){
  ggsave(paste0("figures/", filename),
         device = "png", width = w, height = h, units = "cm", ...
  )
}
