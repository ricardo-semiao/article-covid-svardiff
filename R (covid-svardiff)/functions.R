table_adf <- function(data) {
  result_adf <- data %>%
    filter(date >= date_vaccine$week) %>%
    select(ends_with("l")) %>%
    map_dfr(function(col) {
      col %>%
        tseries::adf.test() %>%
        `[`(c("statistic", "parameter", "p.value"))
    })
  
  result_adf %>%
    set_names(c("Estatística", "Lag", "P-valor")) %>%
    mutate(across(-Lag, ~ formatC(.x, digits = 2, format = "f"))) %>%
    mutate(Variável = vars_names, .before = 1) %>%
    stargazer(summary = FALSE, title = "Teste Dickey-Fuller Aumentado", label = "tb:testadf")
}

pvalue.ast = function(x) {
  c("", "*", "**", "***")[pmap_dbl(tibble(x), ~ sum(. < c(0.01, 0.05, 0.10)) + 1)]
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

ccf_ba <- function(data_list, labels) {
  data_ccf <- imap_dfr(data_list, function(df, name) {
    confInterval <- qnorm((1 - 0.95)/2)/sqrt(nrow(df))
    ccf_values <- ccf(df$cases_s, df$deaths_s, lag.max = (start_vaccine), plot = FALSE)
    tibble(
      name = as.factor(name), Lag = ccf_values$lag[,,1],
      value = ccf_values$acf[,,1], ci = confInterval
    )
  })
  
  ggplot(data_ccf, aes(x = Lag, y = value)) +
    geom_area() + 
    geom_ribbon(aes(ymin = -ci, ymax = ci), linetype = 2, fill = NA) +
    geom_hline(yintercept = 0) +
    geom_vline(xintercept = 0) +
    facet_wrap(vars(name), ncol = 3, labeller = labeller(
      name = set_names(labels, c("Before", "After", "Afterp"))
    )) +
    labs(y = "Valor", x = "Lag (casos)")
}

granger_ba = function(data, p, test, ...){
  data <- map(data, ~ select(.x, !starts_with("vaccines")))
  mods <- map(data, ~VAR(.x, p))
  mods_vcov <- map(mods, ~vcovHC(.x, "HC1"))
  
  data_granger <- map2(mods, mods_vcov, function(mod, vcov) {
    causality(mod, cause = "cases_s", vcov. = vcov)[[test]] %>% 
      `[`(c("statistic", "parameter", "p.value")) %>%
      reduce(c)
  }) %>%
    reduce(rbind) %>%
    as_tibble()
  
  data_granger %>%
    setNames(c("Estatística", "GL 1", "GL 2", "P-valor")) %>%
    mutate(
      Estatística = round(Estatística, 2),
      across(contains("GL"), round),
      `P-valor` = paste(round(`P-valor`, 4), pvalue.ast(`P-valor`))
    ) %>%
    mutate(Variável = c("Antes", "Depois", "Depois+"), .before = 1) %>%
    stargazer(
      summary = FALSE, title = "Causalidade de Granger",
      label = "tb:grangerba", notes = "GL: graus de liberdade",
      rownames = FALSE, ...
    )
}
