table_adf <- function(data, start, vars_names, title = "Teste Dickey-Fuller Aumentado",
  table_names = c("Variáveis", "Estatística", "Lag", "P-valor"), ...) {
  result_adf <- data %>%
    filter(date >= start) %>%
    select(-date) %>%
    map_dfr(function(col) {
      col %>%
        na.omit() %>%
        tseries::adf.test() %>%
        `[`(c("statistic", "parameter", "p.value"))
    })
  
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
    confInterval <- qnorm((1 - 0.95)/2)/sqrt(nrow(df))
    ccf_values <- ccf(df[grepl("cases", names(df))], df[grepl("deaths", names(df))], plot = FALSE, ...)
    tibble(
      name = as.factor(name),
      Lag = ccf_values$lag[,,1],
      value = ccf_values$acf[,,1],
      ci = confInterval
    )
  })
  
  ggplot(data_ccf, aes(x = Lag, y = value)) +
    geom_vline(xintercept = 0, color = "darkgray") +
    geom_hline(yintercept = 0) +
    geom_segment(aes(xend = Lag, yend = 0)) + 
    geom_ribbon(aes(ymin = -ci, ymax = ci), linetype = 2, fill = NA, color = "blue") +
    facet_wrap(vars(name), ncol = 3, labeller = labeller(
      name = set_names(labels, c("Before", "After", "Afterp"))
    )) +
    labs(y = "Valor", x = "Lag (casos)")
}


granger_ba = function(data, p, test, ...){
  pvalue.ast = function(x) {
    c("", "*", "**", "***")[pmap_dbl(tibble(x), ~ sum(. < c(0.01, 0.05, 0.10)) + 1)]
  }
  
  data <- map(data, ~ select(.x, !starts_with("vaccines")))
  mods <- map(data, ~VAR(.x, p))
  mods_vcov <- map(mods, ~vcovHC(.x, "HC1"))
  
  data_granger <- map2(mods, mods_vcov, function(mod, vcov) {
    causality(mod, cause = grep("cases", names(mod$varresult), value = TRUE), vcov. = vcov)[[test]] %>% 
      `[`(c("statistic", "parameter", "p.value")) %>%
      unlist()
  }) %>%
    reduce(rbind) %>%
    as_tibble() %>%
    set_names(if(test == "Granger") {
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
    mutate(Variável = c("Antes", "Depois", "Depois+"), .before = 1) %>%
    stargazer(
      summary = FALSE, title = "Causalidade de Granger",
      label = "tb:grangerba", notes = "GL: graus de liberdade",
      rownames = FALSE, ...
    )
}


ggvar_irf2 <- function(data, ...) {
  ggplot(data, aes(horizon, irf, ymin = Lower, ymax = Upper)) +
    geom_line() +
    geom_ribbon(fill = NA, color = "blue", linetype = 2) +
    geom_hline(yintercept = 0) +
    facet_wrap(vars(model), ncol = 3, ...) +
    labs(title = "", y = "Resposta")
}
