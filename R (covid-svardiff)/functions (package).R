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
    geom_ribbon(aes(ymin = -ci, ymax = ci), linetype = 2, fill = NA, color = pal[1]) +
    facet_wrap(vars(name), ncol = 3, labeller = labeller(
      name = set_names(labels, c("Before", "After", "Afterp"))
    )) +
    labs(y = "Correlação", x = "Lag (casos)")
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


ggvar_select <- function(mod_select, lag.max) {
  mod_select %>%
    imap_dfr(function(sel, mod) {
      crit <- as_tibble(apply(t(sel$criteria), 2, \(x) x/x[1]))
      tibble(lag = 1:lag.max, mod = mod, select(crit, -`FPE(n)`))
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
    `[`(names(.) %in% colnames(mod$datamat))
  
  result <- custom_predict(mod, n.ahead = n - n_vac + 1, dumvar = exogen, ..., nnn = nnn)
  
  rmv <- grep("vaccines", names(result$fcst))
  if (length(rmv) == 1) result$fcst[[rmv]] <- NULL
  result
}

create_cf_2 <- function(mod, ...) {
  exogen <<- create_exogen(data, include, week, slice = (n_vac - mod$p):n) %>%
    `[`(names(.) %in% colnames(mod$datamat))
  data <- data %>% slice((n_vac - mod$p):n)
  
  pred <- map_dfc(1:mod$p, function(x) {
    to <- nrow(data) + x - 1
    map_dfc(data[grep("cases|deaths", names(data))], ~ c(.x, numeric(x))[x:to]) %>%
      set_names(paste0(names(.), ".l", x))
  }) %>%
    tibble(
      const = if (type %in% c("const", "both")) 1,
      trend = if (type %in% c("trend", "both")) mod$p:(nrow(.) + mod$p - 1),
      exogen
    ) %>%
    as.matrix()
  
  coefs <- mod$varresult[[name("death", data)]]$coefficients %>% .[!grepl("vaccines", names(.))]
  cols <- grepl("deaths", colnames(pred))
  
  for (t in 1:nrow(pred)) {
    values <- pred[t,] %>% .[names(.) %in% names(coefs)]
    pred[t, cols] <- c(sum(coefs * values), pred[t, cols][-mod$p])
  }
  
  yse <- get_yse(mod, n.ahead, n_pred = nrow(pred), ...) %>% .[,grep("death", colnames(.))]
  
  result <- list(
    fcst = list(
      as.matrix(tibble(
        fcst = pred[-(1:mod$p), grep("death.+l1$", colnames(pred))],
        lower = fcst - yse, upper = fcst + yse, CI = yse
      ))
    ) %>% set_names(grep("death", vars_order, value = TRUE)),
    model = mod
  )
  
  class(result) <- "varprd"
  result
}
