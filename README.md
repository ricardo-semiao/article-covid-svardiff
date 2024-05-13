# Article: Aggregated Effect of COVID-19 Vacination in Brazil, a SVAR Method

Welcome! This is the repository for the revision of my graduation dissertation "Efeito Agregado das Vacinas do COVID-19 no Brasil, uma metodologia SVAR". The article is in portuguese, but i present the english abstract below. The finished article is [text.pdf](text.pdf). The data comes from Our World in Data and is at [data/](data/). The main script is [main.R](main.R), and the article is built with [text.tex](text.tex).

Bem-vindo! Este é o repositório da revisão da minha dissertação de graduação "Efeito Agregado das Vacinas do COVID-19 no Brasil, uma metodologia SVAR". O artigo finalizado está em [text.pdf](text.pdf). Os dados são provenientes do Our World in Data e estão em [data/](data/). O script principal é [main.R](main.R), e o artigo é construído com [text.tex](text.tex).


## Resumo

Encontrar evidências para a efetividade das vacinas do COVID-19 é muito importante para guiar a ação de \textit{policymakers}, ainda mais dado o quão recente a pandemia é. Esta tese contribui com essa literatura utilizando uma \textit{Structural Vector Auto-Regression} para modelar a relação entre as variáveis, motivada e identificada por um modelo epidemiológio SIR. A captura das dinâmicas temporais entre as variáveis, o relaxamento de hipóteses de taxas de infecção, mortalidade, e vacinação constantes no tempo, e o potencial de gerar diferentes tipos de resultados, são os benefícios principais desse método. Foi descoberto que resposta de mortes à um choque de $1000$ casos foi $60$\% (ou $10,24$, em termos absolutos) menor depois da inclusão das vacinas. Também foi feita uma análise de contrafactuais, indicando que o número de vidas salvas vai de um \textit{lower bound} de $16$ mil a um \textit{upper bound} de $279$ mil, com indícios de maior proximidade do segundo.

**Palavras-chave:** Coronavirus. Vacinação. Efeito Causal. SVAR.


## Abstract
Finding evidence for the effectiveness of COVID vaccines is very important to guide the action of policymakers, especially considering how recent the pandemic is. This thesis contributes to literature by using an SVAR, motivated by and identified by a SIR epidemiological model. The capture of temporal dynamics between variables, the relaxation of assumptions on infection, mortality, and vaccination rates constant over time, as well as the potential of generating different types of results are this method’s main benefits. It was discovered that the response of deaths to a shock of 1000 cases was 81\% (or 18.21, in absolute terms) less than after the inclusion of vaccines. There was also an analysis of counterfactuals, indicating that the number of lives saved goes from lower bound of 31 thousand to a upper bound of 281 thousand (but there is indication that this band is substantially less wide).

**Keywords:** Coronavirus. Vaccination. Causal Effect. SVAR.


## Mudanças Futuras
Ainda há trabalho a ser feito. Além do que listo como "possíveis expansões futuras" no texto, que não farei neste artigo, ainda há algumas melhorias que desejo fazer em um futuro próximo:

- Adicionar bootstraps para os IRF's;
- Estudar o uso de métodos de contrafactuais mais tradicionais do SVAR, como o utilizado pela função R `svars::cf`;
- Atualizar o apêndice de diagnósticos;
- Criar um pacote com as funções do "método" e reorganizar o código.
