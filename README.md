# Covid SVAR Diff
This is the repo for the revision of my graduation dissertation "Efeito Agregado das Vacinas do COVID-19 no Brasil, uma metodologia SVAR".

The finished article is [text.pdf](text.pdf). The data comes from Our World in Data and is at [Data/](Data/). The main script is [main.R](main.R), and the article is built with [text.tex](text.tex).


## Abstract
Finding evidence for the effectiveness of COVID vaccines is very important to guide the action of policymakers, especially considering how recent the pandemic is. This thesis contributes to literature by using an SVAR, motivated by and identified by a SIR epidemiological model. The capture of temporal dynamics between variables, the relaxation of assumptions on infection, mortality, and vaccination rates constant over time, as well as the potential of generating different types of results are this method’s main benefits. It was discovered that the response of deaths to a shock of 1000 cases was 81\% (or 18.21, in absolute terms) less than after the inclusion of vaccines. There was also an analysis of counterfactuals, indicating that the number of lives saved goes from lower bound of 31 thousand to a upper bound of 281 thousand (but there is indication that this band is substantially less wide).

**Keywords:** Coronavirus. Vaccination. Causal Effect. SVAR.


## Future Changes
There are still some work to do. Besides what I list as "possible future expansions" in the text, which I won't be doing in this paper, there are still some improvements I want to do in the near future:

- Add bootstraps for the IRF's;
- Study using more traditional SVAR counterfactual methods, like the one used by the R function `svars::cf`;
- Update the diagnostics appendix;
- Create an package with the "method"'s functions, and reorganize the code.
