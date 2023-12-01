# Covid SVAR Diff
This is the repo for the revision of my graduation dissertation "Efeito Agregado das Vacinas do COVID-19 no Brasil, uma metodologia SVAR".


## Abstract
Finding evidence for the effectiveness of COVID vaccines is very important to guide the action of policymakers, especially considering how recent the pandemic is. This thesis contributes to literature by using an SVAR, motivated by and identifyed by a SIR epidemiological model. The capture of temporal dynamics between variables, the relaxation of assumptions on infection, mortality, and vaccination rates constant over time, as well as the potential of generating different types of results are this method’s main benefits. It was discovered that the response of deaths to a shock of 1000 cases was 81\% (or 18.21, in absolute terms) less than after the inclusion of vaccines. There was also an analysis of counterfactuals, indicating that the number of lives saved goes from \textit{lower bound} of 31 thousand to a \textit{upper bound} of 281 thousand (but there is indication that this band is substantially less wide).

**Keywords:** Coronavirus. Vaccination. Causal Effect. SVAR


## Directories
The "Latex" directory has the .tex files to generate the paper .pdf

The "R (covid-svardiff)" has the R files related to the empirical (and reproducible) analisys


## Future Changes
There are still some work to do. Besides what I list as "possible future expansions" in the text, which I won't be doing in this paper, there are still some improvements I want to do in the near future:

- Add bootstraps for the IRF's;
- Study using more traditional SVAR counterfactual methos, like the one used by the R function `svars::cf`;
- Update the diagnostics appendix;
- Create an package with the "method"'s functions, and reorganize the code.
