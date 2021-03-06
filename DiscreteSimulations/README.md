# Discrete Comparison Function

We have implemented three algorithms:

BRL- Bayesian Record Linkage described in Sadinle, Mauricio. "Bayesian estimation of bipartite matchings for record linkage." Journal of the American Statistical Association 112.518 (2017): 600-612.

CIBRL - Conditionally independent record linkage.

MLBRL - Multilayered recrod linkage application.

Each of the algrothims is an R script that can be implemented with three error parameters **Regionerrorprob, IncomeError, and DOBerrorprob**. These probabilities represent the error in he variables region, income and date of birth, respectively. They are entered as a value between 0-10.

Here is the command line to run the scripts for no error in region variable, 20% error in income variable, and 40% error in date of birth:

BRL - Rscript --verbose BRL_Disc_noDOBDay.R 0 2 4
CIBRL - Rscript --verbose CIBRL_Disc_noDOBDay.R 0 2 4
MLBRL - Rscript --verbose CIBRL_Disc_noDOBDay.R 0 2 4
