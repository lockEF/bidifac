# BIDIFAC
Bidimensional matrix factorization

R code to perform factorization and dimension reduction of bidimensionally linked data matrices [1].  The file BIDIFAC.R gives functions to perform the BIDIFAC factorization, simulate data, and perform missing value imputation using BIDIFAC.  See BIDIFAC_example.R for example usage.  Written by [Jun Young Park](https://www.statisticspark.com/).  

[1] J Park and EF Lock. Integrative Factorization of Bidimensionally Linked Matrices. <em> Biometrics</em>, 2019. [(article link)](https://onlinelibrary.wiley.com/doi/abs/10.1111/biom.13141) [arXiv preprint: 1906.03722](https://arxiv.org/abs/1906.03722).

The functions bidifac.plus.R, bidifac.plus.impute.R, and bidifac.plus.given.R perform the BIDIFAC+ method, described in [2].  See bidifac.plus_example.R for example usage.    

[2] EF Lock, J Park and KA Hoadley.  Bidimensional linked matrix factorization for pan-omics pan-cancer analysis, <em> Annals of Applied Statistics </em> 2022.  [(article link)](https://projecteuclid.org/journals/annals-of-applied-statistics/volume-16/issue-1/Bidimensional-linked-matrix-factorization-for-pan-omics-pan-cancer-analysis/10.1214/21-AOAS1495.short)

The functions in EV_Bidifac_Functions.R perform an empirical variational Bayes version of the decompostion, as described in [3]. See EV_Bidifac_Examples.R for example usage. 

[3] Working preprint: https://ericfrazerlock.com/EVBIDI_Manuscript.pdf .  

