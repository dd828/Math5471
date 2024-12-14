This is for the final project for the course MATH 5472. I totally implement 5 algorithms and seal them into a package:"Math5472Project."

Reference Paper: Matrix Completion and Low-Rank SVD via Fast Alternating Least Squares  [link].

Functions: 

rank_restricted_soft_svd: implement the algorithm 2.1 in the reference paper.

rank_restricted_softImpute_ALS_3.1: implement the algorithm 3.1 in the reference paper.

rank_restricted_softImpute_ALS_5.1: implement the algorithm 5.1 in the reference paper.

als_matrix_completion: implement the modified algorithm 5.2(with the addition of a regularization parameter $\lambda$ )  in the reference paper. 

soft_impute: implement the original algorithm of Mazumder et al. (2010), as layed out in steps (2)–(4). I use rank_restricted_soft_svd with $\lambda=0$
to perform svd for $\hat{X}$.


How to reproduce the results in the report?

Step1: Open R studio and install the package: 
```{r}
devtools::install_github("dd828/Math5472_project")
```

Step2: load the package: 
```{r}
library(Math5472Project).
```
