# plrd
Partially Linear Regression Discontinuity Inference, as proposed by Ghosh, Imbens and Wager (2025).

The development version of this package can be installed using devtools:

```R
devtools::install_github("ghoshadi/plrd")
```
Replication files for Ghosh, Imbens and Wager (2025) are available in
the directory `Experiments`.

Example usage:

```R
library(plrd)
# Simple example of regression discontinuity design
set.seed(42)
n = 1000; threshold = 0
X = runif(n, -1, 1)
W = as.numeric(X >= threshold)
Y = (1 + 2*W)*(1 + X^2) + 1 / (1 + exp(X)) + rnorm(n, sd = .5)
out = plrd(Y, X, threshold)
print(out)
plot(out)
```

#### References
Aditya Ghosh, Guido Imbens and Stefan Wager.
<b>PLRD : Partially Linear Regression Discontinuity Inference.</b>, [arXiv preprint arXiv:2503.09907](https://arxiv.org/abs/2503.09907).


#### Funding
Development of this software was supported by the National Science Foundation under grant number SES-2242876.
