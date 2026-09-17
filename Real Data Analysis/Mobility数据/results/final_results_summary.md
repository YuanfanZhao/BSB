# Women's Mobility: training/test split results (final)

Protocol: the 8445 observations are split once into a training set (70%) and a
held-out test set (30%); all models are fitted on the training set only, and the
Bernstein degree K is chosen inside the training set by 5-fold x 3 repeated
cross-validation with the validation likelihood, giving K = 25.

## (A) Training and test samples

| X | 0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 |
|---|---|---|---|---|---|---|---|---|---|
| Training frequency | 556 | 1172 | 1622 | 1074 | 805 | 274 | 137 | 103 | 169 |
| Training percentage | 9.40 | 19.82 | 27.44 | 18.17 | 13.62 | 4.63 | 2.32 | 1.74 | 2.86 |
| Test frequency | 273 | 445 | 747 | 459 | 321 | 115 | 56 | 40 | 77 |
| Test percentage | 10.78 | 17.57 | 29.49 | 18.12 | 12.67 | 4.54 | 2.21 | 1.58 | 3.04 |

n = 8445, n1 = 5912 (training), n2 = 2533 (test).

## (B) Parameter estimates and training-set fit

| Model | Parameters | Estimate | #par | logLik | AIC | BIC | Pearson X2 | GOF p |
|---|---|---|---|---|---|---|---|---|
| Binomial | p | 0.3198 | 1 | -12280.73 | 24563.46 | 24570.14 | 45196.83 | 0.001 |
| Beta-Binomial | a, b | 2.4661, 5.1846 | 2 | -11483.10 | 22970.20 | 22983.57 | 1208.49 | 0.001 |
| Logit-Normal-Binomial | mu, sigma | -0.8580, 0.8233 | 2 | -11451.23 | 22906.46 | 22919.83 | 1031.15 | 0.001 |
| Kumaraswamy-Binomial | a, b | 1.8854, 5.9790 | 2 | -11512.76 | 23029.53 | 23042.89 | 1286.11 | 0.001 |
| Bern-Bino (K = 25) | lambda | see below | 26 | -11144.94 | 22341.89 | 22515.69 | 47.84 | 0.001 |

Bern-Bino weights (lambda_0 to lambda_25):

| k | 0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 | 11 | 12 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| lambda_k | 0.0088 | 0 | 0 | 0 | 0 | 0 | 0.1601 | 0.7797 | 0 | 0 | 0 | 0 | 0 |

| k | 13 | 14 | 15 | 16 | 17 | 18 | 19 | 20 | 21 | 22 | 23 | 24 | 25 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| lambda_k | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0.0001 | 0.0335 | 0.0001 | 0.0178 |

## (C) Test-set validation

| Model | #par | test logLik | test AIC | test BIC | test MSE | test L1 |
|---|---|---|---|---|---|---|
| Binomial | 1 | -5285.79 | 10573.57 | 10579.41 | 0.001504 | 0.1289 |
| Beta-Binomial | 2 | -4926.18 | 9856.36 | 9868.03 | 0.001186 | 0.1108 |
| Logit-Normal-Binomial | 2 | -4912.77 | 9829.54 | 9841.21 | 0.001133 | 0.1075 |
| Kumaraswamy-Binomial | 2 | -4937.87 | 9879.73 | 9891.41 | 0.001313 | 0.1158 |
| Bern-Bino (K = 25) | 26 | -4769.26 | 9590.52 | 9742.29 | 0.000413 | 0.0655 |

## Definitions

- Pearson chi-square statistic: X2 = sum_x (O_x - E_x)^2 / E_x, where O_x is the
  observed frequency in the training set and E_x = n1 * pmf(x) is the expected
  frequency under the fitted model. The reported p-value is a Monte Carlo
  (parametric bootstrap) p-value with B = 999 replications, identical to the
  chi-square test used in the simulation study: samples are drawn from the fitted
  distribution and the statistic is recomputed.
- test logLik: sum_x c_x^test * log pmf_train(x); the log-likelihood of the
  held-out test sample evaluated under the model fitted on the training set.
- test AIC = 2k - 2*test_logLik, test BIC = k*log(n2) - 2*test_logLik.
- test MSE = mean_x ( empirical_x - pmf_train(x) )^2.
- test L1 = 0.5 * sum_x | empirical_x - pmf_train(x) |, the total variation
  distance between the empirical test distribution and the fitted distribution.
