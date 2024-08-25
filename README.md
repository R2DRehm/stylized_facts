# Stylized Facts Analysis of PerkinElmer Inc. (PKI)

## Overview

**Company**: PerkinElmer, Inc.  
**Sector**: Diagnostics, Life Science Research, Food, Environmental, and Industrial Testing  
**Stock Exchange**: NYSE (Ticker: PKI)

## Summary Table

| Statistic              | Daily Log-Returns | Monthly Log-Returns |
|------------------------|------------------:|--------------------:|
| Mean                   | 0.0006442         | 0.01350            |
| Std. Deviation         | 0.01551177        | 0.06322579         |
| Confidence Interval    | [0.001252, 0.0235]| [0.02347, 0.108]    |
| Skewness               | -0.2474           | -0.2267            |
| Kurtosis               | 8.8021            | 3.1209             |
| Excess Kurtosis        | 8.8021            | 3.1209             |
| Min                    | -0.1272930        | -0.17439           |
| Quantile (5%)          | -0.02324          | -0.0998            |
| Quantile (25%)         | -0.0069993        | -0.0237            |
| Median                 | 0.0008088         | 0.01927            |
| Quantile (75%)         | 0.0087134         | 0.05266            |
| Quantile (95%)         | 0.02349           | 0.1086             |
| Max                    | 0.1059092         | 0.18458            |
| Jarque-Bera Statistic  | 3553.4            | 1.1007             |
| Jarque-Bera p-value    | < 2.2e-16         | 0.5768             |
| Lilliefors Test Statistic | 0.06778       | 0.062639           |
| Lilliefors Test p-value   | < 2.2e-16     | 0.2957             |
| Number of Observations | 2515              | 120                |

---

## Stylized Facts Analysis

### 1. **Prices are Non-Stationary**
- The price of PKI has a stochastic trend, showing exponential and irregular growth since 2009.
- Empirical evidence shows different statistics for the periods 2009-2017 (mean: 3.581, std. dev: 0.394) and 2018-2019 (mean: 4.436, std. dev: 0.093).
- The log prices exhibit a random walk with drift. The ACF of log daily prices shows long memory properties.

### 2. **Returns are Stationary**
- The log-returns fluctuate around 0, with statistical tests rejecting the hypothesis of zero mean at a 5% significance level for both daily and monthly returns.

### 3. **Asymmetry**
- Histograms of returns display a leptokurtic distribution with negative skewness, as evidenced by the summary table statistics.

### 4. **Heavy Tails**
- The kurtosis values (8.80 for daily and 3.12 for monthly returns) indicate heavy tails.
- QQ plots confirm the heavy tails, showing deviations in the tails compared to a normal distribution.

### 5. **High Frequency Non-Gaussianity**
- Jarque-Bera test strongly rejects normality for daily returns.
- The Lilliefors test suggests normality for monthly returns with a p-value of 0.2957.
- Lower frequency returns tend to be closer to a normal distribution.

### 6. **Returns are Not Autocorrelated**
- The ACF of daily and monthly returns are close to zero. The Box-Pierce and Ljung-Box tests do not reject the null hypothesis of zero autocorrelation.

### 7. **Volatility Clustering and Long-Range Dependence**
- Volatility clustering is evident in the squared log-returns, showing periods of high and low volatility.
- The ACF for absolute and squared returns decays slowly, suggesting long-memory properties.

### 8. **Leverage Effect**
- A negative correlation is observed between returns and future volatility, evidenced by the correlation with the VIX index.

---

## ARCH & GARCH Analysis

### Introduction
ARCH and GARCH models are used to forecast volatility:
- **ARCH**: AutoRegressive Conditional Heteroskedasticity
- **GARCH**: Generalized AutoRegressive Conditional Heteroskedasticity

### Model Estimations
- An **ARCH(4)** model excluding AR coefficients is better suited, with all parameters being significant.
- GARCH(1,1) model fits well, with simpler models like IGARCH and EGARCH also providing good estimates.

### Simulation and Stylized Facts
- Simulations of ARCH(1) and GARCH processes exhibit similar properties to real asset returns:
  - **Non-Stationary Prices, Stationary Returns** (Stylized Fact 1)
  - **Insignificant Autocorrelations** (Stylized Fact 2)
  - **Heavy Tails** (Stylized Fact 3)
  - **Volatility Clustering** (Stylized Fact 5)

---

**Note**: For additional statistical details and the R code used in this analysis, please refer to the appendix.
