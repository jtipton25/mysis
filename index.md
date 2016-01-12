---
layout: index
author: "John Tipton"
date: "January 12, 2016"
---

To determine if there was any difference between the two sized nets, we have a few questions of interest. Our main questions of interest are 1) Is there a difference in density of *Mysids* caught between the two nets, 2) Is there a difference in mean length of *Mysids* caught between the two nets, and 3) Does the density when broken down by age and sex class differ between the two net sizes? Note that we will use count and density interchangeable in what follows because we normalize to a one $m^2$ area.


## 1) Is there a difference in density of *Mysids* caught between the two nets?



```r
load("~/mysis/data/mysisCountData.RData")                                 ## load count data
```


```r
ggplot(data=mysisCountData, aes(count)) + geom_histogram() + facet_grid(net ~ date) + 
  ggtitle("Normalized counts by net size and date")
```

![plot of chunk plotCount](https://github.com/jtipton25/mysis/figure/drafts/2016-01-12-Appendix/plotCount-1.png)

By looking at the distributions of counts across time, we see that they are similar between net sizes, but the normality of the counts is questionable. From the histogram, we see little differences between the two net sizes, except for perhaps a spike in catches of about 120 by the small net in August. Because counts cannot be less that zero and often are right skewed, the assumption of normality in the data is questionable. Hence, methods like linear regression, t-tests, and ANOVA that assume normal distributions are not ideal, but, for completeness we begin with a paired t-test.


```r
y1 <- mysisCountData$count[mysisCountData$net == "Small Net"]
y2 <- mysisCountData$count[mysisCountData$net == "Large Net"]
pairttest <- t.test(y1, y2, paired=TRUE)
pairttest
```

```
## 
## 	Paired t-test
## 
## data:  y1 and y2
## t = 0.07, df = 40, p-value = 0.9
## alternative hypothesis: true difference in means is not equal to 0
## 95 percent confidence interval:
##  -25.3  27.1
## sample estimates:
## mean of the differences 
##                     0.9
```

From the paired t-test, we find that there is no significant difference ($t$ = 0.069, $df$ = 39, $p$-value = 0.945) in expected counts.


```r
mns <- by(mysisCountData[, "count"], mysisCountData[, "net"], mean)       ## mean of counts by net size
vars <- by(mysisCountData[, "count"], mysisCountData[, "net"], var)       ## variance of counts by net size
```

We test for a difference in density between the nets formally by fitting a model to compare if the count of *Mysids* caught by the large net is the same as the count of *Mysids* caught by the large net, after normalizing for the area of the two nets and controlling for covariates. The class of models known as generalized linear models (glms) offers a rigorous solution for modeling count data. These models allow for regression models that account for the effects of covariates like in linear regression but do not assume a normal distribution. For count data, there are two natural distributions to use: the Poisson distribution and the negative binomial distribution. The Poisson distribution has the strong assumption that the mean equals the variance, which is often not met in practical datasets. In our data, the count means are $(183.2, 182.3)$ and the count variances are $(1.524 \times 10^4, 1.795 \times 10^4)$, for the small and large nets, respectively. Therefore, we use a negative binomial model
$$\label{neg:bin} y\_i \sim \mbox{NegBin}(\mu\_i, \phi),$$
where $\mu\_i$ is the mean of the negative binomial distribution for sampling event $i$ and $\phi$ is an overdispersion parameter that allows for the mean and variance to be different. We model $\mu\_i$ with covariates using the canonical log link function
$$ 
\mu\_i = log(\mathbf{X}\_i \boldsymbol{\beta}),
$$
where $\mathbf{X}\_i$ is the set of covariates for observation $i$. Then we perform inference on our coefficients $\boldsymbol{\beta}$, where the interpretation of $\beta\_j$ is a percent change in count per unit change in $X\_{ij}$.


```r
summary(nbmod <- glm.nb(count ~ net + date + station, data=mysisCountData))
```

```
## 
## Call:
## glm.nb(formula = count ~ net + date + station, data = mysisCountData, 
##     init.theta = 6.348525199, link = log)
## 
## Deviance Residuals: 
##     Min       1Q   Median       3Q      Max  
## -2.7717  -0.8767  -0.0976   0.5885   2.1331  
## 
## Coefficients: (1 not defined because of singularities)
##                 Estimate Std. Error z value Pr(>|z|)    
## (Intercept)      4.88574    0.17609   27.75  < 2e-16 ***
## netLarge Net     0.00864    0.09068    0.10  0.92412    
## dateAugust       0.57941    0.33138    1.75  0.08038 .  
## dateSeptember   -0.19079    0.11119   -1.72  0.08617 .  
## station DIL-10  -0.73362    0.40498   -1.81  0.07007 .  
## station DIL-2   -0.44799    0.40363   -1.11  0.26705    
## station DIL-3   -0.53352    0.40400   -1.32  0.18664    
## station DIL-4    0.37937    0.40133    0.95  0.34451    
## station DIL-5   -0.25298    0.40291   -0.63  0.53008    
## station DIL-6   -1.06893    0.40714   -2.63  0.00865 ** 
## station DIL-7   -0.83050    0.40553   -2.05  0.04057 *  
## station DIL-8   -0.59766    0.40429   -1.48  0.13933    
## station DIL-9   -0.47799    0.40376   -1.18  0.23648    
## stationDIL-1     0.62639    0.23359    2.68  0.00733 ** 
## stationDIL-10    0.30033    0.23418    1.28  0.19968    
## stationDIL-2     0.50775    0.23378    2.17  0.02986 *  
## stationDIL-3     0.88135    0.23324    3.78  0.00016 ***
## stationDIL-4     1.24702    0.23288    5.35  8.6e-08 ***
## stationDIL-5     0.11220    0.23462    0.48  0.63249    
## stationDIL-6    -0.21706    0.23561   -0.92  0.35690    
## stationDIL-7    -0.29093    0.23588   -1.23  0.21744    
## stationDIL-8    -0.44658    0.23652   -1.89  0.05901 .  
## stationDIL-9          NA         NA      NA       NA    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## (Dispersion parameter for Negative Binomial(6.35) family taken to be 1)
## 
##     Null deviance: 210.992  on 79  degrees of freedom
## Residual deviance:  82.648  on 58  degrees of freedom
## AIC: 932.7
## 
## Number of Fisher Scoring iterations: 1
## 
## 
##               Theta:  6.35 
##           Std. Err.:  1.03 
## 
##  2 x log-likelihood:  -886.68
```

```r
## get coefficient estimates
coefs <- c(nbmod$coefficients, nbmod$theta)           
```

To test for a change in total counts between the two net sizes, we construct a negative binomial regression model that examines the effect of net size while accounting for date of sampling and sampling location. The results shown below show that there is no evidence of an effect on the number of counts observed between the two net sizes ($z$-value = 0.095, $p$-value=0.924).


## 2) Is there a difference in mean length of *Mysids* caught between the two nets?
First we examine some histograms and get an idea of what the length data look like.


```r
load("~/mysis/data/mysisLengthData.RData")
```


```r
ggplot(data=mysisLengthData, aes(y)) + geom_histogram() + facet_grid(net ~ date) +
  ggtitle("Lengths by net size")
```

![plot of chunk unnamed-chunk-1](https://github.com/jtipton25/mysis/figure/drafts/2016-01-12-Appendix/unnamed-chunk-1-1.png)

When looking at the distribution of lengths between the two net sizes, we see that the histograms look quite similar across the different dates, with a general increasing trend through time. What stands out is that the small net catches about a quarter of the number of *Mysids*, as expected by the difference in net sizes, whereas the general shapes of the histograms appear quite similar across net sizes, suggesting there is not a difference in *Mysid* length distribution between nets across time.


```r
ggplot(data=mysisLengthData, aes(y)) + geom_histogram() + facet_grid(net ~ gender) + 
  ggtitle("Lengths by sex and net size")
```

![plot of chunk unnamed-chunk-2](https://github.com/jtipton25/mysis/figure/drafts/2016-01-12-Appendix/unnamed-chunk-2-1.png)

When we plot the distribution of lengths between the net sizes with respect to sex class, the distributions look quite similar as well. This suggests that there is not a lot of difference in length distribution of *Mysids* between the two net sizes when broken down by sex class.


```r
y1 <- by(mysisLengthData$y[mysisLengthData$net == "Small Net"], 
         mysisLengthData$label[mysisLengthData$net == "Small Net"], mean)
y1 <- y1[!is.na(y1)]

y2 <- by(mysisLengthData$y[mysisLengthData$net == "Large Net"], 
         mysisLengthData$label[mysisLengthData$net == "Large Net"], mean)
y2 <- y2[!is.na(y2)]
pairttestlength <- t.test(y1, y2, paired=TRUE)
pairttestlength
```

```
## 
## 	Paired t-test
## 
## data:  y1 and y2
## t = -2, df = 40, p-value = 0.04
## alternative hypothesis: true difference in means is not equal to 0
## 95 percent confidence interval:
##  -0.5677 -0.0183
## sample estimates:
## mean of the differences 
##                  -0.293
```

```r
qplot(y1-y2, xlab="difference in means")
```

![plot of chunk pairedTTestLength](https://github.com/jtipton25/mysis/figure/drafts/2016-01-12-Appendix/pairedTTestLength-1.png)

From the histogram of paired differences, it looks like there is some light evidence that the smaller net catches slightly smaller *Mysids*, although the distribution appears to be centered near 0. From the paired t-test, we find that there is a significant difference ($t$ = -2.158, $df$ = 39, $p$-value = 0.037) in expected counts, although the $p$-value is close to the 0.05 level. Despite the statistically significant difference, the practical difference is quite small as judged by the histograms and small difference in mean *Mysid* size difference (-0.293 mm).



```r
y1trim <- by(mysisLengthData$y[mysisLengthData$net == "Small Net"],
             mysisLengthData$label[mysisLengthData$net == "Small Net"], mean, trim=0.05)
y1trim <- y1trim[!is.na(y1trim)]

y2trim <- by(mysisLengthData$y[mysisLengthData$net == "Large Net"],
             mysisLengthData$label[mysisLengthData$net == "Large Net"], mean, trim=0.05)
y2trim <- y2trim[!is.na(y2trim)]
pairttestlengthtrim <- t.test(y1trim, y2trim, paired=TRUE)
pairttestlengthtrim
```

```
## 
## 	Paired t-test
## 
## data:  y1trim and y2trim
## t = -2, df = 40, p-value = 0.09
## alternative hypothesis: true difference in means is not equal to 0
## 95 percent confidence interval:
##  -0.5320  0.0379
## sample estimates:
## mean of the differences 
##                  -0.247
```

```r
qplot(y1trim-y2trim, xlab="difference in trimmed means")
```

![plot of chunk pairedTTestLengthtrim](https://github.com/jtipton25/mysis/figure/drafts/2016-01-12-Appendix/pairedTTestLengthtrim-1.png)

One could argue that because the larger net catches four times the number of *Mysids*, the large net is more likely to catch *Mysids* at the extremes of the size classes (either really small or really large *Mysids*). To account for this, we perform a paired t-test using a trimmed mean, removing the smallest and largest 5% of the lengths before calculating the sample mean. From the paired t-test using the trimmed means, we find that there is not a significant difference ($t$ = -1.754, $df$ = 39, $p$-value = 0.087) in expected counts, although the $p$-value is close to the 0.05 level. As with the untrimmed means, the practical difference is quite small as judged by the histograms and small difference in mean *Mysid* size difference (-0.247 mm).

To better account for the different sample sizes (the group means are from unbalanced sample sizes), we fit a linear model and a robust linear model (that assumes overdispersion in the data). We start by checking for heterogeneity in the group variances.


```r
vars <- as.vector(by(mysisLengthData$y, mysisLengthData$label, var))
qplot(vars, main = "var by sampling location", xlab="var")
```

![plot of chunk unnamed-chunk-3](https://github.com/jtipton25/mysis/figure/drafts/2016-01-12-Appendix/unnamed-chunk-3-1.png)

```r
bart_test <- bartlett.test(mysisLengthData$y, mysisLengthData$label)
bart_test
```

```
## 
## 	Bartlett test of homogeneity of variances
## 
## data:  mysisLengthData$y and mysisLengthData$label
## Bartlett's K-squared = 1000, df = 80, p-value <2e-16
```

```r
levene_test <- leveneTest(mysisLengthData$y, mysisLengthData$label)
levene_test
```

```
## Levene's Test for Homogeneity of Variance (center = median)
##         Df F value Pr(>F)    
## group   79    12.3 <2e-16 ***
##       9047                   
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

When looking at the distribution of sampling variance by sampling occasion (10 locations, 2 nets, and 4 different sampling dates gives 80 total sampling occasions), we see that there is some difference in variance of lengths by sampling occasion and Bartlett's test ($p$-value <2e-16) and Levene's test ($p$-value <2e-16).

To test if there is a difference in length of net sizes while accounting for heterogeneity, we use a weighted least squares regression where the weights are the inverse of the group variances.


```r
## weighted linear regression
samps <- by(rep(1, length(mysisLengthData$y)), mysisLengthData$label, sum)
w <- rep(0, length(mysisLengthData$y))
for (i in 1:length(unique(mysisLengthData$label))) {
  if (i==1) {
    w[1:samps[1]] <- rep(1 / vars[i], samps[i])
  } else {
    w[(sum(samps[1:(i-1)])+1):sum(samps[1:i])] <- rep(1 / vars[i], samps[i])
  }
}
summary(weightedlm <- lm(y ~ net + date + gender + station, 
                          data=mysisLengthData, weights=w))
```

```
## 
## Call:
## lm(formula = y ~ net + date + gender + station, data = mysisLengthData, 
##     weights = w)
## 
## Weighted Residuals:
##    Min     1Q Median     3Q    Max 
## -3.594 -0.335 -0.021  0.322  3.888 
## 
## Coefficients:
##               Estimate Std. Error t value Pr(>|t|)    
## (Intercept)    15.2348     0.0835  182.43  < 2e-16 ***
## netLarge Net    0.0858     0.0402    2.13   0.0329 *  
## dateAugust      1.0681     0.0435   24.54  < 2e-16 ***
## dateSeptember   1.2951     0.0451   28.69  < 2e-16 ***
## genderJ        -8.1158     0.0695 -116.74  < 2e-16 ***
## genderM        -1.4670     0.0737  -19.91  < 2e-16 ***
## genderU        -3.9890     0.0705  -56.58  < 2e-16 ***
## stationDIL-10  -0.1835     0.0742   -2.47   0.0135 *  
## stationDIL-2   -0.0170     0.0737   -0.23   0.8173    
## stationDIL-3   -0.5895     0.0512  -11.52  < 2e-16 ***
## stationDIL-4   -0.2424     0.0548   -4.42  9.8e-06 ***
## stationDIL-5    0.2334     0.0856    2.73   0.0064 ** 
## stationDIL-6    0.1715     0.0845    2.03   0.0424 *  
## stationDIL-7    0.5037     0.0981    5.14  2.9e-07 ***
## stationDIL-8    0.0248     0.0890    0.28   0.7807    
## stationDIL-9    0.4405     0.0979    4.50  6.9e-06 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Residual standard error: 0.528 on 9111 degrees of freedom
## Multiple R-squared:  0.823,	Adjusted R-squared:  0.823 
## F-statistic: 2.82e+03 on 15 and 9111 DF,  p-value: <2e-16
```

From the weighted least squares regression, we find that there is a statistically significant difference in mean *Mysid* length caught between the two nets ($t$ = 2.134, $p$-value = 0.033). Now, after accounting for the location, date, sex distribution, and heterogeneity of variance we find that the difference in mean length caught between the two net sizes is 0.086 mm.

As a final analysis, we apply a robust regression, using an M-estimator model that accounts for heterogeneity in variance and presence of outlying observations.


```r
summary(rlmmod <- rlm(y ~ net + date + gender + station, data=mysisLengthData, w=w))
```

```
## 
## Call: rlm(formula = y ~ net + date + gender + station, data = mysisLengthData, 
##     weights = w)
## Residuals:
##      Min       1Q   Median       3Q      Max 
## -3.62219 -0.31110 -0.00929  0.32088  3.90777 
## 
## Coefficients:
##               Value    Std. Error t value 
## (Intercept)     15.316    0.078    196.222
## netLarge Net     0.078    0.038      2.090
## dateAugust       1.191    0.041     29.266
## dateSeptember    1.449    0.042     34.354
## genderJ         -8.244    0.065   -126.880
## genderM         -1.413    0.069    -20.510
## genderU         -4.347    0.066    -65.961
## stationDIL-10   -0.192    0.069     -2.767
## stationDIL-2    -0.071    0.069     -1.032
## stationDIL-3    -0.590    0.048    -12.333
## stationDIL-4    -0.296    0.051     -5.770
## stationDIL-5     0.147    0.080      1.843
## stationDIL-6     0.092    0.079      1.159
## stationDIL-7     0.427    0.092      4.653
## stationDIL-8    -0.017    0.083     -0.208
## stationDIL-9     0.288    0.092      3.149
## 
## Residual standard error: 0.467 on 9111 degrees of freedom
```

```r
tvalues <- summary(rlmmod)$coefficients[2, 3]
df <- summary(rlmmod)$df[2]
pvalues <- pt(abs(tvalues), df, lower.tail=FALSE) * 2
```

For the robust linear model, we still find statistically significant differences in mean length between the net sizes ($t$-value = 2.09, $p$-value = 0.037 on 9111 degrees of freedom), with the effect size 0.078 about the same as from the weighted linear model


```r
SD <- sqrt(((sum(mysisLengthData$net == "Large Net") - 1) *
              sd(mysisLengthData$y[mysisLengthData$net == "Large Net"])^2 + 
              (sum(mysisLengthData$net == "Small Net") - 1) *
              sd(mysisLengthData$y[mysisLengthData$net == "Small Net"])^2) / 
             (sum(mysisLengthData$net == "Large Net") + 
                sum(mysisLengthData$net == "Small Net")))

d <- (mean(mysisLengthData$y[mysisLengthData$net == "Large Net"]) - 
        mean(mysisLengthData$y[mysisLengthData$net == "Small Net"])) / SD
```

Although the statistical tests using both the weighted linear model and the robust linear model show statistically significant differences in mean length caught between the net sizes, this is not unexpected because have very large sample sizes and and the regression coefficients are very small relative to the other sources of variation in the model. Of greater interest is whether the observed effect is of practical significance.

The effect size of the mean length difference between the two net sizes is small (0.078mm) and the sample size is large (n = 9127). Given a large sample size, a hypothesis test will show statistical significance unless the population effect size is exactly zero (this explains why all of the $p$-values in the table above are less than 0.05). Therefore, the practical effect of a difference in mean length of 0.078mm on a species with a mean length of 10.351mm is small (this is the smallest effect of all the effects estimated by almost a factor of two) and a difference in means of this size is not of practical interest. Another measure of effect size is Cohen's $d$ which measures the difference in means relative to a pooled standard deviation. For our data, Cohen's $d$=0.075, which implies that the effect of net size is quite small in terms of practical significance.

## 3)  Does the density when broken down by age and sex class differ between the two net sizes?
The next question we wish to explore is whether the total counts when grouped into sex classes of male, female, juvenile, and unknown vary with net size or other covariates. A visual inspection of counts by sex class shows no difference in counts by net size, but does show a change in counts over time for each sex class.


```r
## load data
load("~/mysis/data/mysisSexCountData.RData")   
## Normalize the counts <- divide by 4 because rather than multiply by 4 to lessen the impact of zeros
mysisSexCountData$count[mysisSexCountData$net == "Large Net"] <- 
  mysisSexCountData$count[mysisSexCountData$net == "Large Net"] / 4

ggplot(data=mysisSexCountData, aes(count)) +         
  geom_histogram() + facet_grid(net ~ gender) + 
  ggtitle("Normalized counts by sex and net size")
```

![plot of chunk unnamed-chunk-6](https://github.com/jtipton25/mysis/figure/drafts/2016-01-12-Appendix/unnamed-chunk-6-1.png)

From the histogram of counts broken down by sex class and net size, we see no outstanding differences of counts between net size.


```r
ggplot(data=mysisSexCountData, aes(count)) +
  geom_histogram() + facet_grid(gender ~ date) + 
  ggtitle("Normalized counts by sex and date")
```

![plot of chunk unnamed-chunk-7](https://github.com/jtipton25/mysis/figure/drafts/2016-01-12-Appendix/unnamed-chunk-7-1.png)

From the histogram of sex class by date, we see a pattern of juveniles maturing to males and females as time progresses and an increase in unknowns as juveniles grow in size but are not sexually differentiated.


```r
## Note: this gives warnings in the model as there are now non-integers
summary(nbmod <- glm.nb(count ~ net + date + station + gender, data = mysisSexCountData))
```

```
## 
## Call:
## glm.nb(formula = count ~ net + date + station + gender, data = mysisSexCountData, 
##     init.theta = 2.056586381, link = log)
## 
## Deviance Residuals: 
##    Min      1Q  Median      3Q     Max  
## -2.712  -0.915  -0.248   0.363   2.933  
## 
## Coefficients:
##               Estimate Std. Error z value Pr(>|z|)    
## (Intercept)     1.3982     0.1738    8.04  8.7e-16 ***
## netLarge Net    0.0159     0.0884    0.18  0.85767    
## dateAugust      0.1145     0.1089    1.05  0.29314    
## dateSeptember   0.2653     0.1079    2.46  0.01391 *  
## station2       -0.3489     0.2119   -1.65  0.09972 .  
## station3       -0.0597     0.1928   -0.31  0.75697    
## station4       -0.0177     0.1925   -0.09  0.92685    
## station5        0.4813     0.1892    2.54  0.01095 *  
## station6       -0.2847     0.1950   -1.46  0.14415    
## station7       -0.6360     0.1887   -3.37  0.00075 ***
## station8       -0.6067     0.1988   -3.05  0.00227 ** 
## station9       -0.7225     0.2005   -3.60  0.00031 ***
## station10      -0.3292     0.1954   -1.68  0.09211 .  
## genderJ         1.8525     0.1283   14.44  < 2e-16 ***
## genderM         0.8199     0.1323    6.20  5.7e-10 ***
## genderU         0.8572     0.1320    6.49  8.5e-11 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## (Dispersion parameter for Negative Binomial(2.06) family taken to be 1)
## 
##     Null deviance: 679.28  on 319  degrees of freedom
## Residual deviance: 345.52  on 304  degrees of freedom
## AIC: 2029
## 
## Number of Fisher Scoring iterations: 1
## 
## 
##               Theta:  2.057 
##           Std. Err.:  0.201 
## 
##  2 x log-likelihood:  -1995.020
```

To correct for the warning in the model above, we work with densities rounded to the nearest whole number (the normalization from counts to densities introduces decimals that need to be rounded). Notice that the results (estimates and $p$-values) do not change in any meaninful manner by rounding.


```r
## Now fit the model with rounded counts
mysisSexCountData$count[mysisSexCountData$net == "Large Net"] <- 
  round(mysisSexCountData$count[mysisSexCountData$net == "Large Net"])
## Note: this no longer gives a warning
summary(nbmod <- glm.nb(count ~ net + date + station + gender, data = mysisSexCountData))
```

```
## 
## Call:
## glm.nb(formula = count ~ net + date + station + gender, data = mysisSexCountData, 
##     init.theta = 2.018481451, link = log)
## 
## Deviance Residuals: 
##    Min      1Q  Median      3Q     Max  
## -3.006  -0.911  -0.250   0.353   2.902  
## 
## Coefficients:
##               Estimate Std. Error z value Pr(>|z|)    
## (Intercept)     1.3862     0.1752    7.91  2.5e-15 ***
## netLarge Net    0.0165     0.0890    0.18  0.85336    
## dateAugust      0.1163     0.1097    1.06  0.28910    
## dateSeptember   0.2696     0.1087    2.48  0.01310 *  
## station2       -0.3543     0.2137   -1.66  0.09738 .  
## station3       -0.0627     0.1945   -0.32  0.74720    
## station4       -0.0125     0.1940   -0.06  0.94882    
## station5        0.4822     0.1908    2.53  0.01148 *  
## station6       -0.2859     0.1966   -1.45  0.14576    
## station7       -0.6230     0.1900   -3.28  0.00104 ** 
## station8       -0.5994     0.2003   -2.99  0.00276 ** 
## station9       -0.7049     0.2018   -3.49  0.00048 ***
## station10      -0.3210     0.1969   -1.63  0.10315    
## genderJ         1.8598     0.1292   14.39  < 2e-16 ***
## genderM         0.8312     0.1332    6.24  4.3e-10 ***
## genderU         0.8618     0.1330    6.48  9.2e-11 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## (Dispersion parameter for Negative Binomial(2.02) family taken to be 1)
## 
##     Null deviance: 674.26  on 319  degrees of freedom
## Residual deviance: 346.12  on 304  degrees of freedom
## AIC: 2032
## 
## Number of Fisher Scoring iterations: 1
## 
## 
##               Theta:  2.018 
##           Std. Err.:  0.198 
## 
##  2 x log-likelihood:  -1997.927
```

To test for a change in counts broken down by sex category between the two net sizes, we construct a negative binomial regression model that examines the effects of net size, date of sampling, sampling location, and sex class. The results show that there is no effect of net size on counts when controlling for sampling date, sampling location, and sex class ($z$-value = 0.185, $p$-value = 0.853). 

We can also examine if the counts within a sex class are different between the net sizes. To begin, we compare the counts of juveniles between the two nets.


```r
juvenileCountData <- data.frame(count= mysisSexCountData$count[mysisSexCountData$gender == "J"],
                                date=mysisSexCountData$date[mysisSexCountData$gender == "J"], 
                                net=mysisSexCountData$net[mysisSexCountData$gender == "J"], 
                                station=mysisSexCountData$station[mysisSexCountData$gender == "J"])
maleCountData <- data.frame(count= mysisSexCountData$count[mysisSexCountData$gender == "M"],
                                date=mysisSexCountData$date[mysisSexCountData$gender == "M"], 
                                net=mysisSexCountData$net[mysisSexCountData$gender == "M"], 
                                station=mysisSexCountData$station[mysisSexCountData$gender == "M"])
femaleCountData <- data.frame(count= mysisSexCountData$count[mysisSexCountData$gender == "F"],
                                date=mysisSexCountData$date[mysisSexCountData$gender == "F"], 
                                net=mysisSexCountData$net[mysisSexCountData$gender == "F"], 
                                station=mysisSexCountData$station[mysisSexCountData$gender == "F"])

ggplot(data=juvenileCountData, aes(count)) +
  geom_histogram() + facet_grid(net ~ date) + 
  ggtitle("Normalized counts of juveniles by net size and date")
```

![plot of chunk unnamed-chunk-10](https://github.com/jtipton25/mysis/figure/drafts/2016-01-12-Appendix/unnamed-chunk-10-1.png)

The histogram shows no visual differences in the juvenile count between net sizes but shows a decrease in juvenile count through time, as expected. To test this formally, we construct a negative binomial model that controls for sampling date and station.


```r
## count model for juveniles
summary(nbmodjuvenile <- glm.nb(count ~ net + date + station, data = juvenileCountData))
```

```
## 
## Call:
## glm.nb(formula = count ~ net + date + station, data = juvenileCountData, 
##     init.theta = 6.242505272, link = log)
## 
## Deviance Residuals: 
##    Min      1Q  Median      3Q     Max  
## -2.374  -0.893  -0.143   0.722   1.943  
## 
## Coefficients:
##               Estimate Std. Error z value Pr(>|z|)    
## (Intercept)     3.8623     0.1741   22.19  < 2e-16 ***
## netLarge Net   -0.0779     0.1066   -0.73   0.4650    
## dateAugust     -0.7571     0.1309   -5.78  7.3e-09 ***
## dateSeptember  -1.3412     0.1413   -9.49  < 2e-16 ***
## station2       -0.3099     0.2505   -1.24   0.2160    
## station3       -0.2378     0.2269   -1.05   0.2945    
## station4        0.3134     0.2211    1.42   0.1564    
## station5        0.8651     0.2176    3.98  7.0e-05 ***
## station6       -0.5967     0.2324   -2.57   0.0102 *  
## station7       -1.1012     0.2272   -4.85  1.3e-06 ***
## station8       -1.2495     0.2479   -5.04  4.6e-07 ***
## station9       -1.1727     0.2456   -4.77  1.8e-06 ***
## station10      -0.8850     0.2383   -3.71   0.0002 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## (Dispersion parameter for Negative Binomial(6.24) family taken to be 1)
## 
##     Null deviance: 371.536  on 79  degrees of freedom
## Residual deviance:  89.702  on 67  degrees of freedom
## AIC: 588.4
## 
## Number of Fisher Scoring iterations: 1
## 
## 
##               Theta:  6.24 
##           Std. Err.:  1.47 
## 
##  2 x log-likelihood:  -560.38
```

```r
## count model for males
summary(nbmodmale <- glm.nb(count ~ net + date + station, data = maleCountData))
```

```
## 
## Call:
## glm.nb(formula = count ~ net + date + station, data = maleCountData, 
##     init.theta = 5.253833744, link = log)
## 
## Deviance Residuals: 
##    Min      1Q  Median      3Q     Max  
## -1.997  -0.889  -0.174   0.558   2.594  
## 
## Coefficients:
##               Estimate Std. Error z value Pr(>|z|)    
## (Intercept)     2.1519     0.2140   10.06   <2e-16 ***
## netLarge Net   -0.0723     0.1266   -0.57    0.568    
## dateAugust      0.0391     0.1577    0.25    0.804    
## dateSeptember   0.3358     0.1522    2.21    0.027 *  
## station2       -0.4230     0.3112   -1.36    0.174    
## station3       -0.0227     0.2746   -0.08    0.934    
## station4       -0.3133     0.2831   -1.11    0.268    
## station5        0.3041     0.2674    1.14    0.256    
## station6        0.1076     0.2715    0.40    0.692    
## station7       -0.4303     0.2719   -1.58    0.114    
## station8       -0.1333     0.2776   -0.48    0.631    
## station9       -0.6654     0.2967   -2.24    0.025 *  
## station10      -0.0364     0.2750   -0.13    0.895    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## (Dispersion parameter for Negative Binomial(5.25) family taken to be 1)
## 
##     Null deviance: 102.047  on 79  degrees of freedom
## Residual deviance:  79.331  on 67  degrees of freedom
## AIC: 477.4
## 
## Number of Fisher Scoring iterations: 1
## 
## 
##               Theta:  5.25 
##           Std. Err.:  1.36 
## 
##  2 x log-likelihood:  -449.43
```

```r
## count model for females
summary(nbmodfemale <- glm.nb(count ~ net + date + station, data = femaleCountData))
```

```
## 
## Call:
## glm.nb(formula = count ~ net + date + station, data = femaleCountData, 
##     init.theta = 3.234583798, link = log)
## 
## Deviance Residuals: 
##    Min      1Q  Median      3Q     Max  
## -2.010  -1.188  -0.187   0.469   2.350  
## 
## Coefficients:
##               Estimate Std. Error z value Pr(>|z|)    
## (Intercept)    -0.0124     0.3686   -0.03   0.9732    
## netLarge Net    0.2885     0.1796    1.61   0.1081    
## dateAugust      0.3763     0.2271    1.66   0.0975 .  
## dateSeptember   1.0463     0.2097    4.99    6e-07 ***
## station2        0.0474     0.5002    0.09   0.9245    
## station3        1.0706     0.4230    2.53   0.0114 *  
## station4       -0.1662     0.4892   -0.34   0.7341    
## station5        1.0967     0.4222    2.60   0.0094 ** 
## station6        0.8771     0.4294    2.04   0.0411 *  
## station7        0.7549     0.4195    1.80   0.0719 .  
## station8        0.7329     0.4349    1.69   0.0919 .  
## station9        0.4528     0.4477    1.01   0.3118    
## station10       1.0367     0.4240    2.44   0.0145 *  
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## (Dispersion parameter for Negative Binomial(3.23) family taken to be 1)
## 
##     Null deviance: 128.787  on 79  degrees of freedom
## Residual deviance:  87.168  on 67  degrees of freedom
## AIC: 371.6
## 
## Number of Fisher Scoring iterations: 1
## 
## 
##               Theta:  3.23 
##           Std. Err.:  1.04 
## 
##  2 x log-likelihood:  -343.58
```

From our model, we see that time of year and location are important in influencing the juvenile counts, but net size is not ($z$-value = -0.731, $p$-value = 0.465). Similar negative binomial regression models for counts of males and females show similar effects ($z$-value = -0.571, $p$-value = 0.568 and $z$-value = 1.607, $p$-value = 0.108 for males and females, respectively). For further details and `R` code of the analysis, see the attached appendix and the online supplement at http://jtipton25.github.io/mysis/.


## Additional extreme value analysis (not in the manuscript) for exploring if the nets sample the extremely large *Mysids* differently

I prefer to look at the histograms and see that the distributions between net sizes for all of the different subsets of the data are consistent rather than do a formal statistical test (which is much more complicated). A quick and dirty extreme value analysis suggests that the model parameters for the generalized extreme value distributions are quite similar for the two net sizes, suggesting that both nets sample the large *Mysids* equally.


```r
small_gev <- gev.fit(xdat=mysisLengthData$y[mysisLengthData$net == "Small Net"])
```

```
## $conv
## [1] 0
## 
## $nllh
## [1] 4951
## 
## $mle
## [1] 8.2704 2.9999 0.0302
## 
## $se
## [1] 0.0875 0.0686 0.0292
```

```r
large_gev <- gev.fit(xdat=mysisLengthData$y[mysisLengthData$net == "Large Net"])
```

```
## $conv
## [1] 0
## 
## $nllh
## [1] 20087
## 
## $mle
## [1] 8.4750 3.1620 0.0271
## 
## $se
## [1] 0.0470 0.0372 0.0155
```
By comparing the estimates for the tail parameter of the extreme value distribution estimate, we find that there is overlap between the two net sizes. The 95% confidence interval for the small net tail parameter is (-0.028,0.089) and the large net tail parameter is (-0.004,0.058). Thus both net sizes do not show evidence of heavy tail distributions.
