---
title: "Effect of net size on estimates of abundance, size, age and sex ratio of Mysis diluviana Journal of Great Lakes Research - Statistical Analysis"
author: "John Tipton"
date: "January 3, 2016"
output: html_document
---


```r
## We start by loading packages and define helper functions
options(digits=3)
set.seed(123)
library(knitr)
# knitr::opts_chunk$set(fig.path='../images/mysisAnalysis-')
knitr::opts_chunk$set(echo=TRUE)
knitr::opts_chunk$set(eval=TRUE)
knitr::opts_chunk$set(tidy=TRUE)
knitr::opts_chunk$set(include=TRUE)
knitr::opts_chunk$set(message=FALSE)
knitr::opts_chunk$set(warning=FALSE)
library(ggplot2)
```

```
## Warning: package 'ggplot2' was built under R version 3.2.3
```

```r
library(Gmisc)
library(MASS)
library(stringr)
library(mgcv)
library(lme4)
library(ismev)
library(car)
```

```
## Warning: package 'car' was built under R version 3.2.3
```

```r
library(betareg)
```

To determine if there was any difference between the two sized nets, we have a few questions of interest. Our main questions of interest are 1) Is there a difference in density of *Mysids* caught between the two nets, 2) Is there a difference in mean length of *Mysids* caught between the two nets, and 3) Does the density when broken down by age and sex class differ between the two net sizes? Note that we will use count and density interchangeable in what follows because we normalize to a one $m^2$ area.


1) Is there a difference in density of *Mysids* caught between the two nets?



```r
load("~/mysis/data/mysisCountData.RData")  ## load count data
```


```r
ggplot(data = mysisCountData, aes(count)) + geom_histogram() + facet_grid(net ~ 
    date) + ggtitle("Normalized counts by net size and date")
```

![plot of chunk plotCount](/figure/drafts/mysisAnalysis/plotCount-1.png)

By looking at the distributions of counts across time, we see that they are similar between net sizes, but the normality of the counts is questionable. From the histogram, we see little differences between the two net sizes, except for perhaps a spike in catches of about 120 by the small net in August. Because counts cannot be less that zero and often are right skewed, the assumption of normality in the data is questionable. Hence, methods like linear regression, t-tests, and ANOVA that assume normal distributions are not ideal, but, for completeness we begin with a paired t-test.


```r
y1 <- mysisCountData$count[mysisCountData$net == "Small Net"]
y2 <- mysisCountData$count[mysisCountData$net == "Large Net"]
pairttest <- t.test(y1, y2, paired = TRUE)
```

From the paired t-test, we find that there is no significant difference ($t$ = 0.069, $df$ = 39, $p$-value = 0.945) in expected counts.


```r
mns <- by(mysisCountData[, "count"], mysisCountData[, "net"], mean)  ## mean of counts by net size
vars <- by(mysisCountData[, "count"], mysisCountData[, "net"], var)  ## variance of counts by net size
```

We test for a difference in density between the nets formally by fitting a model to compare if the count of *Mysids* caught by the large net is the same as the count of *Mysids* caught by the large net, after normalizing for the area of the two nets and controlling for covariates. The class of models known as generalized linear models (glms) offers a rigorous solution for modeling count data. These models allow for regression models that account for the effects of covariates like in linear regression but do not assume a normal distribution. For count data, there are two natural distributions to use: the Poisson distribution and the negative binomial distribution. The Poisson distribution has the strong assumption that the mean equals the variance, which is often not met in practical datasets. In our data, the count means are $(183.2, 182.3)$ and the count variances are $(1.524 &times; 10<sup>4</sup>, 1.795 &times; 10<sup>4</sup>)$, for the small and large nets, respectively. Therefore, we use a negative binomial model
$$
\label{neg:bin}
y_i \sim \mbox{NegBin}(\mu_i, \phi),
$$
where $\mu_i$ is the mean of the negative binomial distribution for sampling event $i$ and $\phi$ is an overdispersion parameter that allows for the mean and variance to be different. We model $\mu_i$ with covariates using the canonical log link function
$$ 
\mu_i = log(\mathbf{X}_i \boldsymbol{\beta}),
$$
where $\mathbf{X}_i$ is the set of covariates for observation $i$. Then we perform inference on our coefficients $\boldsymbol{\beta}$, where the interpretation of $\beta_j$ is a percent change in count per unit change in $X_{ij}$.


```r
summary(nbmod <- glm.nb(count ~ net + date + station, data = mysisCountData))
## get coefficient estimates
coefs <- c(nbmod$coefficients, nbmod$theta)
```

To test for a change in total counts between the two net sizes, we construct a negative binomial regression model that examines the effect of net size while accounting for date of sampling and sampling location. The results shown below show that there is no evidence of an effect on the number of counts observed between the two net sizes ($z$-value = 0.095, $p$-value=0.924).


```r
tableData <- cbind(summary(nbmod)$coeff[1:6, ])
colnames(tableData) <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
rownames(tableData) <- c("Intercept", "Net Size", "August", "September", "Net Size * August", 
    "Net Size * September")
kable(tableData, digits = 4, format = "html")
```



2) Is there a difference in mean length of *Mysids* caught between the two nets?
First we examine some histograms and get an idea of what the length data look like.


```r
load("~/mysis/data/mysisLengthData.RData")
```


```r
ggplot(data = mysisLengthData, aes(y)) + geom_histogram() + facet_grid(net ~ 
    date) + ggtitle("Lengths by net size")
```

![plot of chunk unnamed-chunk-1](/figure/drafts/mysisAnalysis/unnamed-chunk-1-1.png)

When looking at the distribution of lengths between the two net sizes, we see that the histograms look quite similar across the different dates, with a general increasing trend through time. What stands out is that the small net catches about a quarter of the number of *Mysids*, as expected by the difference in net sizes, whereas the general shapes of the histograms appear quite similar across net sizes, suggesting there is not a difference in *Mysid* length distribution between nets across time.


```r
ggplot(data = mysisLengthData, aes(y)) + geom_histogram() + facet_grid(net ~ 
    gender) + ggtitle("Lengths by sex and net size")
```

![plot of chunk unnamed-chunk-2](/figure/drafts/mysisAnalysis/unnamed-chunk-2-1.png)

When we plot the distribution of lengths between the net sizes with respect to sex class, the distributions look quite similar as well. This suggests that there is not a lot of difference in length distribution of *Mysids* between the two net sizes when broken down by sex class.


```r
y1 <- by(mysisLengthData$y[mysisLengthData$net == "Small Net"], mysisLengthData$label[mysisLengthData$net == 
    "Small Net"], mean)
y1 <- y1[!is.na(y1)]

y2 <- by(mysisLengthData$y[mysisLengthData$net == "Large Net"], mysisLengthData$label[mysisLengthData$net == 
    "Large Net"], mean)
y2 <- y2[!is.na(y2)]
pairttestlength <- t.test(y1, y2, paired = TRUE)
qplot(y1 - y2, xlab = "difference in means")
```

![plot of chunk pairedTTestLength](/figure/drafts/mysisAnalysis/pairedTTestLength-1.png)

From the histogram of paired differences, it looks like there is some light evidence that the smaller net catches slightly smaller *Mysids*, although the distribution appears to be centered near 0. From the paired t-test, we find that there is a significant difference ($t$ = -2.158, $df$ = 39, $p$-value = 0.037) in expected counts, although the $p$-value is close to the 0.05 level. Despite the statistically significant difference, the practical difference is quite small as judged by the histograms and small difference in mean *Mysid* size difference (-0.293 mm).



```r
y1trim <- by(mysisLengthData$y[mysisLengthData$net == "Small Net"], mysisLengthData$label[mysisLengthData$net == 
    "Small Net"], mean, trim = 0.05)
y1trim <- y1trim[!is.na(y1trim)]

y2trim <- by(mysisLengthData$y[mysisLengthData$net == "Large Net"], mysisLengthData$label[mysisLengthData$net == 
    "Large Net"], mean, trim = 0.05)
y2trim <- y2trim[!is.na(y2trim)]
pairttestlengthtrim <- t.test(y1trim, y2trim, paired = TRUE)
qplot(y1trim - y2trim, xlab = "difference in trimmed means")
```

![plot of chunk pairedTTestLengthtrim](/figure/drafts/mysisAnalysis/pairedTTestLengthtrim-1.png)

One could argue that because the larger net catches four times the number of *Mysids*, the large net is more likely to catch *Mysids* at the extremes of the size classes (either really small or really large *Mysids*). To account for this, we perform a paired t-test using a trimmed mean, removing the smallest and largest 5% of the lengths before calculating the sample mean. From the paired t-test using the trimmed means, we find that there is not a significant difference ($t$ = -1.754, $df$ = 39, $p$-value = 0.087) in expected counts, although the $p$-value is close to the 0.05 level. As with the untrimmed means, the practical difference is quite small as judged by the histograms and small difference in mean *Mysid* size difference (-0.247 mm).

To better account for the different sample sizes (the group means are from unbalanced sample sizes), we fit a linear model and a robust linear model (that assumes overdispersion in the data). We start by checking for heterogeneity in the group variances.


```r
vars <- as.vector(by(mysisLengthData$y, mysisLengthData$label, var))
qplot(vars, main = "var by sampling location", xlab = "var")
```

![plot of chunk unnamed-chunk-3](/figure/drafts/mysisAnalysis/unnamed-chunk-3-1.png)

```r
bart_test <- bartlett.test(mysisLengthData$y, mysisLengthData$label)
levene_test <- leveneTest(mysisLengthData$y, mysisLengthData$label)
```

When looking at the distribution of sampling variance by sampling occasion (10 locations, 2 nets, and 4 different sampling dates gives 80 total sampling occasions), we see that there is some difference in variance of lengths by sampling occasion and Bartlett's test ($p$-value = 7.43 &times; 10<sup>-175</sup>) and Levene's test ($p$-values = 5.361 &times; 10<sup>-145</sup>).

To test if there is a difference in length of net sizes while accounting for heterogeneity, we use a weighted least squares regression where the weights are the inverse of the group variances.


```r
## weighted linear regression
samps <- by(rep(1, length(mysisLengthData$y)), mysisLengthData$label, sum)
w <- rep(0, length(mysisLengthData$y))
for (i in 1:length(unique(mysisLengthData$label))) {
    if (i == 1) {
        w[1:samps[1]] <- rep(1/vars[i], samps[i])
    } else {
        w[(sum(samps[1:(i - 1)]) + 1):sum(samps[1:i])] <- rep(1/vars[i], samps[i])
    }
}
summary(weightedlm <- lm(y ~ net + date + gender + station, data = mysisLengthData, 
    weights = w))
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
summary(lmmod <- lm(y ~ net + date + gender + station, data = mysisLengthData))
```

```
## 
## Call:
## lm(formula = y ~ net + date + gender + station, data = mysisLengthData)
## 
## Residuals:
##    Min     1Q Median     3Q    Max 
## -15.14  -1.08  -0.06   1.02  10.84 
## 
## Coefficients:
##               Estimate Std. Error t value Pr(>|t|)    
## (Intercept)    15.7599     0.0877  179.68  < 2e-16 ***
## netLarge Net    0.1188     0.0428    2.77   0.0056 ** 
## dateAugust      0.9915     0.0458   21.64  < 2e-16 ***
## dateSeptember   1.0257     0.0484   21.21  < 2e-16 ***
## genderJ        -8.8020     0.0670 -131.40  < 2e-16 ***
## genderM        -1.5957     0.0719  -22.19  < 2e-16 ***
## genderU        -4.2663     0.0710  -60.08  < 2e-16 ***
## stationDIL-10  -0.1645     0.0824   -2.00   0.0459 *  
## stationDIL-2   -0.0105     0.0713   -0.15   0.8833    
## stationDIL-3   -0.5019     0.0674   -7.44  1.1e-13 ***
## stationDIL-4   -0.1939     0.0600   -3.23   0.0012 ** 
## stationDIL-5    0.1631     0.0790    2.06   0.0391 *  
## stationDIL-6    0.0950     0.0808    1.18   0.2398    
## stationDIL-7    0.3619     0.0870    4.16  3.2e-05 ***
## stationDIL-8   -0.0300     0.0900   -0.33   0.7389    
## stationDIL-9    0.2367     0.0825    2.87   0.0041 ** 
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Residual standard error: 1.63 on 9111 degrees of freedom
## Multiple R-squared:  0.83,	Adjusted R-squared:  0.83 
## F-statistic: 2.96e+03 on 15 and 9111 DF,  p-value: <2e-16
```

```r
summary(rlmmod <- rlm(y ~ net + date + gender + station, data = mysisLengthData, 
    w = w))
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
pvalues <- pt(abs(tvalues), df, lower.tail = FALSE) * 2
```


```r
tableData <- cbind(summary(rlmmod)$coeff, pvalues)
```

For the robust linear model, we still find statistically significant differences in mean length between the net sizes ($t$-value = 2.09, $p$-value = 0.037 on 9111 degrees of freedom), with the effect size 0.078 about the same as from the weighted linear model


```r
SD <- sqrt(((sum(mysisLengthData$net == "Large Net") - 1) * sd(mysisLengthData$y[mysisLengthData$net == 
    "Large Net"])^2 + (sum(mysisLengthData$net == "Small Net") - 1) * sd(mysisLengthData$y[mysisLengthData$net == 
    "Small Net"])^2)/(sum(mysisLengthData$net == "Large Net") + sum(mysisLengthData$net == 
    "Small Net")))

d <- (mean(mysisLengthData$y[mysisLengthData$net == "Large Net"]) - mean(mysisLengthData$y[mysisLengthData$net == 
    "Small Net"]))/SD
```

Although the statistical tests using both the weighted linear model and the robust linear model show statistically significant differences in mean length caught between the net sizes, this is not unexpected because have very large sample sizes and and the regression coefficients are very small relative to the other sources of variation in the model. Of greater interest is whether the observed effect is of practical significance.

The effect size of the mean length difference between the two net sizes is small (0.078mm) and the sample size is large (n = 9127). Given a large sample size, a hypothesis test will show statistical significance unless the population effect size is exactly zero (this explains why all of the p-values in the table above are less than 0.05). Therefore, the practical effect of a difference in mean length of 0.078mm on a species with a mean length of 10.351mm is small (this is the smallest effect of all the effects estimated by almost a factor of two) and a difference in means of this size is not of practical interest. Another measure of effect size is Cohen's $d$ which measures the difference in means relative to a pooled standard deviation. For our data, Cohen's $d$=0.075, which implies that the effect of net size is quite small in terms of practical significance.

3)  Does the density when broken down by age and sex class differ between the two net sizes?
The next question we wish to explore is whether the total counts when grouped into sex classes of male, female, juvenile, and unknown vary with net size or other covariates. A visual inspection of counts by sex class shows no difference in counts by net size, but does show a change in counts over time for each sex class.


```r
## load data
load("~/mysis/data/mysisSexCountData.RData")
## Normalize the counts <- divide by 4 because of the zeros for the small net
mysisSexCountData$count[mysisSexCountData$net == "Large Net"] <- mysisSexCountData$count[mysisSexCountData$net == 
    "Large Net"]/4

ggplot(data = mysisSexCountData, aes(count)) + geom_histogram() + facet_grid(net ~ 
    gender) + ggtitle("Normalized counts by sex and net size")
```

![plot of chunk unnamed-chunk-6](/figure/drafts/mysisAnalysis/unnamed-chunk-6-1.png)

From the histogram of counts broken down by sex class and net size, we see no outstanding differences of counts between net size.


```r
ggplot(data = mysisSexCountData, aes(count)) + geom_histogram() + facet_grid(gender ~ 
    date) + ggtitle("Normalized counts by sex and date")
```

![plot of chunk unnamed-chunk-7](/figure/drafts/mysisAnalysis/unnamed-chunk-7-1.png)

From the histogram of sex class by date, we see a pattern of juveniles maturing to males and females as time progresses and an increase in unknowns as juveniles grow in size but are not sexually differentiated.


```r
## Note: this gives warnings in the model as there are now non-integers
summary(nbmod <- glm.nb(count ~ net + date + station + gender, data = mysisSexCountData))
## Now fit the model with rounded counts
mysisSexCountData$count[mysisSexCountData$net == "Large Net"] <- round(mysisSexCountData$count[mysisSexCountData$net == 
    "Large Net"])
## Note: this no longer gives a warning
summary(nbmod <- glm.nb(count ~ net + date + station + gender, data = mysisSexCountData))
```

To test for a change in counts broken down by sex category between the two net sizes, we construct a negative binomial regression model that examines the effects of net size, date of sampling, sampling location, and sex class. The results show that there is no effect of net size on counts when controlling for sampling date, sampling location, and sex class ($z$-value = 0.185, $p$-value = 0.853). 

We can also examine if the counts within a sex class are different between the net sizes. To begin, we compare the counts of juveniles between the two nets.


```r
juvenileCountData <- data.frame(count = mysisSexCountData$count[mysisSexCountData$gender == 
    "J"], date = mysisSexCountData$date[mysisSexCountData$gender == "J"], net = mysisSexCountData$net[mysisSexCountData$gender == 
    "J"], station = mysisSexCountData$station[mysisSexCountData$gender == "J"])
maleCountData <- data.frame(count = mysisSexCountData$count[mysisSexCountData$gender == 
    "M"], date = mysisSexCountData$date[mysisSexCountData$gender == "M"], net = mysisSexCountData$net[mysisSexCountData$gender == 
    "M"], station = mysisSexCountData$station[mysisSexCountData$gender == "M"])
femaleCountData <- data.frame(count = mysisSexCountData$count[mysisSexCountData$gender == 
    "F"], date = mysisSexCountData$date[mysisSexCountData$gender == "F"], net = mysisSexCountData$net[mysisSexCountData$gender == 
    "F"], station = mysisSexCountData$station[mysisSexCountData$gender == "F"])

ggplot(data = juvenileCountData, aes(count)) + geom_histogram() + facet_grid(net ~ 
    date) + ggtitle("Normalized counts of juveniles by net size and date")
```

![plot of chunk unnamed-chunk-9](/figure/drafts/mysisAnalysis/unnamed-chunk-9-1.png)

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

