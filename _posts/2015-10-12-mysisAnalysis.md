---
title: "Mysis analysis"
author: "John Tipton"
date: "September 9, 2015"
output: html_document
---
This document shows the analyses performed in Johnson et. al. We start by loading packages and define helper functions



## Model for Total counts for each net
We begin the analysis by loading the data used for measuring total count. At each sampling event $i$ we recorded the number of mysis shrimp caught and important covariates including date, site, net size, etc. 

```r
load("~/mysis/data/mysisCountData.RData")                                 ## load count data
mns <- by(mysisCountData[, "count"], mysisCountData[, "net"], mean)       ## mean of counts by net size
vars <- by(mysisCountData[, "count"], mysisCountData[, "net"], var)       ## variance of counts by net size
```

After loading the data, we can look at histograms of data to get a better understanding of what variables effect mysis shrimp count. First we examine the total counts broken down by net size and sampling month


```r
ggplot(data=mysisCountData, aes(count)) + geom_histogram() + facet_grid(net ~ date) + 
  ggtitle("Histogram of normalized counts by net size and date")
```

<img src="../images/mysisAnalysis-plotCount-1.png" title="plot of chunk plotCount" alt="plot of chunk plotCount" width="300px" style="display: block; margin: auto;" />

From these histograms, we see little differences between the two net sizes, except for perhaps a spike in catches of about 120 by the small net in August. We test this formally by fitting a model to compare if the count of shrimp caught by the large net is the same as the count of shrimp caught by the large net, after normalizing for the area of the two nets and controlling for covariates. Because counts cannot be less that zero and often are right skewed, the assumption of Normality in the data is questionable. Hence, we don't want to use Gaussian methods like linear regression and ANOVA. Luckily, the class of models known a generalized linear models (or glms) offers a solution. These models allow for regression models that control for the effects of covariates like linear regression but do not assume the normal distribution. For count data, there are two natural distributions to use: the Poisson distribution and the negative binomial distribution. The Poisson distribution has the strong assumption that the mean equals the variance, which is often not met in practical datasets. In our data, the count means are $(182.3, 183.2)$ and the count variances are $(1.795 &times; 10<sup>4</sup>, 1.524 &times; 10<sup>4</sup>)$, for the small and large nets, respectively. Therefore we use a negative binomial model
$$
\label{neg:bin}
y_i \sim \mbox{NegBin}(\mu_i, \phi),
$$
where $\mu$ is the mean of the negative binomial distribution and $\phi$ is an overdispersion parameter that allows for the mean and variance to be different. We model $\mu$ with covariates using a log link function
$$ 
\mu_i = log(\mathbf{X}_i \boldsymbol{\beta}),
$$
where $\mathbf{X}_i$ is the set of covariates for observation $i$. Then we perform inference on our coefficients $\boldsymbol{\beta}$, where the interpretation of $\beta_j$ is a percent change in count per unit change in $X_{ij}$.

```r
## construct covariates of net size, date and their interaction
X <- model.matrix(~ net * date - 1, data=mysisCountData)  
## fit negative binomial model
summary(nbmod <- glm.nb(y ~ ., data = data.frame(y=mysisCountData$count, X)))
## get coefficient estimates
coefs <- c(nbmod$coefficients, nbmod$theta)           
```

To test for a change in total counts between the two net sizes, we construct a negative binomal regression model that examines the effects of net size, date of sampling and the interaction between net size and sampling date. The results shown below show that there is no evidence of an effect on the number of counts observed between the two net sizes (p=0.955).



```r
tableData <- cbind(summary(nbmod)$coeff[1:6, ])
colnames(tableData) <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
rownames(tableData) <- c("Intercept", "Net Size", "August", "September", 
                         "Net Size * August", "Net Size * September")
kable(tableData, digits = 4, format="pandoc")
```

                        Estimate   Std. Error   z value   Pr(>|z|)
---------------------  ---------  -----------  --------  ---------
Intercept                 5.3132        0.139   38.1112      0.000
Net Size                 -0.0111        0.197   -0.0565      0.955
August                   -0.3420        0.242   -1.4138      0.157
September                -0.0812        0.242   -0.3360      0.737
Net Size * August         0.2459        0.342    0.7193      0.472
Net Size * September     -0.2734        0.342   -0.7997      0.424

## Model for counts by sex

The next question we wish to explore is whether the total counts when grouped into sex classes of male, female, juvenile, and unknown vary with net size or other covariates. To do this we use (\ref{eq:negBin}) above. We begin by loading the data where counts are divided into sex classes.


```r
## load data
load("~/mysis/data/mysisSexCountData.RData")   
## construct model matrix for covariates
X <- model.matrix(~ net * gender + net * date + gender * date-1, 
                  data=mysisSexCountData)
## fit model
nbmod <- glm.nb(y ~ ., data = data.frame(y=mysisSexCountData$count, X))
```

To explore this data in more detail we examine two sets of histograms, the counts of each sex class broken down by date in figure1 and the counts of each sex class broken down by net size in figure2. From these we see...


```r
ggplot(data=mysisSexCountData, aes(count)) +         
  geom_histogram() + facet_grid(net ~ gender) + 
  ggtitle("Histogram of normalized counts by sex and net size")
ggplot(data=mysisSexCountData, aes(count)) +
  geom_histogram() + facet_grid(gender ~ date) + 
  ggtitle("Histogram of normalized counts by sex and date")
```

<img src="../images/mysisAnalysis-plot SexCount-1.png" title="plot of chunk plot SexCount" alt="plot of chunk plot SexCount" width="300px" style="display: block; margin: auto;" /><img src="../images/mysisAnalysis-plot SexCount-2.png" title="plot of chunk plot SexCount" alt="plot of chunk plot SexCount" width="300px" style="display: block; margin: auto;" />

To test for a change in counts broken down by gender category between the two net sizes, we construct a negative binomal regression model that examines the effects of net size, date of sampling, the interaction between net size and sampling date, and the interaction between sampling date and gender. The results shown below show that there is no evidence of an interaction between net size and sampling date and no evidence of an interaction between net size and gender. 


```r
tableData <- cbind(summary(nbmod)$coeff)
colnames(tableData) <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
rownames(tableData) <- c("Intercept", "Net Size", "Females", "Juveniles", 
                         "Unknown Gender", "August", "September", 
                         "Net Size * Females", "Net Size * Juveniles", 
                         "Net Size * Unknowns", "Net Size * August",
                         "Net Size * September",
                         "Females * August", "Juveniles * August", 
                         "Unknown Gender * August",
                         "Females * September", "Juveniles * September", 
                         "Unknown Gender * September")
kable(tableData, digits = 4, format="pandoc")
```

                              Estimate   Std. Error   z value   Pr(>|z|)
---------------------------  ---------  -----------  --------  ---------
Intercept                       3.3587        0.158    21.285     0.0000
Net Size                        0.0678        0.196     0.347     0.7287
Females                        -1.2668        0.220    -5.755     0.0000
Juveniles                       1.5758        0.212     7.436     0.0000
Unknown Gender                 -0.5729        0.216    -2.657     0.0079
August                          0.1567        0.239     0.655     0.5127
September                       0.1468        0.239     0.615     0.5387
Net Size * Females              0.3123        0.252     1.240     0.2151
Net Size * Juveniles            0.0037        0.245     0.015     0.9880
Net Size * Unknowns             0.0879        0.247     0.355     0.7224
Net Size * August              -0.3157        0.215    -1.467     0.1425
Net Size * September            0.2531        0.214     1.180     0.2380
Females * August                0.3251        0.309     1.052     0.2929
Juveniles * August             -0.7032        0.300    -2.342     0.0192
Unknown Gender * August         1.0094        0.303     3.332     0.0009
Females * September             0.5445        0.306     1.779     0.0753
Juveniles * September          -1.5523        0.300    -5.166     0.0000
Unknown Gender * September      0.9672        0.302     3.204     0.0014

Thus, we refit the model without the interaction between net size and sampling date. To do this we start by loading the new covariate matrix $\mathbf{X}$.



On this reduced model, the results shown below show that there is no evidence of an effect on the number of counts observed between the two net sizes (p=0.636).

```r
tableData <- cbind(summary(nbmod)$coeff)
colnames(tableData) <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
rownames(tableData) <- c("Intercept", "Net Size", "Females", "Juveniles", 
                         "Unknown Gender", "August", "September", 
                         "Net Size * Females", "Net Size * Juveniles", 
                         "Net Size * Unknowns", 
                         "Females * August", "Juveniles * August", 
                         "Unknown Gender * August",
                         "Females * September", "Juveniles * September", 
                         "Unknown Gender * September")
kable(tableData, digits = 4, format="pandoc")
```

                              Estimate   Std. Error   z value   Pr(>|z|)
---------------------------  ---------  -----------  --------  ---------
Intercept                       3.4341        0.152   22.5222     0.0000
Net Size                       -0.0834        0.176   -0.4739     0.6356
Females                        -1.2681        0.222   -5.7202     0.0000
Juveniles                       1.5772        0.214    7.3834     0.0000
Unknown Gender                 -0.5726        0.217   -2.6358     0.0084
August                          0.0092        0.216    0.0425     0.9661
September                       0.2911        0.215    1.3543     0.1756
Net Size * Females              0.3160        0.254    1.2454     0.2130
Net Size * Juveniles            0.0019        0.247    0.0079     0.9937
Net Size * Unknowns             0.0877        0.249    0.3516     0.7252
Females * August                0.3390        0.311    1.0891     0.2761
Juveniles * August             -0.7198        0.303   -2.3784     0.0174
Unknown Gender * August         1.0143        0.305    3.3231     0.0009
Females * September             0.5426        0.308    1.7596     0.0785
Juveniles * September          -1.5839        0.303   -5.2294     0.0000
Unknown Gender * September      0.9655        0.304    3.1741     0.0015

## Model for proportion of juveniles














