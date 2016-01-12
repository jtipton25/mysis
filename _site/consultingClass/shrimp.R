shrimp <- read.csv("~/consultingClass/Net comparison dataset2.csv")
y_data <- subset(shrimp, !duplicated(shrimp$Sample.Label))
hist(y_data$Count)
y_data$Count[y_data$Net.Size..m. == 0.5] <- 4 * y_data$Count[y_data$Net.Size..m. == 0.5]
hist(y_data$Count)

require(ggplot2)
require(sandwich)
require(msm)
require(MASS)
require(reshape2)

# summary(lin_mod <- lm(y_data$Count ~ y_data$Net.Size..m + y_data$Depth..m. + y_data$Station + y_data$Collection.Date))
# plot(lin_mod)

summary(lin_mod <- lm(y_data$Count ~ y_data$Net.Size..m))
plot(lin_mod)

# summary(pos_mod <- glm(y_data$Count ~ y_data$Net.Size..m + y_data$Depth..m. + y_data$Station + y_data$Collection.Date, family="poisson"))
# plot(pos_mod)

summary(pos_mod <- glm(y_data$Count ~ y_data$Net.Size..m, family="poisson"))
plot(pos_mod)

summary(pos_mod <- glm(y_data$Count ~ y_data$Net.Size..m * y_data$Depth..m. , family="poisson"))

# summary(qpos_mod <- glm(y_data$Count ~ y_data$Net.Size..m + y_data$Depth..m. + y_data$Station + y_data$Collection.Date, family="quasipoisson"))
# plot(qpos_mod)

summary(qpos_mod <- glm(y_data$Count ~ y_data$Net.Size..m, family="quasipoisson"))
layout(matrix(1:4, 2, 2))
plot(qpos_mod)

summary(nb_mod <- glm.nb(y_data$Count ~ y_data$Net.Size..m + y_data$Depth..m. + y_data$Station + y_data$Collection.Date))
layout(matrix(1:4, 2, 2))
plot(nb_mod)

summary(nb_mod1 <- glm.nb(y_data$Count ~ y_data$Net.Size..m + y_data$Depth..m. + y_data$Station))
layout(matrix(1:4, 2, 2))
plot(nb_mod1)

summary(nb_mod2 <- glm.nb(y_data$Count ~ as.factor(y_data$Net.Size..m)))
layout(matrix(1:4, 2, 2))
plot(nb_mod2)

summary(nb_mod3 <- glm.nb(y_data$Count ~ as.factor(y_data$Net.Size..m) + y_data$Depth..m.))
plot(nb_mod3)



pchisq(2 * (logLik(nb_mod1) - logLik(nb_mod2)), df = 1, lower.tail = FALSE)



library(glmnet)

glmnet(y_data$Count ~ as.factor(y_data$Net.Size..m), family = 'poisson')


y_data_frame = as.data.frame(y_data)

g <- ggplot(y_data_frame, aes(x = Count, fill = Gender)) + geom_density(alpha = 0.2)
g

h <- ggplot(y_data_frame, aes(x = Count, fill = as.factor(Net.Size..m.))) + geom_density(alpha = 0.2, adjust = 1.5)
h


levels(shrimp$Gender) <- c(levels(shrimp$Gender), "M")
shrimp$Gender[shrimp$Gender == "F1"] <- "F"
shrimp$Gender[shrimp$Gender == "F2"] <- "F"
shrimp$Gender[shrimp$Gender == "F3"] <- "F"
shrimp$Gender[shrimp$Gender == "M1"] <- "M"
shrimp$Gender[shrimp$Gender == "M2"] <- "M"
shrimp$Gender[shrimp$Gender == "M3"] <- "M"

tmp <- array(0, c(80, 2, 4))
for(i in 1:80){
  for(j in 1:4){
    for(k in 1:2){
      if(k == 1){ ## record counts
        if(j == 1){ ## males
          tmp[i, k, j] <- sum(shrimp$Gender[shrimp$Sample.Label == unique(shrimp$Sample.Label)[i]] == "M")
        } else if(j == 2){ ## females
          tmp[i, k, j] <- sum(shrimp$Gender[shrimp$Sample.Label == unique(shrimp$Sample.Label)[i]] == "F")
        } else if(j == 3){ ## juveniles
          tmp[i, k, j] <- sum(shrimp$Gender[shrimp$Sample.Label == unique(shrimp$Sample.Label)[i]] == "J")
        } else {
          tmp[i, k, j] <- sum(shrimp$Gender[shrimp$Sample.Label == unique(shrimp$Sample.Label)[i]] == "U")
        }
      } else { ## k == 2, record net size
        tmp[i, k, j] <- shrimp$Net.Size..m.[shrimp$Sample.Label == unique(shrimp$Sample.Label)[i]][1]
      }
    }
  }
}
tmp[, 1, ][tmp[, 2, ] == 0.5] <- 4 * tmp[, 1, ][tmp[, 2, ] == 0.5]
Counts <- as.data.frame(tmp[, 1, ])
colnames(Counts) <- c("F", "M", "J", "U")
Counts <- melt(Counts)
netSizes <- as.data.frame(tmp[, 2, ])
netSizes <- melt(netSizes)
tmp <- cbind(Counts, netSizes$value)

colnames(tmp) <- c("Gender", "Count", "Size")

i <- ggplot(tmp, aes(x = Count, fill = as.factor(Size))) + facet_wrap( ~ Gender) + geom_density(alpha = 0.2, adjust = 2.5)
i


summary(mod_gen <- glm.nb(Count ~ as.factor(Size) * as.factor(Gender), data = tmp))
plot(mod_gen)

shrimp$Gender

levels(shrimp$Gender) <- c("F", "M", "J", "U")
unique(shrimp$Gender)


