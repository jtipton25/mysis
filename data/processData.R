## File to process the CSV files to RData files

##
## Created 12_15_2015
## Edited 12_15_2015
##

##
## Count Data
##

data_one_way <- read.csv("~/mysis/data/ttestData.csv", skip = 3, nrows = 40, header = TRUE)
# data_one_way$full ## count per m^2 for full net
# data_one_way$half ## count per m^2 for half net
count <- c(data_one_way$full.5,  ## count for full net
           4 * data_one_way$half.5) ## count for half net adjusted for radius
date <- rep((data_one_way$month), 2)
date[date == 7] <- "July"
date[date == 8] <- "August"
date[date == 9] <- "September"
date <- factor(date, levels=c("July", "August", "September"))
net <- factor(c(rep("Large Net", 40), rep("Small Net", 40)), levels=c("Small Net", "Large Net"))
station <- factor(data_one_way$Station)
mysisCountData <- data.frame(count=count, date=date, net=net, station=station)

## Save RData file
save(mysisCountData, file="~/mysis/data/mysisCountData.RData")   
## remove all data to avoid conflicts
rm(list=ls())

## 
## Length Data
##

data_length <- read.csv("~/mysis/data/lengthData.csv", skip = 4, header = TRUE)
net <- data_length$Size..m.
net[net== 0.5] <- "Small Net"
net[net== 1.0] <- "Large Net"
net <- factor(net, levels=c("Small Net", "Large Net"))
y <- data_length$Length..mm.
## remove whitespace from string, then truncate into 
## broad gender classes, not subclasses
library(stringr)
gender <- as.factor(substr(str_trim(as.character(data_length$Gender)), 0, 1))

date <- as.numeric(data_length$Date)
date[date == 2] <- 1      ## collapse into monthly values
date[date == 3] <- 2
date[date == 4] <- 3
date <- as.factor(date)
plot(as.numeric(date), type='l', main = "Is this a data error? - yes?")
## correct date mislabeling in the raw data
date <- as.numeric(data_length$Date)
date[7330:7519] <- 4
date[7559:7773] <- 4
date[date == 2] <- 1      ## collapse into monthly values
date[date == 3] <- 2
date[date == 4] <- 3
# date <- as.factor(date)
plot(as.numeric(date), type='l', main = "Data looks better")
date[date == 1] <- "July"
date[date == 2] <- "August"
date[date == 3] <- "September"
date <- factor(date, levels=c("July", "August", "September"))
station <- factor(data_length$Station)
label <- factor(data_length$Sample.Label)
mysisLengthData <- data.frame(y=y, date=date, net=net, gender=gender,
                              station=station, label=label)
## Save RData file
save(mysisLengthData, file='~/mysis/data/mysisLengthData.RData')
## remove data to elimate conflicts
rm(list=ls())

