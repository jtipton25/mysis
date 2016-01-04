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

## 
## Juvenile ratio data
##

fullData <- read.csv("~/mysis/data/lengthData.csv", skip=4)

males <- rep(0, 80)
females <- rep(0, 80)
juveniles <- rep(0, 80)
unknowns <- rep(0, 80)
size <- rep(0, 80)
date_tmp <- rep(0, 80)
station <- rep(0, 80)

idx <- 1
for(i in unique(fullData$Sample.Label)){
  tmp_data = subset(fullData, fullData$Sample.Label == i)
  tmp = substr(as.character(tmp_data$Gender), 0, 1)
  males[idx] <- sum(tmp == "M")
  juveniles[idx] <- sum(tmp == "J")
  females[idx] <- sum(tmp == "F")
  unknowns[idx] <- sum(tmp == "U")
  size[idx] <- tmp_data$Size..m.[1]
  date_tmp[idx] <- tmp_data$Date[1]
  station[idx] <- tmp_data$Station[1]
  idx <- idx + 1
}
prop <- juveniles / (juveniles + females + males + unknowns)

## Correct for mislabeling
date_tmp[70] <- 4
date_tmp[72] <- 4
date <- rep(0, 80)
date[date_tmp == 1] <- "July"
date[date_tmp == 2] <- "July"
date[date_tmp == 3] <- "August"
date[date_tmp == 4] <- "September"

date <- factor(date, levels=c("July", "August", "September"))
size[size == 0.5] <- "Small Net"
size[size == 1.0] <- "Large Net"
net <- factor(size, levels=c("Small Net", "Large Net"))
station <- factor(station)
mysisJuvenileData <- data.frame(prop=prop, date=date, net=net, station=station)

## Save RData file
save(mysisJuvenileData, file='~/mysis/data/mysisJuvenileData.RData')
## remove data to elimate conflicts
rm(list=ls())

## 
## Sex ratio data
##

fullData <- read.csv("~/mysis/data/lengthData.csv", skip=4)

males <- rep(0, 80)
females <- rep(0, 80)
unknowns <- rep(0, 80)
size <- rep(0, 80)
date_tmp <- rep(0, 80)
station <- rep(0, 80)

idx <- 1
for(i in unique(fullData$Sample.Label)){
  tmp_data = subset(fullData, fullData$Sample.Label == i)
  tmp = substr(as.character(tmp_data$Gender), 0, 1)
  males[idx] <- sum(tmp == "M")
  females[idx] <- sum(tmp == "F")
  unknowns[idx] <- sum(tmp == "U")
  size[idx] <- tmp_data$Size..m.[1]
  date_tmp[idx] <- tmp_data$Date[1]
  station[idx] <- tmp_data$Station[1]
  idx <- idx + 1
}
prop <- females / (females + males)
prop_corrected <- (females + 0.5 * unknowns) / 
  (females + males + unknowns)
prop[prop == 0] <- 0.01
prop_corrected[prop_corrected == 0] <- 0.01
## Correct for mislabeling
date_tmp[70] <- 4
date_tmp[72] <- 4
date <- rep(0, 80)
date[date_tmp == 1] <- "July"
date[date_tmp == 2] <- "July"
date[date_tmp == 3] <- "August"
date[date_tmp == 4] <- "September"

date <- factor(date, levels=c("July", "August", "September"))
size[size == 0.5] <- "Small Net"
size[size == 1.0] <- "Large Net"
net <- factor(size, levels=c("Small Net", "Large Net"))
station <- factor(station)
mysisSexData <- data.frame(prop=prop, date=date, net=net, station=station)

## Save RData file
save(mysisSexData, file='~/mysis/data/mysisSexData.RData')
## remove data to elimate conflicts
rm(list=ls())

## 
## Sex count data
##

fullData <- read.csv("~/mysis/data/lengthData.csv", skip=4)

males <- rep(0, 80)
females <- rep(0, 80)
juveniles <- rep(0, 80)
unknowns <- rep(0, 80)
size <- rep(0, 80)
date_tmp <- rep(0, 80)
station <- rep(0, 80)

idx <- 1
for(i in unique(fullData$Sample.Label)){
  tmp_data = subset(fullData, fullData$Sample.Label == i)
  tmp = substr(as.character(tmp_data$Gender), 0, 1)
  males[idx] <- sum(tmp == "M")
  juveniles[idx] <- sum(tmp == "J")
  females[idx] <- sum(tmp == "F")
  unknowns[idx] <- sum(tmp == "U")
  size[idx] <- tmp_data$Size..m.[1]
  date_tmp[idx] <- tmp_data$Date[1]
  station[idx] <- tmp_data$Station[1]
  idx <- idx + 1
}

## construct count vector
count <- c(males, females, juveniles, unknowns)
gender <- factor(rep(1:4, each=80))
## Correct for mislabeling
date_tmp[70] <- 4
date_tmp[72] <- 4
date <- rep(0, 80)
date[date_tmp == 1] <- "July"
date[date_tmp == 2] <- "July"
date[date_tmp == 3] <- "August"
date[date_tmp == 4] <- "September"

date <- factor(date, levels=c("July", "August", "September"))
size[size == 0.5] <- "Small Net"
size[size == 1.0] <- "Large Net"
net <- factor(size, levels=c("Small Net", "Large Net"))
station <- factor(station)
count[rep(net, times=4) == 0.5] <- count[rep(net, times=4) == 0.5] * 4
mysisSexCountData <- data.frame(count=count, date=rep(date, times=4), 
                                net=rep(net, times=4), gender=gender, 
                                station=rep(station, times=4))

## Save RData file
save(mysisSexCountData, file='~/mysis/data/mysisSexCountData.RData')
## remove data to elimate conflicts
rm(list=ls())
