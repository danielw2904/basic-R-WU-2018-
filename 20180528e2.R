rm(list = ls())
library(eurostat)
une_s <- search_eurostat("unemployment", type = "table")
une_c <- une_s$code[une_s$title == 'Total unemployment rate']
une_d <- get_eurostat(id = une_c, filters = )
View(une_d)
une_d <- une_d[une_d$unit == 'PC_ACT', ]
countries <- c("AT", "DE", "DK")
une_d <- une_d[une_d$geo %in% countries, ]

# Tidy Data
# Long and Wide
library(tidyverse)
une_wide <- spread(une_d, key = geo, value = values)
une_wide$time <- format(une_wide$time, "%Y")
une_wide[, 4:6]

library(xtable)
une_t <- xtable(une_wide[, 4:7], caption = "Unemployment rate")
print(une_t, type = 'latex', include.rownames = FALSE)

library(knitr)
kable(une_wide[, 4:7], format = "markdown")

ggplot(une_wide) +
  geom_point(aes(x = time, y = AT, color = 'AT')) +
  geom_point(aes(x = time, y = DE, color = 'DE')) +
  geom_point(aes(x = time, y = DK, color = 'DK')) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_color_manual("Country", 
                     breaks = c("AT", "DE", "DK"),
                     values = c('red', 'green', 'blue')) +
  labs(y = "une")

(une_long <- gather(une_wide[, 4:7], key = "Country", value = "une", -time)) # Gather all BUT time
gather(une_wide[, 4:7], key = 'Country', value = 'une', AT, DE, DK) # Gather AT, DE and DK

ggplot(une_long) +
  geom_point(aes(x = time, y = une, color = Country)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

# Renaming
(une_long2 <- rename(une_long, unemployment = une))

# Sorting
arrange(une_long2, desc(unemployment), desc(time))

# Aggregate
aggregate(une_long$une, by = list(time = une_long$time), FUN = mean)

# Medium Data
library(data.table)
group <- sample(LETTERS, size = 1e7, replace = TRUE)
values <- rnorm(1e7, 10, 100) 
dt <- data.table(group, values)
tables()

dt[, .(Mean = mean(values), SD = sd(values), N = .N), by = group]

dt2 <- dt
dt2[, Sqs := values^2]
dt2[, lapply(.SD, mean), by = group]

setkey(dt, group)

library(microbenchmark)
microbenchmark(
dt['A'],
times = 10)

df <- data.frame(dt)
microbenchmark(
head(df[df$group == 'A', ], 5),
times = 10)

microbenchmark(
  tapply(df$values, df$group, mean),
  times = 10)

microbenchmark(
  dt[, mean(values), by = group],
  times = 10
)

# Data import
# SPSS, Stata, SAS
library(haven)

data(mtcars)
head(mtcars)

setwd("~/pCloudDrive/Daniel/studiVertretung/rKurs")
write_dta(data = mtcars, "./mtcars.dta")
write_sas(mtcars, "./mtcars.sas")
write_sav(mtcars, "./mtcars.sav")

(df <- read_dta("mtcars.dta"))
(df <- read_sas("mtcars.sas", catalog_file = NULL))
(df <- read_sav("mtcars.sav"))

# CSV built in
mtcars2 <- mtcars
mtcars2[1,1] <- "NOT"
mtcars2[, 2] <- as.character(mtcars2[, 2])
write.table(mtcars2, "mtcars.csv", sep = ";", dec = ",", quote = TRUE)

head(read.csv("mtcars.csv", sep = ';', dec = ',', 
         header = TRUE, na.strings = 'NOT', 
         quote = '"', stringsAsFactors = FALSE))

# Excel Not recommended
library(readxl)
read_excel("Test_Data.xlsx", sheet = "Sheet1", na = "NA" , col_names = TRUE)

### Loops ###

# Simulate a random walk
set.seed(12)
start1 <- Sys.time()
mu <- 5
N <- 1e6
series <- matrix(NA, nrow = N)
series[1] <- mu 
for(i in 2:N){
series[i, 1] <- series[(i-1), ] + rnorm(1, sd = 4)
}
end1 <- Sys.time()

set.seed(12)
start2 <- Sys.time()
path <- rnorm(1e6-1, sd = 4)
series2 <- cumsum(c(5, path))
end2 <- Sys.time()

cat("Loop:", end1-start1, "; Vectorized:", end2-start2)

plot(ts(series2), main = "Loop and vectorized", col = 'green')
lines(ts(series), lty='dotted', col = 'red')

# Basic Markov Chain

# Define transition probabilities
p11 <- 0.4
p12 <- 1-p11

p21 <- 0.5
p22 <- 1-p21

cat("d1:", (d1 <- p21/(p12 + p21)))
cat("d2:", (d2 <- p12/(p12 + p21)))

# Number ahead forecast (enough to stabilize)
nahead <- 20
# How often should this be done
nsave <- 10000
Sstore <- matrix(NA, nsave)
ProbStore <- matrix(NA, 2, nsave)
for (i in 1:nsave){
 S <- "a"
  for(j in 1:nahead){
    if (S == "b"){
      S <- sample(c("a", "b"), size = 1, prob = c(p21, p22))
    }else if (S == "a"){
      S <- sample(c("a", "b"), size = 1, prob = c(p11, p12))
    }
  }
  Sstore[i, ] <- S
  ProbStore[1, i] <- sum(Sstore=="a", na.rm = TRUE)/i
  ProbStore[2, i] <- sum(Sstore=="b", na.rm = TRUE)/i}

plot(ProbStore[1,],type = "l")
lines(ProbStore[2,],type = "l", col = "red")

# Final Probability with 500 burn-in
Sstore <- Sstore[500:nsave,]
(probA <- sum(Sstore=="a")/(nsave-500))
(probB <- sum(Sstore=="b")/(nsave-500))

# Funcitonal programming
# Example from 
# http://adv-r.had.co.nz/Functional-programming.html

myFun <- function(x, y, z){
  sum <- x + y + z
  return(sum)
}

myFun(1,2,3)

stationaryDist <- function(gamma){
  p12 <- gamma[1,2]
  p21 <- gamma[2,1]
  cat("d1:", (d1 <- p21/(p12 + p21)))
  cat("d2:", (d2 <- p12/(p12 + p21)))
}

G <- matrix(c(0.4, 0.5, 0.6, 0.5), nrow = 2)
stationaryDist(G)

une_d <- get_eurostat(id = une_c)
View(une_d)
une_d <- une_d[une_d$unit == 'PC_ACT', ]
countries <- c("AT", "DE", "DK")
une_d <- une_d[une_d$geo %in% countries, ]

euplot <- function(data, title="", y=""){
  require(ggplot2)   # Warning instead of error 
  ggplot(data = data, aes(x = time, y = values, color = geo, group = geo, shape =geo)) +
    ggtitle(title) +
    geom_path() +
    geom_point(size = 2) +
    labs(y=y, x="year") +
    scale_shape_discrete(name = "Country") + scale_color_discrete(name = "Country")+
    labs(caption = "Source: Eurostat 2017, own calculations")
}

euplot(une_d, title = 'Unemployment', y = 'percent of working population')

euplot(une_d) + 
  ggtitle("Unemployment Extra") +
  labs(y = "Percentage")

summary2 <- function(x){
  print(c(Mean = mean(x), Range = range(x), Median = median(x), IQR = IQR(x), LogLik = sum(log(x))))
}
summary2(mtcars$mpg)

# Vectorize
# *apply
# Example based on 
# http://adv-r.had.co.nz/Functional-programming.html

for(i in 1:ncol(mtcars)){
  summary2(mtcars[, i])
}

# Apply to all columns
apply(mtcars, 2, summary2)

# Apply to a list
mtSummaryL <- lapply(mtcars, summary2)
mtSummarySimplified <- sapply(mtcars, summary2)

# lapply in user function
# see https://stackoverflow.com/questions/3057341/how-to-use-rs-ellipsis-feature-when-writing-your-own-function
myFun2 <- function(...){
  input <- list(...)
  print("That's what the input looks like:")
  print(input)
  print("That's the output:")
  lapply(input, mean)
}

myFun2(1:10, 2:11, 3:12)

# Apply by group
(wtByCyl <- tapply(X = mtcars$wt, INDEX = mtcars$cyl, FUN = summary2))

library(microbenchmark)
D <- matrix(rnorm(1e4, 10, 20), ncol = 1e2)

sumSave <- matrix(NA, ncol = 6, nrow = ncol(D))
loop <- microbenchmark(
for(i in 1:ncol(D)){
  sumSave[i, ] <- summary2(D[,i])
}, times = 10)

gc()
ap <- microbenchmark(
  sumSave2 <- apply(D, 2, summary2),
  times = 10)

loop
ap

funs <- list(function(x) x^2, function(x)x/2, function(x)x/4)
m <- list(1:10, 1:10, 1:10)

# (c) Peter Knaus
applier <- function(FUN, arg){
  # Diese funktion ist arg fun
  FUN(arg)
}
mapply(applier, funs, mtcars[,1:3])

mtcars[, 1:3] <- lapply(mtcars[,1:3], as.character)
str(mtcars)
mtcars[, 1:3] <- lapply(mtcars[,1:3], as.numeric)
str(mtcars)
