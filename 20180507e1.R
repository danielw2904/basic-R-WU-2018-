### Welcome to the R-tutorial ###
 ##### 07.05.2018 - Basics #####

# A # tells R that the rest of the line is a comment
# R will ignore any text/command preceded by a #
# Multiple ### have the same effect
# Use this to comment on your code so you and others know what is happening

## R - a cheap calculator ##

4+5
7-3

# Space does not matter
5 *   9

# Neither horizontally nor vertically

9/6
# The modulo operator
13%%6

# Power
5^7

# Exponent
exp(1)

# Natural (!) log
log(1)

log(exp(1))
exp(log(1))

## This was fun but not very usefull ##
### Lets talk about *functions*

# All functions follow the same pattern
# name(argument1, argument2, ...)

# name: exponent, argument: 10
exp(10)

# Every operation is a function
`+`(2, 4)

# Functions can be nested and combined
exp(log(10))
2*exp(10^2 + log(6))

# But I don't know what "exp" does
?exp # This also works for all functions

## I'd like to save my result

x <- exp(1)
x
# Is equivalent to
(x <- exp(1))

# Naming is case sensitive
X <- exp(2)
X
x
y <- log(x)
y
z <- x + y
z

## WARNING: Variables can be overwritten
x
x <- 10
x

z
z <- x + y
z

## WARNING !!!: *Everything* can be overwritten
`+` <- `-`
1 + 1

log <- exp
log(1)

# Delete stuff you messed up or do not need anymore
rm(z)
z # It is gone

rm(`+`)
1 + 1 # Works again!
rm(log)
log(1)
z <- 1
# But functions are not overwritten by variables
log <- 4
log(1)
log

# The program is executed from top to bottom
# NEVER jump back in your code
# Don't do it!

## Let's go to another  dimension ##
### Vectors ###

x <- c(log(1), log(2), log(3))
x
## A lot of functions take vectors as input
a <- 1:3
a
y <- log(a)
y
sqrt(a)
exp(a)

sum(a)
mean(a)
summary(a)
range(a)

# Operations on each element
# one vector, one scalar
a^2
a + 1
a * 4
# Element wise
# Two matching vectors
b <- 4:6
a + b

c <- 4:7
a + c # Error!
# WARNING: R will recycle the shorter vector
# whenever the length of the longer is a multiple of the length of the shorter
d <- 4:9
length(a)
length(d)
length(d) %% length(a)
a + d

# Accessing elements of a vector
v <- c(8.5, 7.2, 12.012, 0.234, 2.123)
v[1]
v[3:5]
v[c(3,1,5)]
v[a]
v[1:3]
v[c(1,2,3)]
v[c(TRUE, FALSE)]

# Vector algebra
x <- rnorm(100, mean = 5, sd = 2)
head(x, 50)
str(x)
summary(x)
y <- 2*x + rnorm(100)

(m <- cbind(a, b))
(r <- rbind(a, b))

typeof(m)
class(m)
# "solve" =>  inverse
# "t" => transpose
# What is this?
solve(t(x) %*% x) %*% t(x) %*% y
lm(y~0+x) # More on this later!

# Enough with the numbers!
# Let's look at other data types

a <- c("a", "b", "c", "d", "e")
b <- LETTERS[1:5]
paste(a, b)
paste0(a, b)
# Accessing works the same
paste(a[c(1, 3, 5)], b[c(4, 2, 1)])
# Checking the type
typeof(a)
str(x)
typeof(x)
# A scalar can be a character

c <- c("a", 1)
typeof(c)
c[2]
c[2] + 1 # Error!

# Characters are surrounded by ""
char <- "3.14159"
typeof(char)
str(char)
# *Sometimes* we can coerce one type to another
# using as.numeric, as.character, as.logical, as.factor,...
sca <- as.numeric(char)
str(sca)
sca + 1

actualChar <- "R is Great!"
as.numeric(actualChar)

# Logical
t <- TRUE
ifelse(t, yes =  1, no =  2)
f <- FALSE
ifelse(f, yes = 1, no = 2)

v <- c(TRUE, FALSE, FALSE, TRUE, T, F)
v
# Don't use T of F as variable names
T
F
T <- 4
T
rm(T)

h <- as.numeric(v)
h

# Logical ?=> Numeric ?=> Character
as.logical(h)
as.character(h)
c <- as.character(v)
as.logical(c)

# Factors
f <- c("Female", "Male", "Female", "Female", "Male") # Observation
as.factor(f) # Coerce observation to factor
f <- factor(f, levels = c("Female", "Male", "Transgender")) # Create a factor variable using the vector f. 
f
levels(f)
f[6] <- "Transgender"
f[7] <- "Coffee"
str(f)

## Logical comparison
# is equal
f == "Female"
f == "Coffee"
f < 3
f
v
v == FALSE
v < 1
v == 0

numVec <- c(1:10) 
numVec <= 5
numVec >= 7
numVec == 5
!(numVec<3)
numVec != 10
# Returns if true
numVec[numVec < 5]
# Is equivalent to
lessThan5 <- numVec <5
lessThan5
numVec[lessThan5]

## EXAMPLE
female <- c(1, 0, 0, 1, 0, NA, 1, 1, 0)
gender <- factor(female, levels = c(1, 0, 2), 
                 labels = c("Female", "Male", "Transgender"))
gender
str(gender)
gender[length(gender) + 1] <- "Transgender"
gender
gender[length(gender) + 1] <- 24
# WARNING: coercing factors to numeric will start the to count from 1 and not 0
# This matter e.g. if we want to use factors in regression analysis
as.numeric(gender)
# In general its a good idea to let factors be factors

### Regression Example
# Generate some Data
set.seed(7052018)
gender <- sample(c("F", "T", "M"), size = 1000, replace = T)
gender <- factor(gender, levels = c("F", "T", "M"))
head(gender, n = 20)
income <- rnorm(1000, mean = 1000, sd = 700)
# Add some income for females
# Correct:
income[gender == "F"] <- income[gender == "F"] + 500
lm(income~gender)
# Incorrect:
gender2 <- as.numeric(gender)
str(gender2)
lm(income~gender2)

# In general there are more types for specific cases
# e.g. Dates and time 

#lubridate
?as.Date
heute <- "07.05.2018" # get year??
str( heute)
heute <- as.Date(heute, format = "%d.%m.%Y", tz = "CET")
str(heute)
heute
tutorialDates <- c("2018-05-14", "2018-05-28", "2018-06-04")
tutorialDates <- c(heute, as.Date(tutorialDates))
tutorialDates
mean(tutorialDates)
max(tutorialDates)
min(tutorialDates)
summary(tutorialDates)
Sys.Date()
tutorialDates[tutorialDates>Sys.Date()]

# Vectors are nice but we really want matrices
# Matrices behave much the same way as vectors

A <- matrix(c(1:8, 10), nrow = 3, ncol = 3, byrow = TRUE)
A
A[2, 3] # A[row, column]
A[2, ]
A[ , 2]

diag(A) # Get diagonal
eigen(A)
solve(A)

diag(3) # Create I_3
diag(1:3) # Create diag with 1, 2, 3 as elements
diag(diag(A)) # Create diagonal of diagonal of A

B <- A %*% t(A)
B
chol(B)
t(chol(B)) %*% chol(B)

##################
### EXERCISE 2 ###
##################

# Go even higher
array(1:24, c(2, 3, 4))

# All of these can only hold a single class of data
x <- c(1:10, "a")
class(x)
as.numeric(x)

# For Matrices class is "matrix" therfore use typeof()
Y <- matrix(c(1:6, "a", "b", "c"), nrow = 3)
typeof(Y)
rownames(Y) <- paste0("Obs", 1:3)
colnames(Y) <- paste0("Var", 1:3) 
Y
dim(Y)
nrow(Y)
ncol(Y)
Z <- matrix(as.numeric(Y), nrow = 3)
Z

# Solution: data.frame
df <- data.frame(var1 = 1:10, var2 = letters[1:10], stringsAsFactors = FALSE)
df
class(df)
typeof(df) # See below!
# Subsetting
df[, 1]
df[2]
df[c(1:5, 10), ]
# Get Every other
df[c(TRUE, FALSE),]
x <- df['var2']
x
class(x)
y <- df[['var2']]
y
class(y)

df$var1
df$var1 <- 10:19
df

df2 <- data.frame(gender = c("F", "M", "F", "M", "F", "F"), 
                  income = rnorm(6, 1200, 100))
df2
# Use only to access, not to overwrite
df2$g
df2['g']
df2[['g']]
df2['gender']
df2[['gender']]
df2$inc
df$var

(df3 <- df2[order(df2$gender), ])

# Lists can hold anything
data.frame(nums = 1:10, chars = letters[1:9])
L <- list(nums = 1:10, chars = letters[1:9])
L
L2 <- list(all = list(rep(L, 5)), another = matrix(1:19, nrow = 3))
L2$all[[1]][1]
L2$another

## Data Visualization

t <- seq(0,1,length.out = 100)
x <- t*sin(2*pi*t)
y <- t*cos(2*pi*t)

plot(x,y, type = "l",
     col = "red",
     main = "Last Math Test",
     sub = "... would've been easier with R",
     xlab = "x(t)", ylab = "y(t)",
     xlim = c(-1,0.5),
     ylim = c(-.75, 1))
plot(x,y)


plot(t, type = 'l', col = 'blue',
     ylim = c(-0.75, 1), lty = 2)

lines(t^2, col = 'green', lty = 3)

lines(y, col = 'red', lty = 4)

legend(x = "bottomright",legend = c("t", "t^2", "y"), 
       col = c("blue", "green", "red"), lty=c(2,3,4))

plot(df2$income~df2$gender)
boxplot(df2$income~df2$gender, ylab = "Income", xlab = "Gender")
boxplot(df2$income~df2$gender, horizontal = TRUE)
plot(df2$income, type = 'p')
plot(df2$income)
plot(df2$income, type = 'b')

# Why does this subsetting work?
plot(sort(df2$income[df2$gender=="F"]), col = 'green', axes = FALSE, 
     ylim = c((min(df2$income)-10), (max(df2$income)+10)),
     xlab = "",
     ylab = "Income")
points(sort(df2$income[df2$gender=="M"]), col = 'red')
axis(side = 1, c(1,3,4))
axis(side = 2, seq(round(min(df2$income),3), max(df2$income), length.out = 5))
box()
abline(h = seq(round(min(df2$income),3), max(df2$income), length.out = 5), lty='dashed',
       col = 'lightgray')
abline(v = 1:4, lty='dashed', col = 'lightgray')
legend(x = "topleft", legend = c("F", "M"), col = c("green", "red"), pch = 1)
text(x = c(1:4), y = sort(df2$income[df2$gender=="F"]), labels= 'F', pos = 1)
text(x = c(1:2), y = sort(df2$income[df2$gender=="M"]), labels = "M", pos = 3)

(counts <- table(df2$gender))

barplot(counts, col = c('green', 'blue'))

df2$area <- sample(c("city", "rural"), size = 6, replace = TRUE)
(ncounts <- table(df2$gender,df2$area ))
mosaicplot(ncounts, color = TRUE, main = "Area")
png(filename = "/home/daniel/greatplot.png")
plot(density(df2$income[df2$gender=='F']),
     xlim = c(700, 1600), 
     ylim = c(0, 0.03),
     col = 'blue', 
     main = 'Density of Income by Gender',
     lty = 2)
lines(density(df2$income[df2$gender=='M']), col = 'red')
legend(x='topleft', 
       legend=c("F", "M"), 
       lty=c(2,1), 
       col = c('green', 'red'))
dev.off()
par(mfrow=c(3,1))
plot(density(df2$income), main = "All")
plot(density(df2$income[df2$gender=="M"]), main = "M")
plot(density(df2$income[df2$gender=="F"]), main = "F")
par(mfrow=c(1,1))
##################
### Exercise 3 ###
##################

### Packages ###

# Set seed for reproducible results
set.seed(1234)
# install.packages(c("tidyverse", "AER"))
library(tidyverse)
df3 <- data.frame(gender = factor(sample(c("M", "F"), size = 100, replace = TRUE)),
                  parentIncomeAt15 = rnorm(100, 500, 100),
                  area = factor(sample(c("city", "rural"), size = 100, replace = TRUE)),
                  education = factor(sample(c("Primary", "Secondary", "Tertiary"), size = 100, replace = TRUE)),
                  stringsAsFactors = FALSE)

df3$income <- 500 + 
  10 * df3$parentIncomeAt15 -
  0.01 * df3$parentIncomeAt15^2 + # infliction point at 500
  900 * ifelse(df3$education=="Tertiary", 1, 0) + 
  500 * ifelse(df3$education == "Secondary", 1, 0) +
  100 * ifelse(df3$gender == "F", 1, 0) +
   50 * ifelse(df3$area == "city", 1, 0) +
  200 * ifelse(df3$education == "Tertiary" & df3$gender == "F", 1, 0) +
  rnorm(100, sd = 10)

ggplot(df3, aes(x = income)) +
  geom_histogram(binwidth = 30, aes(y = ..count..)) 

ggplot(df3) +
  geom_density(aes(x = income, fill = area, linetype = education), alpha = 0.6) +
  labs(title = 'Density of Income', 
       subtitle = 'by area and education', 
       caption = "Source: Daniel's made up data") 
  

ggplot(df3)+
  geom_boxplot(aes(x = education, y = income)) +
  facet_wrap(~gender)

ggplot(df3)+
  geom_line(aes(x = 1:100, y = sort(income), color = area))+
  geom_point(aes(x = 1:100, y = sort(income), color = area, shape = area))+
  labs(x="", y='Income', title='Income', subtitle="By Area") 

ggplot(df3) +
  geom_point(aes(x = parentIncomeAt15, y = income, color=gender)) +
  labs(x = "Parents' income at 15", title = "Intergenerational mobility") +
  facet_wrap(~area) +
  theme(legend.title = element_blank()) +
  scale_x_continuous(breaks = c(200, 400, 600, 700))

ggplot(df3) +
  geom_point(aes(x = parentIncomeAt15, y = income, color=gender, shape = gender)) +
  labs(x = "Parents' income at 15", title = "Intergenerational mobility") +
  facet_grid(education~area)

### Econometrics
# Distributions
rnorm(10, mean = 10, sd = 1) # random numbers
dx <- density(rnorm(1000), bw = "nrd0", adjust = 1, kernel = "epanechnikov") # density estimate
dx
plot(dx)

pnorm(5, mean = 5, sd = 10) # CDF
plot(function(x)pnorm(x, 5, 1), xlim = c(0, 10))

dnorm(5, mean = 5, sd = 1) # PDF
plot(function(x)dnorm(x, 5, 1), xlim = c(0, 10))

# Same for other distributions
rt(10, df = 10)
rf(10, df1 = 10, df2 = 5)
runif(10)
rgamma(10, shape = 4)
rchisq(10, df = 1)
rlnorm(10)
?distributions

# Create some data
x <- runif(100, 0, 10)
y <- 5 + 3*x + rnorm(100)

# Formula
model <- y ~ x
class(model)

# OLS
ols1 <- lm(model)
ols1

summary(ols1)
AIC(ols1)
BIC(ols1)

# Testing Model assumptions

# Heteroskedasticity
residuals1 <- residuals(ols1)
fitted1 <- fitted(ols1)
plot(residuals1~fitted1)
# Normality
qqnorm(residuals1)
shapiro.test(residuals1)

# Autocorrelation
acf(residuals1)

# Optional: What are we doing?
X <- cbind(1,x)
xxinv <- solve(crossprod(X))
beta <- xxinv %*% crossprod(X, y)
df <- length(y) - ncol(X)
resids <- y - X %*% beta
Sigmasq <- as.numeric(crossprod(resids)/df)
vcov <- xxinv * Sigmasq
sds <- sqrt(diag(vcov))
tval <- beta/sds
pval <- 2*pt(abs(tval), df = df, lower.tail = FALSE)
beta
sds
pval

# White standard errors
u2 <- diag(residuals1^2)
X <- cbind(1,x)
xxinv <- solve(crossprod(X))
vcov <- xxinv %*% crossprod(X, u2) %*% X %*% xxinv
colnames(vcov) <- c("intercept", "x")
robse <- sqrt(diag(vcov))
names(robse) <- c("intercept", "x")
robse
tval_robust <- beta/robse
p_robust <- 2*pt(abs(tval_robust), df = df, lower.tail = FALSE)
tval_robust
p_robust

# External Packages
library(AER)
coeftest(ols1, vcov = vcov)

library(stargazer)
stargazer(ols1, se = list(robse2), type = "text")

(ct <- coeftest(ols1, vcov = vcovHC, type = "HC1"))
stargazer(ct, type = "text")

# Heteroskedasticity

x <- runif(1000, 10, 100)
y <- 40 +  5 * x + rnorm(1000, sd = x)
mod <- lm(y~x)
plot(y~x)
lines(fitted(mod)~x, col = 'red')
plot(residuals(mod)~fitted(mod))
abline(h = c(60, -60), col = 'red')

summary(mod)
ct <- coeftest(mod, vcov = vcovHC, type = "HC1")
stargazer(ct, type = 'text')

# "Real Model"

model <- income ~ parentIncomeAt15 + gender + area + education + gender*education
ols <- lm(model, data = df3)
summary(ols)
ols2 <- update(ols, . ~ . + I(parentIncomeAt15^2))
summary(ols2)

library(AER)
(ct <- coeftest(ols2, vcov = vcovHC, type = "HC3"))
library(stargazer)
stargazer(ct, type = "text")

#install.packages('latex2exp')
library(latex2exp)
ggplot(df3) +
  geom_smooth(aes(x = parentIncomeAt15, y = income),
              method = 'loess', level = 0.99) +
  labs(y = "Parents' income at age 15", title = "Intergenerational mobility", 
       subtitle = TeX('$income = \\alpha + \\beta_1 parent Income + \\beta_2 (parent Income)^2 + \\epsilon $'))+
  annotate("text", x = 500, y = 2600, 
           label = as.character(TeX('Income $ = \\alpha + \\beta_1 parent Income + \\beta_2 (parent Income)^2 + \\epsilon $')),
           parse = TRUE) +
  labs(caption = TeX('Income $ = \\alpha + \\beta_1 parent Income + \\beta_2 (parent Income)^2 + \\epsilon $'))

# Probit/logit
data("Affairs")
affairs <- Affairs
head(affairs)
affairs$hadAffair <- ifelse(affairs$affairs > 0, 1, 0)
mod <- hadAffair ~ gender + age + yearsmarried + children + religiousness + education + rating
out <- glm(mod, family = binomial(link = 'logit'), data = affairs)
summary(out)

# Poisson
# E[Y|X] = exp(x'B)
modCount <- update(mod, affairs ~ . )
outCount <- glm(modCount, family = poisson, data = affairs)
summary(outCount)
dispersiontest(outCount)

# Tobit
outTobit <- tobit(modCount, left = 0, data = affairs)
summary(outTobit)

# Time Series Econometrics

library(eurostat)
library(forecast)
youthUne <- get_eurostat("tipslm80", filters = list(geo = "AT", unit = "PC_ACT"))
ggplot(youthUne) +
  geom_line(aes(x = time, y = values))
yu_ts <- youthUne$values[order(youthUne$time)]
yu_ts <- ts(yu_ts, start = 1995, end = 2017, frequency = 1)
armod <- Arima(yu_ts, order = c(1,0,0))
summary(armod)
plot(forecast(armod, 30))

data("UKDriverDeaths")
plot(UKDriverDeaths, main = 'Driver Deaths')

season <- seasonaldummy(UKDriverDeaths)
trend <- 1:length(UKDriverDeaths)
# De-trend and de-season
mod <- lm(UKDriverDeaths ~ season + trend )
summary(mod)
res <- residuals(mod)
plot(res, type = "l")
Acf(res)
Pacf(res)

armod <- Arima(res, order = c(1,0,0))
summary(armod)
plot(forecast(armod, 10))

mod2 <- auto.arima(UKDriverDeaths)
summary(mod2)
plot(forecast(mod2, 30))

# some more visualization

une_s <- search_eurostat("unemployment", type = "table")
une_une_c <- une_s$code[une_s$title == 'Total unemployment rate']
une_d <- get_eurostat(id = une_c, filters = )
View(une_d)
une_d <- une_d[une_d$unit == 'PC_ACT', ]
countries <- c("AT", "DE", "DK")
une_d <- une_d[une_d$geo %in% countries, ]

une_wide <- spread(une_d, key = geo, value = values)
une_wide$time <- format(une_wide$time, "%Y")
library(xtable)
une_t <- xtable(une_wide[, 4:7], caption = "Unemployment rate")
print(une_t, type = 'latex', include.rownames = FALSE)

library(knitr)
kable(une_wide[, 4:7], format = 'markdown')

ggplot(une_wide) +
  geom_point(aes(x = time, y = AT)) +
  geom_point(aes(x = time, y = DE), color = 'red')
 
une_long <- gather(une_wide[, 4:7], key = "Country", value = "une", -time) 
une_long
ggplot(une_long) +
  geom_point(aes(x = time, y = une, color = Country)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
