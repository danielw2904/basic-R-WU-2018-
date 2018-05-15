library(AER)
library(multiwayvcov)

data(petersen)
mod <- lm(y~x, data = petersen)
stargazer(mod, type = 'text')
cvcov <- cluster.vcov(mod, petersen$firmid)

library(stargazer)

stargazer(coeftest(mod, cvcov), type = 'text')

install.packages("cowplot")
library(cowplot)
