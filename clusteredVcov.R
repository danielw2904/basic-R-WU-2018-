library(AER)
library(multiwayvcov)

data(petersen)
mod <- lm(y~x, data = petersen)

cvcov <- cluster.vcov(mod, petersen$firmid)

coeftest(mod, cvcov)
