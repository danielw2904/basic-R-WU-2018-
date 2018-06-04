tutorialDates <- c("2018-05-14", "2018-05-28", "2018-06-04")

tutorialDates <- as.Date(tutorialDates, format="%Y-%m-%d")

format(tutorialDates, "%Y")

library(lubridate)
year(tutorialDates)
# Wow that's great