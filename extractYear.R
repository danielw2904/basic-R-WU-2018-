tutorialDates <- c("2018-05-14", "2018-05-28", "2018-06-04")
# SOme work here
tutorialDates <- as.Date(tutorialDates, format="%Y-%m-%d")

format(tutorialDates, "%Y")

library(lubridate)
year(tutorialDates)
