## Download and prepare the data for Acemoglu and Johnson (2007)
#
# Citation:
# Acemoglu, Daron, and Simon Johnson. "Disease and development: the effect of
# life expectancy on economic growth." Journal of political Economy 115, no. 6
# (2007): 925-985.
#
# Available at:
# https://www.journals.uchicago.edu/doi/full/10.1086/529000
library(data.table)

#### Download the data and unzip
# This probably has to be done by hand since it needs journal access
dir.create(here::here("data/AcemogluJohnson2007"), recursive = TRUE)
url <- "https://www.journals.uchicago.edu/doi/suppl/10.1086/529000/suppl_file/32236data.zip"
download.file(url, destfile = here::here("data/AcemogluJohnson2007/32236data.zip"))

unzip(here::here("data/AcemogluJohnson2007/32236data.zip"), exdir = here::here("data/AcemogluJohnson2007"))

#### Read and process into first differences
aj <- haven::read_dta(here::here("data/AcemogluJohnson2007/MASTERSETFORREGRESSIONS.dta")) %>%
  data.table()

aj2 <- na.omit(aj[sjbasesamplenoncomm == 1 & (sample40 == 1 & sample80 == 1) & year %in% c(1940, 1980), .(
  loggdppcmadd, loglifeexpect,
  compsjmhatit, globmort1000,
  loggdppcmadd301950, loggdppcmadd301960, loggdppcmadd301970, loggdppcmadd301980,
  ctry, ctrycluster, year, yr1940, yr1980)])

# Manual FD
ajfd <- aj2[, .(
  loggdppcmadd = loggdppcmadd - shift(loggdppcmadd),
  loglifeexpect = loglifeexpect - shift(loglifeexpect),
  compsjmhatit = shift(compsjmhatit),
  globmort1000 = shift(globmort1000),
  loggdppcmadd301980),
  by = "ctry"][!is.na(loggdppcmadd)]

#### Test outcome to ensure it matches
df <- ajfd[, .(loggdppcmadd, loglifeexpect, loggdppcmadd301980, compsjmhatit)]
lm1 <- lm(loglifeexpect ~ loggdppcmadd301980 + compsjmhatit, data = df)
df$y1_hat <- fitted(lm1)
lm2 <- lm(loggdppcmadd ~ loggdppcmadd301980 + y1_hat, data = df)
stopifnot((coef(lm2)[3] - -1.589)^2 < 1e-6)

#### Save dataset
AcemogluJohnson2007 <- ajfd
usethis::use_data(AcemogluJohnson2007, overwrite = TRUE)
