## load SWMPr
## install.packages('SWMPr')
library(SWMPr)

## Most of the data retrieval functions will not work unless the IP address is registered
## See http://cdmo.baruch.sc.edu/webservices.cfm

# retrieve metadata for all sites
site_codes()

# retrieve metadata for a single site
site_codes_ind('apa')

# all parameters for a station, most recent
all_params('hudscwq')

# get all parameters within a date range
all_params_dtrng('hudscwq', dtrng = c('09/01/2013', '10/01/2013'))

# get single parameter within a date range
all_params_dtrng('hudscwq', dtrng = c('09/01/2013', '10/01/2013'),
  param = 'do_mgl')

# single parameter for a station, most recent
single_param('hudscwq', param = 'do_mgl')

# import all paramaters for the station
# three most recent records
exdat <- all_params('apadbwq', Max = 3, trace = F)
exdat

# import sample data from package
data(apadbwq)
dat <- apadbwq

# view all attributes of dat
attributes(dat)

# view a single attribute of dat
attr(dat, 'station')

# view available methods for swmpr class
methods(class = 'swmpr')

# qaqc screen for a swmpr object, retain only '0'
qaqc(dat)

# retain all data regardless of flag
qaqc(dat, qaqc_keep = NULL)

# retain only '0' and '-1' flags
qaqc(dat, qaqc_keep = c(0, -1))

# view the number of observations in each QAQC flag
qaqcchk(dat)

# import data
data(apaebmet)
dat <- apaebmet

# select two parameters from dat
subset(dat, select = c('rh', 'bp'))

# subset records greater than or equal to a date
subset(dat, subset = '2013-01-01 0:00', operator = '>=')

# subset records within a date range, select two parameters
subset(dat, subset = c('2012-07-01 6:00', '2012-08-01 18:15'),
  select = c('atemp', 'totsorad'))

# import, qaqc removal
data(apadbwq)
dat <- qaqc(apadbwq)

# convert time series to two hour invervals
# tolerance of +/- 30 minutes for matching existing data
setstep(dat, timestep = 120, differ = 30)

# get nut, wq, and met data as separate objects
data(apacpnut)
data(apacpwq)
data(apaebmet)
swmp1 <- apacpnut
swmp2 <- apacpwq
swmp3 <- apaebmet

# combine nut and wq data by union
comb(swmp1, swmp2, method = 'union')

# combine nut and wq data by intersect
comb(swmp1, swmp3, method = 'intersect')

# combine nut, wq, and met data by nut time series, two hour time step
comb(swmp1, swmp2, swmp3, timestep = 120, method = 'apacpnut')

# get data
data(apadbwq)
dat <- apadbwq

# subset for daily decomposition
dat <- subset(dat, subset = c('2013-07-01 00:00', '2013-07-31 00:00'))

# daily decomposition of DO and plot
dc_dat <- decomp(dat, param = 'do_mgl', frequency = 'daily')
plot(dc_dat)

# get data
data(apacpnut)
dat <- apacpnut
dat <- qaqc(dat, qaqc_keep = NULL)

# decomposition of chl
decomp_cj(dat, param = 'chla_n')

## import water quality and weather data
data(apadbwq)
data(apaebmet)

## qaqc, combine
wq <- qaqc(apadbwq)
met <- qaqc(apaebmet)
dat <- comb(wq, met)

## estimate metabolism
res <- ecometab(dat, trace = FALSE)
plot_metab(res)

## import data
data(apacpnut)
dat <- qaqc(apacpnut)

## plot
plot_summary(dat, param = 'chla_n', years = c(2007, 2013))

## import data
data(apacpwq)
dat <- qaqc(apacpwq)

## plot
overplot(dat, select = c('depth', 'do_mgl', 'ph', 'turb'),
  subset = c('2013-01-01 0:0', '2013-02-01 0:0'), lwd = 2)

# plot the stations at Jacques Cousteau reserve
map_reserve('jac')

