## ----setup, cache = F, echo = F------------------------------------------
library(knitr)
# set global chunk options
opts_chunk$set(fig.align = 'center', message = F, dev = 'pdf', dev.args = list(family = 'serif'), fig.pos = '!ht', warning = F, background = 'white', highlight = FALSE, prompt = TRUE, size = 'small')
options(replace.assign=TRUE,width=70,digits=1)

## ----eval = F, message = F-----------------------------------------------
## install.packages('SWMPr')
## library(SWMPr)

## ----eval = T, echo = F, message = F, cache = F--------------------------
devtools::load_all('M:/docs/SWMPr')

## ----results = 'asis', echo = FALSE--------------------------------------
funcs <- c('all\\_params', 'all\\_params\\_dtrng', 'import\\_local', 'single\\_param', 'site\\_codes', 'site\\_codes\\_ind')
funcs <- paste0('\\texttt{', funcs, '}')

descrips <- c(
  'Retrieve records starting with the most recent at a given station, all parameters.  Wrapper to \\texttt{exportAllParamsXMLNew} function on web services.',
  'Retrieve records of all parameters within a given date range for a station.  Optional argument for a single parameter. Wrapper to \\texttt{exportAllParamsDateRangeXMLNew}.', 
  'Import files from a local path.  The files must be in a specific format, such as those returned from the \\ac{CDMO} using the zip downloads option.',
  'Retrieve records for a single parameter starting with the most recent at a given station.  Wrapper to \\texttt{exportSingleParamXMLNew} function on web services.', 
  'Metadata for all stations, wrapper to \\texttt{exportStationCodesXMLNew} function on web services.',
  'Metadata for all stations at a single site, wrapper  to \\texttt{NERRFilterStationCodesXMLNew} function on web services.'
  )

to_tab <- data.frame(Functions = funcs, Description = descrips, stringsAsFactors = F)

library(Hmisc)

latex(
  to_tab[, 'Description', drop = F],
  file = '',
  caption = "Retrieval functions available from the SWMPr package. Full documentation for each function is in the help file (e.g., execute \\texttt{?all\\_params} for individual functions or \\texttt{help.search(`retrieve', package = `SWMPr')} for all).",
  rowlabel = 'Function',
  colheads = 'Description',
  rowname = to_tab$Functions,
  caption.loc = 'top',
  col.just=c("p{3.5in}"), 
  label = 'tab:retrieve', 
  table.env = FALSE
  )

## ----eval = F------------------------------------------------------------
## # retrieve metadata for all sites
## site_codes()
## 
## # retrieve metadata for a single site
## site_codes_ind('apa')

## ----eval = F------------------------------------------------------------
## # all parameters for a station, most recent
## all_params('hudscwq')
## 
## # get all parameters within a date range
## all_params_dtrng('hudscwq', c('09/10/2012', '02/8/2013'))
## 
## # get single parameter within a date range
## all_params_dtrng('hudscwq', c('09/10/2012', '02/8/2013'),
##   param = 'do_mgl')
## 
## # single parameter for a station, most recent
## single_param('hudscwq', 'do_mgl')

## ----eval = F------------------------------------------------------------
## # import local data for apaebmet
## 
## # this is an example path with the decompressed csv files
## path <- 'C:/my_path/'
## 
## # import, do not include file extension
## import_local(path, 'apadbwq')

## ------------------------------------------------------------------------
# import all paramaters for the station
# three most recent records
exdat <- all_params('apadbwq', Max = 3, trace = F)
exdat

## ----eval = F, cache = T-------------------------------------------------
## # import sample data from package
## data(apadbwq)
## dat <- apadbwq
## 
## # view all attributes of dat
## attributes(dat)
## 
## # view a single attribute of dat
## attr(dat, 'station')

## ----results = 'asis', echo = FALSE--------------------------------------
# create the attributes table
data(apadbwq)
dat <- apadbwq

# view all attributes of dat
atts <- paste0('\\texttt{', names(attributes(dat)), '}')
atts <- gsub('_', '\\_', atts, fixed = T)
atts_class <- c('character', 'integer', 'character', 'character', 'character', 'logical', 'POSIXct', 'character', 'character')
atts_desc <- c(
  'Column names of the entire data set, inherited from the \\texttt{data.frame} object class',
  'Row names of the data set, inherited from the \\texttt{data.frame} object class',
  'Class of the data object indicating \\texttt{swmpr} and \\texttt{data.frame}',
  'Station identifier used by \\ac{NERRS} as a string with 7 or 8 characters',
  "Character vector of column names for data parameters, e.g., \\texttt{`do\\_mgl'}",
  'Indicates if \\ac{QAQC} columns are present in the raw data',
  'Start and end dates for the raw data',
  'Timezone of the station using the city/country format\\textsuperscript{\\textit{a}}',
  'Class of the \\texttt{datetimestamp} column, usually POSIXct unless data have been aggregated'
)

to_tab <- data.frame(Attributes = atts, Class = atts_class, Description = atts_desc, stringsAsFactors = F)
foot <- c('\\textsuperscript{\\textit{a}}\\footnotesize Time zones that do not observe daylight savings are used for \\texttt{swmpr} objects and may not be cities in the United States.  For example, ``America/Jamaica" is used for Eastern Standard Time.')

library(Hmisc)

latex(
  to_tab[, -1, drop = F],
  file = '',
  caption = "Attributes of a \\texttt{swmpr} object that describe characteristics of the data.",
  rowlabel = 'Attributes',
  colheads = c('Class', 'Description'),
  rowname = to_tab$Attributes,
  caption.loc = 'top',
  insert.bottom = foot,
  col.just=c("p{0.75in}", "p{3.25in}"), 
  label = 'tab:attributes', 
  table.env = FALSE
  )

## ----eval = F, cache = T-------------------------------------------------
## # view available methods for swmpr class
## methods(class = 'swmpr')

## ----results = 'asis', echo = FALSE--------------------------------------

funcs <- c('comb', 'qaqc', 'qaqcchk', 'rem\\_reps', 'setstep', 'subset')
funcs <- paste0('\\texttt{', funcs, '}')

descrips <- c(
  'Combines \\texttt{swmpr} objects to a common time series using setstep, such as combining the weather, nutrients, and water quality data for a single station. Only different data types can be combined.',
  'Remove \\ac{QAQC} columns and remove data based on \\ac{QAQC} flag values for a \\texttt{swmpr} object.  Only applies if \\ac{QAQC} columns are present. ', 
  'View a summary of the number of observations in a \\texttt{swmpr} object that are assigned to different \\ac{QAQC} flags used by \\ac{CDMO}.  The output can be used to inform further processing.', 
  'Remove replicate nutrient data that occur on the same day.  The default is to average replicates.', 
  'Format data from a \\texttt{swmpr} object to a continuous time series at a given timestep.  The function is used in \\texttt{comb} and can also be used with individual stations.',
  'Subset by dates and/or columns for a \\texttt{swmpr} object.  This is a method passed to the generic \\texttt{subset} function provided in the base package.'
  )

to_tab <- data.frame(Functions = funcs, Description = descrips, stringsAsFactors = F)

library(Hmisc)

latex(
  to_tab[, 'Description', drop = F],
  file = '',
  caption = "Organizing functions available from the SWMPr package. Full documentation for each function is in the help file (e.g., execute \\texttt{?comb} for individual functions or \\texttt{help.search(`organize', package = `SWMPr')} for all).",
  rowlabel = 'Function',
  colheads = 'Description',
  rowname = to_tab$Functions,
  caption.loc = 'top',
  col.just=c("p{3.5in}"), 
  label = 'tab:organize', 
  table.env = FALSE
  )

## ----eval = F, cache = T-------------------------------------------------
## # qaqc screen for a swmpr object, retain only '0'
## qaqc(dat)
## 
## # retain all data regardless of flag
## qaqc(dat, qaqc_keep = NULL)
## 
## # retain only '0' and '-1' flags
## qaqc(dat, qaqc_keep = c(0, -1))

## ----eval = F, cache = T-------------------------------------------------
## # view the number of observations in each QAQC flag
## qaqcchk(dat)

## ----eval = F, cache = T-------------------------------------------------
## # get nutrient data
## data(apacpnut)
## dat <- apacpnut
## dat <- qaqc(dat)
## 
## # remove replicate nutrient data
## rem_reps(dat)
## 
## # use different function to aggregate replicates
## func <- function(x) max(x, na.rm = T)
## rem_reps(dat, FUN = func)

## ----eval = F, cache = T-------------------------------------------------
## # import data
## data(apaebmet)
## dat <- apaebmet
## 
## # select two parameters from dat
## subset(dat, select = c('rh', 'bp'))
## 
## # subset records greater than or equal to a date
## subset(dat, subset = '2013-01-01 0:00', operator = '>=')
## 
## # subset records within a date range
## subset(dat, subset = c('2012-07-01 6:00', '2012-08-01 18:15'))
## 
## # subset records within a date range, select two parameters
## subset(dat, subset = c('2012-07-01 6:00', '2012-08-01 18:15'),
##   select = c('atemp', 'totsorad'))
## 
## # remove rows/columns that do not contain data
## subset(dat, rem_rows = T, rem_cols = T)

## ----eval = F, cache = T-------------------------------------------------
## # import, qaqc removal
## data(apadbwq)
## dat <- qaqc(apadbwq)
## 
## # convert time series to two hour invervals
## # tolerance of +/- 30 minutes for matching existing data
## setstep(dat, timestep = 120, differ = 30)
## 
## # convert a nutrient time series to a continuous time series
## # then remove empty rows and columns
## data(apacpnut)
## dat_nut <- apacpnut
## dat_nut <- setstep(dat_nut, timestep = 60)
## subset(dat_nut, rem_rows = T, rem_cols = T)

## ----eval = F, cache = T-------------------------------------------------
## # get nut, wq, and met data as separate objects
## data(apacpnut)
## data(apacpwq)
## data(apaebmet)
## swmp1 <- apacpnut
## swmp2 <- apacpwq
## swmp3 <- apaebmet
## 
## # combine nut and wq data by union
## comb(swmp1, swmp2, method = 'union')
## 
## # combine nut and wq data by intersect
## comb(swmp1, swmp3, method = 'intersect')
## 
## # combine nut, wq, and met data by nut time series, two hour time step
## comb(swmp1, swmp2, swmp3, timestep = 120, method = 'apacpnut')

## ----results = 'asis', echo = FALSE--------------------------------------

funcs <- c('aggreswmp', 'aggremetab', 'ecometab', 'decomp', 'decomp\\_cj', 'hist', 'lines', 'map\\_reserve', 'na.approx', 'plot', 'plot\\_metab', 'plot\\_summary', 'smoother')
funcs <- paste0('\\texttt{', funcs, '}')

descrips <- c(
  'Aggregate \\texttt{swmpr} objects for different time periods - years, quarters, months,  weeks, days, or hours.  The aggregation function defaults to the mean.',
  'Aggregate metabolism data from a \\texttt{swmpr} object.  This is primarly used within \\texttt{plot\\_metab} but may be useful for simple summaries of raw metabolism data.',
  'Estimate ecosystem metabolism for a combined water quality and weather dataset using the open-water method.',
  'Decompose a \\texttt{swmpr} time series into trend, seasonal, and residual components.  This is a simple wrapper to \\texttt{decompose} \\cite{Kendall83}.  Decomposition of monthly or daily trends is possible.',
  'Decompose a \\texttt{swmpr} time series into grandmean, annual, seasonal, and events components.  This is a simple wrapper to \\texttt{decompTs} in the wq package \\cite{Jassby14}.  Only monthly decomposition is possible.',
  'Plot a histogram for a \\texttt{swmpr} object.',
  'Add lines to an existing plot created with \\texttt{plot}.',
  'Create a map of all stations in a reserve using the ggmap package.',
  'Linearly interpolate missing data (\\texttt{NA} values) in a \\texttt{swmpr} object. The maximum gap size that is interpolated is defined by the arguments.', 
  'Plot a univariate  time series for a \\texttt{swmpr} object.  The parameter name must be specified.',
  'Plot ecosystem metabolism estimates after running \\texttt{ecometab} on a \\texttt{swmpr} object.',  
  'Create summary plots of seasonal/annual trends and anomalies for a water a single paramter of interest.',
  'Smooth \\texttt{swmpr} objects with a moving window average.  Window size and sides (e.g., centered) can be specified, passed to \\texttt{filter}.'
  )

to_tab <- data.frame(Functions = funcs, Description = descrips, stringsAsFactors = F)

library(Hmisc)

latex(
  to_tab[, 'Description', drop = F],
  file = '',
  caption = "Analysis functions available from the SWMPr package.  Full documentation for each function is in the help file (e.g., execute \\texttt{?aggreswmp} for individual functions or \\texttt{help.search(`analyze', package = `SWMPr')} for all).",
  rowlabel = 'Function',
  colheads = 'Description',
  rowname = to_tab$Functions,
  caption.loc = 'top',
  col.just=c("p{3.5in}"), 
  label = 'tab:analyze', 
  table.env = FALSE
  )

## ----eval = T, cache = T-------------------------------------------------
# get data, keep all observations
data(apacpnut)
dat <- qaqc(apacpnut, qaqc_keep = NULL)

# aggregate by quarters
agg_dat <- aggreswmp(dat, by = 'quarters')
nrow(agg_dat)

# aggregate by quarters, remove rows with NA values
# note the reduction in the number of rows
agg_dat2 <- aggreswmp(dat, by = 'quarters', na.action = na.omit)
nrow(agg_dat2)

## ----smooth_ex, eval = T, fig.height = 3, cache = T, fig.cap = "Raw and smoothed dissolved oxygen data for a two-week period after using the \\texttt{smoother} function."----
# import data, qaqc and subset
data(apadbwq)
dat <- qaqc(apadbwq)
dat <- subset(dat, select = 'do_mgl', 
  subset = c('2012-07-09 00:00', '2012-07-24 00:00')
  )

# smooth
dat_smooth <- smoother(dat, window = 50, params = 'do_mgl')

# plot raw and smoothed
plot(dat)
lines(dat_smooth, col = 'red', lwd = 2)

## ----interp_ex, eval = T, fig.height = 6, cache = T, fig.scap = "Examples illustrating use of the \\texttt{na.approx} function to fill gaps of different sizes in a dissolved oxygen time series for a four day period.", fig.cap = "Examples illustrating use of the \\texttt{na.approx} function to fill gaps of different sizes in a dissolved oxygen time series for a four day period."----
# get data, qaqc and subset
data(apadbwq)
dat <- qaqc(apadbwq)
dat <- subset(dat, select = 'do_mgl', 
  subset = c('2013-01-22 00:00', '2013-01-26 00:00'))

# interpolate, maxgap of 10 records
fill1 <- na.approx(dat, params = 'do_mgl', maxgap = 10)

# interpolate maxgap of 30 records
fill2 <- na.approx(dat, params = 'do_mgl', maxgap = 30)

# plot for comparison
par(mfrow = c(3, 1))
plot(dat, main = 'Raw')
plot(fill1, col = 'red', main = 'Interpolation - maximum gap of 10 records')
lines(dat)
plot(fill2, col = 'red', main = 'Interpolation - maximum gap of 30 records')
lines(dat)

## ----decomp_ex1, eval = T, fig.height = 6, cache = T, fig.cap = "An additive decomposition of dissolved oxygen into a trend, seasonal, and random component using the \\texttt{decomp} function."----
# get data
data(apadbwq)
swmp1 <- apadbwq

# subset for daily decomposition
dat <- subset(swmp1, subset = c('2013-07-01 00:00', '2013-07-31 00:00'))

# decomposition and plot
test <- decomp(dat, param = 'do_mgl', frequency = 'daily')
plot(test)

## ----decomp_ex2, eval = T, fig.height = 6, warning = F, cache = T, fig.cap = "Additive decomposition of a multi-year chlorophyll time series into the grandmean, annual, seasonal, and events components using the \\texttt{decomp\\_cj} function."----
# get data
data(apacpnut)
dat <- apacpnut
dat <- qaqc(dat, qaqc_keep = NULL)

# decomposition of chl
decomp_cj(dat, param = 'chla_n')

## ----summary_ex, fig.height = 7, fig.width = 13, message = F, cache = T, fig.cap = "Summaries of a multi-year chlorophyll time series using the \\texttt{plot\\_summary} function.  Summaries include monthly distributions (means on top left, quantiles on bottom left), monthly histograms (center), monthly means by year (top right), deviation from monthly means (middle right), and annual trends as deviations from the grand mean (bottom right)"----
## import data
data(apacpnut)
dat <- qaqc(apacpnut)

## plot
plot_summary(dat, param = 'chla_n', years = c(2007, 2013))

## ----metab_ex, eval = TRUE, cache = TRUE, fig.height = 4, fig.width = 8, warning = FALSE, fig.cap = "Monthly means (95\\% confidence) of ecosystem metabolism estimates (net ecosystem metabolism, gross production, and total respiration) for combined water quality and weather data for two years at Apalachicola Bay, Florida."----
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

## ----map_ex, fig.height = 5, message = F, cache = T, fig.cap = "Locations of all sites at the Jacques Cousteau reserve using the \\texttt{map\\_reserve} function."----
# plot the stations at Jacques Cousteau reserve
map_reserve('jac')

## ----results = 'asis', echo = FALSE--------------------------------------

funcs <- c('calckl', 'metab\\_day', 'param\\_names', 'parser', 'swmpr', 'time\\_vec')
funcs <- paste0('\\texttt{', funcs, '}')

descrips <- c(
  'Estimate the reaeration coefficient for air-sea gas exchange.  Used in the \\texttt{ecometab} function.',
  'Identify the metabolic day for each approximate 24 period in an hourly time series.  Used in the \\texttt{ecometab} function.',
  'Returns column names as a list for the parameter types (nutrients, weather, or water quality).  Includes \\ac{QAQC} columns with \\texttt{f\\_} prefix. Used in the data retrieval functions.',
  'Parses HTML returned from \\ac{CDMO} data requests.  Used in the retrieval functions.',
  'Creates a \\texttt{swmpr} object class.  Used in the data retrieval functions.',
  'Converts time vectors to \\texttt{POSIXct} objects with the appropriate time zone for a site.  Used in the data retrieval functions.'
  )

to_tab <- data.frame(Functions = funcs, Description = descrips, stringsAsFactors = F)

library(Hmisc)

latex(
  to_tab[, 'Description', drop = F],
  file = '',
  caption = "Miscellaneous functions available from the SWMPr package.  Most are used within the main functions above but may be useful for customized evaluations of \\ac{SWMP} data.  Full documentation for each function is in the help file (e.g., execute \\texttt{?calckl} at the command line).",
  rowlabel = 'Function',
  colheads = 'Description',
  rowname = to_tab$Functions,
  caption.loc = 'top',
  col.just=c("p{3.5in}"), 
  label = 'tab:misc', 
  table.env = FALSE
  )


## ----metab_plo, eval = T, cache = T, echo = F, fig.cap = 'Aggregated estimates of net metabolism, gross production, and total respiration for two sites at each \\ac{NERRS} reserve.  Values are daily integrated estimates as mean annual values averaged across all years with 95\\% confidence intervals.  Two sites were chosen from each reserve that had the longest available time series. Sites were assigned to regions based on approximate geographic coordinates.', fig.height = 6, fig.width = 8, out.width = '\\textwidth'----
# # packages to use
# library(SWMPr)
# library(httr)
# library(XML)
# library(foreach)
# library(doParallel)
# 
# # names of files on server
# files_s3 <- httr::GET('https://s3.amazonaws.com/swmpalldata/')$content
# files_s3 <- rawToChar(files_s3)
# files_s3 <- htmlTreeParse(files_s3, useInternalNodes = T)
# files_s3 <- xpathSApply(files_s3, '//contents//key', xmlValue)
# files_s3 <- gsub('\\.RData$', '', files_s3) 
# 
# # find only active water quality and weather sites
# # IP address must be registed with CDMO 
# meta <- site_codes()
# sel <- meta$status %in% 'Active' & !grepl('nut$', meta$station_code)
# meta <- meta[sel, ]
# 
# # get wq sites from meta, then filter by those on server
# wq_sites <- grep('wq$', meta$station_code, value = T)
# wq_sites <- files_s3[gsub('\\.RData$', '', files_s3) %in% wq_sites]
# 
# # setup parallel backend for processing
# cl <- makeCluster(8)
# registerDoParallel(cl)
# strt <- Sys.time()
# 
# # process
# metabs <- foreach(wq_site = wq_sites) %dopar% {
#  
#   library(SWMPr)
#   
#   # progress
#   sink('log.txt')
#   cat(wq_site, which(wq_site == wq_sites), 'of', length(wq_sites), '\n')
#   print(Sys.time() - strt)
#   sink()
#   
#   # find corresponding wx station
#   met_site <- substr(wq_site, 1, 3)
#   met_site <- grep(paste0('^', met_site, '.*met'), files_s3, value = T)
#   met_site <- gsub('\\.RData$', '', met_site)
# 
#   # continue if wx data found
#   if(length(met_site) > 0){
#     
#     # if > 1 wx site, pick the one with more obs
#     if(length(met_site) > 1){
# 
#       met_list <- vector('list', length(met_site))
#       names(met_list) <- met_site
#       for(met in met_site){
#       
#         met_tmp <- import_remote(met)
#         met_list[[met]] <- met_tmp
#       
#       }
#     
#       # pick the one with most obs
#       met_most <- which.max(unlist(lapply(met_list, nrow)))
#       met <- met_list[met_most][[1]]
# 
#     # otherwise load only one
#     } else {
#     
#       # met
#       met <- import_remote(met_site)
#     
#     }
#     
#     ##
#     # load the wq file
#     
#     # wq
#     wq <- import_remote(wq_site)
#     
#     ## 
#     # combine, estimate metabolism, reduce data volume
#     dat <- comb(wq, met, method = attr(met, 'station'))
#     dat <- ecometab(dat)
#     dat <- dat[1, ]
#     
#     dat
#     
#   }
#   
# }
# names(metabs) <- wq_sites
# save(metabs, file = 'metabs.RData')

# # take all metab estimates and summarize by year
# # done for two stations at each reserve with the longest time series
# data(metabs)
# 
# metabs <- lapply(metabs, attr, 'metabolism')
# metabs <- do.call('rbind', metabs)
# metabs$site <- gsub('\\.[0-9]*$', '', as.character(row.names(metabs)))
# row.names(metabs) <- 1:nrow(metabs)
# 
# library(dplyr)
# library(tidyr)
# metabs <- select(metabs, date, Pg, Rt, NEM, site) %>% 
#   mutate(date = as.numeric(strftime(date, '%Y'))) %>% 
#   gather('metab', 'value', 2:4) %>% 
#   group_by(date, site, metab) %>% 
#   summarize(value = mean(value, na.rm = T)) %>% 
#   group_by(site, metab) %>% 
#   na.omit %>% 
#   mutate(reserve = substr(site, 1, 3)) %>% 
#   spread(metab, value)
# 
# # select two sites from each reserve, those with longest records
# metabs <- split(metabs, metabs$reserve)
# metabs <- lapply(metabs, 
#   function(x){
#     
#     num_yrs <- aggregate(date ~ site, x, length)
#     num_yrs <- num_yrs[order(num_yrs$date, decreasing = T), ]
#     
#     if(nrow(num_yrs) >1){
#       site_sel <- num_yrs$site[c(1, 2)]
#     } else {
#       site_sel <- num_yrs$site[1]
#     }
#     
#     out <- x[x$site %in% site_sel, ]
#     return(out)
#     
#   }
#   
# )
# metabs <- do.call('rbind', metabs)
# 
# metabs_sum <- metabs
# save(metabs_sum, file = 'data/metabs_sum.RData')

data(metabs_sum)

library(dplyr)
library(tidyr)

# summarize for annual means, get conf ints
metabs_plo <- gather(metabs_sum, 'metab', 'value', 4:6) %>% 
  group_by(site, metab) %>% 
  summarize(
    margs = qt(1 - 0.05/2, length(value) - 1) * sd(value)/sqrt(length(value)),
    value = mean(value)
    ) %>% 
  mutate(
    value = 0.032 * value,
    margs = 0.032 * margs,
    upper = value + margs, 
    lower = value - margs
    )

# get region
data(stat_locs)
names(stat_locs)[names(stat_locs) %in% 'station_code'] <- 'site'
metabs_plo$site <- gsub('wq$', '', metabs_plo$site)
metabs_plo <- left_join(metabs_plo, stat_locs, by = 'site')
region <- rep('Pacific', nrow(metabs_plo))
region[metabs_plo$longitude > -110] <- 'Gulf'
region[metabs_plo$latitude > 28 & metabs_plo$longitude > -84] <- 'Atlantic'

# combine, select, sort
metabs_plo <- select(metabs_plo, -c(latitude, longitude, gmt_off))
metabs_plo$region <- region
metabs_plo <- arrange(metabs_plo, region, site)

# plot prep
to_plo <- data.frame(metabs_plo)

#reassign factors for ranking in plot
sort_val <- to_plo[to_plo$metab == 'Pg', c('region', 'value', 'site')]
sort_val <- order(sort_val$region, sort_val$value)
to_plo$site <- factor(to_plo$site)
to_plo$site <- factor(
  to_plo$site,
  levels = levels(to_plo$site)[sort_val],
  labels = levels(to_plo$site)[sort_val]
  )
new_labs <- c('Mean net metabolism', 'Mean gross production', 'Mean total respiration')
to_plo$metab <- factor(to_plo$metab, levels = c('NEM', 'Pg', 'Rt'), labels = new_labs)

ylabs <- expression(paste('g ', O[2],' ', m^-2,' ', d^-1))
library(ggplot2)
p <- ggplot(to_plo, 
    aes(x = site, y = value, group = metab, fill = metab)
  ) + 
  geom_bar(stat = 'identity', width = 1, colour = 'black') + 
  scale_y_continuous(ylabs) + 
  facet_grid(metab ~ region, scales = 'free') + 
  scale_x_discrete(name = element_blank()) + 
  theme_bw() + 
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, 
      size = 7),
    legend.position = 'none'
    ) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0)
print(p)


## ----eval = T, results = 'asis', cache = T, echo = F---------------------

data(metabs_sum)

library(dplyr)
library(tidyr)

# summarize by time periods
metabs_tab <- gather(metabs_sum, 'metab', 'value', 4:6)
metabs_tab$yr_grp <- '1995-2004'
metabs_tab$yr_grp[metabs_tab$date > 2004] <- '2005-2014'
metabs_tab <- mutate(metabs_tab, value = 0.032 * value) %>% 
  group_by(site, metab, yr_grp) %>% 
  arrange(yr_grp, date, metab) %>% 
  summarize(ave = round(mean(value, na.rm = TRUE), 2)) %>% 
  spread(yr_grp, ave) 
trends <- apply(metabs_tab, 1, function(x){
  
  val <- as.numeric(x[['2005-2014']])
  trnd <- val > as.numeric(x[['1995-2004']])
  
  if(is.na(trnd)) return(val)
  if(trnd) return(paste0('{\\bf ', val, '}'))
  if(!trnd) return(paste0('{\\it ', val, '}'))
  
})
metabs_tab[, '2005-2014'] <- trends
metabs_tab <- gather(metabs_tab, 'period', 'value', 3:4) %>% 
  unite('var', metab, period, sep = ' ') %>% 
  spread(var, value) %>% 
  mutate(site = gsub('wq$', '', site))

# get region
data(stat_locs)
names(stat_locs)[names(stat_locs) %in% 'station_code'] <- 'site'
metabs_tab <- left_join(metabs_tab, stat_locs, by = 'site')
region <- rep('Pacific', nrow(metabs_tab))
region[metabs_tab$longitude > -110] <- 'Gulf'
region[metabs_tab$latitude > 28 & metabs_tab$longitude > -84] <- 'Atlantic'

# combine, select, sort
metabs_tab <- select(metabs_tab, -c(latitude, longitude, gmt_off)) %>% 
  mutate(region = region) %>% 
  arrange(region, site)

nrgroups <- as.numeric(table(metabs_tab$region))
# table

cap <- 'Trends in metabolism for two sites at each of the \\ac{NERRS} reserves.  Values are averages of mean annual estimates for each period of observation (1994-2004 and 2005-2014). Bold values indicate an increase from the first period, whereas italic values indicate a decrease. Sites were assigned to regions based on approximate geographic coordinates.' 
foots <- '\\footnotesize{\\textsuperscript{{\\it a}}NEM: net ecosystem metabolism, Pg: gross production, Rt: total respiration, all values in g O$_2$ m$^{-2}$ d$^{-1}$ as annual averages.}' 

library(Hmisc)
latex(
  metabs_tab[, !names(metabs_tab) %in% c('site', 'region')], 
  file = '', 
  caption = cap, 
  cgroup = c('NEM\\textsuperscript{{\\it a}}', 'Pg', 'Rt'), 
  n.cgroup = c(2, 2, 2), 
  colheads = rep(c('1995-2004', '2005-2014'), 3),
  rowlabel = 'Site', 
  rowname = metabs_tab$site,
  digits = 2,
  insert.bottom = foots,
  rgroup = unique(metabs_tab$region),
  n.rgroup = nrgroups,
  label = 'tab:metab_tab',
  size = 'footnotesize'
  )
  

