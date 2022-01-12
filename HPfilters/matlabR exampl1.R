## MatlabR
library(matconv)
library(matlabr)
options(matlab.path = "/Applications/Polyspace/Matlab.app/bin/maci64")

have_matlab()
get_matlab()

library(mFilter)
library(xts)
library(tsbox) # Convert time series types


code = c("x = 10", 
         "y=20;", 
         "z=x+y", 
         "a = [1 2 3; 4 5 6; 7 8 10]",
         "save('test.txt', 'x', 'a', 'z', '-ascii')")
res = run_matlab_code(code)

### Code 2 ----

#set a variable in R and save in a csv file
x <- 10
write.table(x, file='~/x.csv', sep=",", row.names=FALSE, col.names=FALSE)

#make a vector where each element is a line of MATLAB code
#matlab code reads in our variable x, creates two variables y and z, 
#and write z in a csv file
matlab.lines <- c(
  "x = csvread('~/x.csv')",
  "y=20",
  "z=x+y",
  "csvwrite('~/z.csv', z)")

#create a MATLAB script containing all the commands in matlab.lines
writeLines(matlab.lines, con="~/myscript.m")

#run our MATLAB script
system("matlab -nodisplay -r \"run('~/myscript.m'); exit\"")

# read in the variable z we created in MATLAB
z <- read.table("~/z.csv")
z

#remove the temporary files we used to pass data between R and MATLAB
system("rm ~/x.csv ~/z.csv")

# Code 3: example one sided hp filter ----
## Load data
## Set working directory
library('rstudioapi')
setwd(dirname(getActiveDocumentContext()$path))
getwd()

country="US"
enddate ='2021-01-01'
filepath1 = ('../../HPCredit/Data Collection/1.Latest/Paper2')
filepath2 = sprintf('/Credit_BISgap_%s.txt',country)
filepath = paste(filepath1, filepath2, sep='')

df <- read.table(filepath, header=TRUE, sep=',')
# df = subset(df, date > as.Date(startdate)) # Limit series data to after 1990
df <- df %>%
  filter(grepl("^(A)", cg_dtype))
df = subset(df, date <= as.Date(enddate)) # Limit series data to before 2021
varlist = c("date", "obs_value")
df = df[varlist]
credit_BIS <- xts(df[,-c(1)], order.by=as.Date(df[,"date"], "%Y-%m-%d"))
credit_BIS <- ts_ts(credit_BIS)

## set a variable in R and save in a csv file
x <- credit_BIS
write.table(x, file='~/x.csv', sep=",", row.names=FALSE, col.names=FALSE)

## Parse with Matlab code
  matlab.lines <- c(
  "working_dir = ['/Users/namnguyen/Documents/GitHub/HPI-Credit-Trasitory-Forecast/HPfilters'];",
  "cd(working_dir);",
  "Y= csvread('~/x.csv')",
  "source('one_sided_hp_filter_kalman.m');",
  "[Ytrend, Yobs] = one_sided_hp_filter_serial(Y));",
  "csvwrite(['~/z.csv'],[Ytrend,Yobs]);")

  res = run_matlab_code(matlab.lines)
## Import back into R
  
system("rm ~/x.csv ~/z.csv")