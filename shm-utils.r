# Copyright (c) 2015 National ICT Australia Limited (NICTA).
# Thierry.Rakotoarivelo@nicta.com.au
#
# Useful R functions for SHM Analysis
#

# About installing the libraries:
# - need on ubuntu: mesa-common-dev , libglu1-mesa-dev , libpq-dev
#     update.packages(checkBuilt = TRUE, ask = FALSE)
#     p = c(Hmisc','RSQLite','TSclust','rwt')
#     install.packages(pkgs=p,dependencies=TRUE)
library('Hmisc')
library("RSQLite")
library("TSclust")
library('rwt')

a_day_in_sec = 60*60*24

# Returns Unix-time timestamps, i.e. number of second since Unix-time
#
day <- function(d,m,y,hour=0,minute=0) {
  s <- sprintf("%s-%s-%s %02d%02d", d,m,y, hour,minute)
  return(as.integer(as.POSIXct(s, format="%d-%m-%Y %H%M")))
}

# Returns Unix-time timestamps in us, i.e. number of microsecond since Unix-time
#
day.us <- function(d,m,y,hour=0,minute=0) {
  s <- sprintf("%s-%s-%s %02d%02d", d,m,y, hour,minute)
  return(as.integer(as.POSIXct(s, format="%d-%m-%Y %H%M"))*1000000)
}

# Convert a number of second or a number of microsecond to a Time
#
epoch.convert <- function(time) {
  # If the input is greater than 10^15, then assume it is in microsecond
  if (time>1000000000000000) { as.POSIXct(time/1000000, origin="1970-01-01") }
  else { as.POSIXct(time, origin="1970-01-01") }
}

# Format a number to the specified amount of digit
#
specify_decimal <- function(x, k) format(round(x, k), nsmall=k)

# Perform 1D Wavelet Denoising with bins of 8*512 samples
#
wavelet.denoise.8 <- function(s, thres=4, quiet=TRUE, autopadding=TRUE) {
  o = default.dwt.option
  o$threshold = thres
  h <- daubcqf(12)$h.0
  w = 8*512 # IMPORTANT: length(s) must be > w!
  k = floor(length(s)/w)
  start = 1
  end = w
  res = c()
  # When 'autopadding' and length(s) is not a multiple of w, add some paddings
  if ((autopadding) & (isTRUE(all.equal(r, as.integer(r)))==FALSE)) {
    k = ceiling(r)
    s = c(s,rep(mean(s),((k*w)-length(s))))
  }
  # When 'quiet', start Redirect stdout/stderr to dev/null
  if (quiet) {
    zz <- file("/dev/null", open = "wt")
    sink(zz)
    sink(zz, type = "message")
  }
  for (i in c(1:k)) {
    ss = s[start:end]
    print(paste(i,length(ss)))
    ret.dwt = denoise.dwt(ss, h,o)
    xd <- ret.dwt$xd
    res = c(res,xd)
    start = end + 1
    end = start + w -1
  }
  if (quiet) {
    # Stop Redirect stdout/stderr to dev/null
    sink(type = "message")
    sink()
    close(zz)
  }
  # When 'autopadding' and length(s) is not a multiple of w, remove the paddings
  if ((autopadding) & (isTRUE(all.equal(r, as.integer(r)))==FALSE)) {
    res=res[1:original_length]
  }
  return(res)
}

# For backward compatibility with older analysis scripts
#
wavelet.denoise.8.autopadding <- function(s, thres=4) {
  return(wavelet.denoise.8(s=s,thres=thres,autopadding=TRUE))
}

# Perform 1D Wavelet Denoising on a single event of 500 samples
#
wd.event <- function(s, thres=4,quiet=TRUE) {
  library('rwt')
  o = default.dwt.option
  o$threshold = thres
  h <- daubcqf(12)$h.0
  s_padded = c(s,rep(mean(s),12))
  # When 'quiet', start Redirect stdout/stderr to dev/null
  if (quiet) {
    zz <- file("/dev/null", open = "wt")
    sink(zz)
    sink(zz, type = "message")
  }
  s_denoised = as.vector((denoise.dwt(s_padded, h,o))$xd)
  if (quiet) {
    # Stop Redirect stdout/stderr to dev/null
    sink(type = "message")
    sink()
    close(zz)
  }
  return(s_denoised[1:500])
}
