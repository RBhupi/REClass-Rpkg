% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/wt-class.R
\name{getScaleBreak}
\alias{getScaleBreak}
\title{compute scale break for convection and stratiform regions.}
\usage{
getScaleBreak(res_km, conv_scale_km)
}
\arguments{
\item{res_km}{resolution of the image.}

\item{conv_scale_km}{expected size of spatial variations due to convection.}
}
\value{
dyadic integer scale break in pixels.
}
\description{
WT will be computed upto this scale and features will be designated as convection.
}
