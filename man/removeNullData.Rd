% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/removeNullData.R
\name{removeNullData}
\alias{removeNullData}
\title{Function removeNullData}
\usage{
removeNullData(user_DataObj, write_log)
}
\arguments{
\item{user_DataObj}{Class DataObj, to be stripped. Must be modifiable.}

\item{write_log}{Function to write output to a log file.}
}
\description{
Function to strip out sgRNA & samples with no count data.
Called by default within DataObjClass.  Do not specify test_tumor_subtype_cols
in DataObjClass
until AFTER erroneous samples are removed, or indicing will be messed up.
}
