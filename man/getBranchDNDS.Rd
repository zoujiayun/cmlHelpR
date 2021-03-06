% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/getBranchdNdS.R
\name{getBranchDNDS}
\alias{getBranchDNDS}
\title{Parse dN/dS branch values from output files}
\usage{
getBranchDNDS(dir_path, models, ext = ".out")
}
\arguments{
\item{dir_path}{Path to parallel_codeml output directory}

\item{models}{Models to get branch dN/dS values for}

\item{ext}{Extension of CODEML output files}
}
\description{
Given the path to the parallel_codeml output directory, this function will parse the branch dN/dS values from the
CODEML output. It returns a list and nested dataframe object of dN/dS values for each file.
}
\keyword{helper}
