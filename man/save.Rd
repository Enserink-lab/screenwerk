% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/save.R
\name{save}
\alias{save}
\title{Generic save function}
\usage{
save(x, .saveto, .fileformat, .sep, .sets, .labels, .split, .format)
}
\arguments{
\item{x}{an object of class 'consolidatedData'.}

\item{.saveto}{\code{character}; a path to a folder location where the object is saved to.}

\item{.fileformat}{\code{character}; a string or list determining the file format of the raw data files.}

\item{.sep}{\code{character}; a string or list determining the field separator character.}

\item{.sets}{\code{numeric}; a value determining the number of sets to create, or logical; if TRUE, the number of sets will depend on the labels provided.}

\item{.labels}{\code{vector}; a set of labels to be used for each individual set.}

\item{.split}{\code{logical}; if TRUE, the function will split each set into individual files. The default is FALSE, in which a single file is created.}

\item{.format}{\code{character}; a predefined identifier, which will decide  the format and layout of the dispensing file (see Details, for more information).}
}
\description{
The function \emph{\code{save}} exports data sets to files.
}
\details{
This function provids the generic function \code{\link{save}} for objects of other class.
}
