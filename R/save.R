#' Generic save function
#'
#' @description
#' The function \emph{\code{save}} exports data sets to files.
#' 
#' @param x an object of class 'consolidatedData'.
#' @param .saveto \code{character}; a path to a folder location where the object is saved to.
#' @param .fileformat \code{character}; a string or list determining the file format of the raw data files.
#' @param .sep \code{character}; a string or list determining the field separator character.
#' @param .sets \code{numeric}; a value determining the number of sets to create, or logical; if TRUE, the number of sets will depend on the labels provided.
#' @param .labels \code{vector}; a set of labels to be used for each individual set.
#' @param .split \code{logical}; if TRUE, the function will split each set into individual files. The default is FALSE, in which a single file is created.
#' @param .format \code{character}; a predefined identifier, which will decide  the format and layout of the dispensing file (see Details, for more information).
#' 
#' @details This function provids the generic function \code{\link{save}} for objects of other class.
#' 
#' 
#' @export
save <- function(x, .saveto, .fileformat, .sep, .sets, .labels, .split, .format) UseMethod("save")
