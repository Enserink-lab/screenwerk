#' Saving the consolidated data to a file
#'
#' @description
#' The function \emph{\code{save}} exports the consolidated data to a file format, which can be used for downstream analysis with other tools.
#' 
#' @param x an object of class 'consolidatedData'.
#' @param .saveto \code{character}; a path to a folder location where the object is saved to.
#' @param .fileformat \code{character}; a string or list determining the file format of the raw data files.
#' @param .format \code{character}; a predefined identifier, which will decide the format and layout of the data to be exported (see Details, for more information).
#' @param .sep \code{character}; a string or list determining the field separator character.
#' @param ... further arguments passed to or from other methods.
#' 
#' @details This function extends the generic function \code{\link{save}} for objects of class 'consolidatedData'. It exports the consolidated data file for use outside the scope of screenwerk. 
#' 
#' With the argument \code{.format} the data can be exported for specific analytical applications in the required format. For now, the list of supported formats include: BREEZE
#' The file type can be set with the argument \code{.fileformat}, which will set the file extension accordingly.
#' In addition, the field separator (deliminator) can be set with the argument \code{.sep}.
#' 
#' @seealso \code{\link{consolidateData}}
#'
#' @examples
#' \donttest{\dontrun{
#' # Saving the consolidated data to a csv file
#' save(consolidatedData, .saveto = "../myExportFiles/files", .fileformat = ".csv", .sep = ";")
#' 
#' # Exporting the consolidated data for BREEZE
#' save(consolidatedData, .saveto = "../myExportFiles/files", .fileformat = ".csv", .format = "breeze")
#' }}
#'
#' @keywords drug screen drug combination dispensing save
#' 
#' 
#' @importFrom utils write.csv
#' 
#' 
#' @export

save.consolidatedData <- function(x, .saveto, ..., .fileformat, .format = FALSE, .sep){
  
  # Check, if the dispensing data has been provided as an object of class S3:consolidatedData
  if(missing(x)){stop("Consolidated data missing! Please provide a consolidated data set.", call. = TRUE)}
  if(class(x) != "consolidatedData"){stop("Data not of class 'consolidatedData'!", call. = TRUE)}
  
  # Check, if a folder location has been provided
  if(missing(.saveto)){ message("Folder location not specified. Saving to default working directory:\n", getwd()) }
  if(missing(.saveto)){ .saveto <- getwd() }
  # Create folder, if it does not exist
  if(!file.exists(file.path(.saveto))){ dir.create(file.path(.saveto), showWarnings = FALSE, recursive = TRUE) }
  if(!utils::file_test("-d", .saveto)){stop("in '.saveto'. Argument needs to be a valid folder location.", call. = TRUE)}
  
  
  # Check, if a field separator was defined
  if(missing(.sep)){ 
    message("Field separator not sepecified. Using default: 'comma (,)'")
    .sep = "," } 

  # Check, if the format have been specified
  if(missing(.fileformat)){
    message("File format not sepecified. Using default: '.csv'")
    .format <- ".csv"
  } else {
    # Clean up file format
    .fileformat <- sub("\\.", "", .fileformat)
  }
  
  # Check, if the format have been specified
  if(any(missing(.format), is.null(.format))){ .format <- FALSE 
  } else if(all(.format != "BREEZE", !grepl("breeze", .format, ignore.case = TRUE)) &
            all(.format != "none", !grepl("no", .format, ignore.case = TRUE))){
    stop("in '.format'. Data format not supported! Please select a valid data format. Supported formats: 'breeze'. \nFor more information and help, see 'R Documentation'.", call. = TRUE)
  }
  
  if(any(.format == "BREEZE", grepl("breeze", .format, ignore.case = TRUE)) & .sep == ";"){
    message("Field separator not compatible with BREEZE format. Using default: 'comma (,)'")
  }
  
  
  
  # Assign object to internal variable
  consolidatedData = x
  
  
  cat('\n', "> Reading consolidated data...", sep = "")
  
  # Extract the consolidated data
  .consolidatedData <- consolidatedData[["consolidated"]]
  rownames(.consolidatedData) <- NULL
  
  
  # Check for export format
  if(any(.format == "breeze", grepl("breeze", .format, ignore.case = TRUE))){
    
    # Export data to breeze format
    dfs.breeze <- .consolidatedData[c("Destination.Well", "Plate.Number", "Drug", "Drug.Concentration", "Sample", "CPS")]
    dfs.breeze <- dfs.breeze[with(dfs.breeze, order(`Sample`, `Plate.Number`,  Drug)),]
    names(dfs.breeze) <- c("WELL", "PLATE", "DRUG_NAME", "CONCENTRATION", "SCREEN_NAME", "WELL_SIGNAL")
    
    cat('\r', "> Saving consolidated file for BREEZE... ", strrep(" ", 100), sep = "")
    
    utils::write.table(dfs.breeze, file = file.path(.saveto, paste("BREEZE data file", ".", .fileformat, sep = "")), sep = ",", dec = ".", row.names = FALSE, quote = FALSE)
    
  } else {
    
    cat('\r', "> Saving consolidated file... ", strrep(" ", 100), sep = "")
    
    utils::write.table(.consolidatedData, file = file.path(.saveto, paste("Consolidated data file", ".", .fileformat, sep = "")), row.names = FALSE, quote = FALSE, sep = .sep)
    
  }
  
  cat('\r', "Finished saving consolidated file for ", dispensingData$dispensingID, ".", strrep(" ", 100), '\n', sep = "")
  
}
