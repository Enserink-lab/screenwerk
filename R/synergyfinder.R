#' Calculating synergies with synergyfinder
#' 
#' @description
#' An essential component of the modular library \pkg{screenwerk}, and imperative for the analysis of the experimental data from a drug combination screen.
#' \emph{\code{synergyfinder}} is a function that calculates synergies between two drugs in a drug combination screen based on either one of the following models: Loewe, Bliss, HSA, or ZIP
#' 
#' @param processedData an object of class 'processedData'.
#' @param synergymodel \code{character}; a predefined identifier, representing the synergy model. 
#' @param .saveoutput \code{logical}; if TRUE, then the output from code{bayesynergy} will be saved to file. The default is FALSE, in which the output is not saved.
#' @param .plot \code{logical}; if TRUE, then plots from code{bayesynergy} will be saved to file. The default is FALSE, in which no plots will be generated.
#' @param .saveto \code{string}; path to a folder location where the results are saved to.
#' 
#' 
#' @details The function \code{synergyfinder} is used to assess the interaction between two drugs across their dose ranges.
#' The analysis is based on the package \code{vignette("synergyfinder", package = "synergyfinder")} as described in the original paper by Zhenga, 2022, Genomics, Proteomics & Bioinformatics (see references).
#' 
#' 
#' Files are saved either to the specified location or the default working environment, with the corresponding folder structure: 'results/data/synergyfinder' and 'results/graphs/synergyfinder'
#' 
#' 
#' @return Returns an class S3 object of type list with the synergy scores for each drug pair.
#' 
#' @references Zheng S, Wang W, Aldahdooh J, Malyutina A, Shadbahr T, Tanoli Z, Pessia A, Tang J (2022). SynergyFinder Plus: Toward Better Interpretation and Annotation of Drug Combination Screening Datasets. Genomics, Proteomics & Bioinformatics, DOI: 10.1016/j.gpb.2022.01.004
#' 
#' @seealso \code{vignette("synergyfinder", package = "synergyfinder")}
#' 
#' @examples
#' \donttest{\dontrun{
#' # Run bayesnergy and save both the output and plots
#' synergyfinder(processedData, synergymodel = "ZIP", .saveoutput = TRUE, .plot = TRUE, .saveto = "path/to/folder/")
#' 
#' # Run bayesynergy without generating plots or saving the output
#' synergyfinder(processedData, synergymodel = "ZIP")
#' }}
#'
#' @keywords drug screen analysis dose response curve matrix
#' 
#' @importFrom synergyfinder ReshapeData CalculateSynergy PlotDoseResponseCurve Plot2DrugHeatmap Plot2DrugContour Plot2DrugSurface PlotDoseResponse
#' 
#' @export

synergyfinder <- function(processedData, synergymodel = c("ZIP", "HSA", "Bliss", "Loewe"), .saveoutput, .plot, .saveto){
  
  # Check, if synergyfinder is installed
  if (!require("synergyfinder")) { stop("Dependency missing! R-package 'synergyfinder' not installed.", call. = TRUE) }
  
  # Check, if the data has been provided as an object of class S3:processedData
  if(missing(processedData)){stop("Data missing! Please provide a data set as an object of class 'processedData' ", call. = TRUE)}
  if(class(processedData) != "processedData"){
    stop("Provided data not of class 'processedData'!", call. = TRUE)
  } else {
    # Extract the required data from the class S3 object
    analysisData <- processedData[["analysisData"]]
    synDataset <- processedData[["splitDataset"]]
    synDRM <- processedData[["doserespMatrix"]]
  }
  
  # Check, which synergymodel should be used
  if(missing(synergymodel)){ synergymodel = "ZIP" } 
  else if(!synergymodel %in% c("ZIP", "HSA", "Bliss", "Loewe")){
    stop("in 'synergymodel'. Argument needs to be one of the following models: 'ZIP', 'HSA', 'Bliss', 'Loewe'", call. = TRUE)
  }
  
  synergymodel <- match.arg(synergymodel)
  
  
  # Check, if synergyfinder output should be saved
  if(missing(.saveoutput)){ .saveoutput <- FALSE } 
  else if(!is.logical(.saveoutput)){
    stop("in '.saveoutput'. Argument needs to be a logical value, either TRUE or FALSE.", call. = TRUE)
  }
  
  # Check, if synergyfinder plots should be generated
  if(missing(.plot)){ .plot <- FALSE } 
  else if(!is.logical(.plot)){
    stop("in '.plot'. Argument needs to be a logical value, either TRUE or FALSE.", call. = TRUE)
  }
  
  # Check, if a folder location has been provided
  if(missing(.saveto)){ .saveto <- getwd() }
  # Create folder, if it does not exist
  if(!file.exists(file.path(.saveto))){ dir.create(file.path(.saveto), showWarnings = FALSE, recursive = TRUE) }
  if(!utils::file_test("-d", .saveto)){ stop("in '.saveto'. Argument needs to be a valid folder location.", call. = TRUE) }
  
  
  
  cat("Running synergy analysis with synergyfinder ", paste(packageVersion("synergyfinder")), ".", '\n\n', sep = "")
  
  
  # SECTION A: REFORMAT DATA SET TO BAYESYNERGY REQUIREMENTS  ##############################
  # Merge individual data sets before dividing data sets based on individual drugs
  
  synDFS = list()
  
  for (samplename in names(synDataset)){
    
    cat("Assembling data set for ", samplename, ".", sep = "")
    
    
    # Split the synergy data for each sample (cell lines) to create a dose response matrix (drm)
    synDFS[[samplename]][["drm"]] <- split(synDataset[[samplename]][["combinationData"]], synDataset[[samplename]][["combinationData"]]$Drug)
    
    
    for(drugname in names(synDFS[[samplename]][["drm"]])){
      
      # Add the corresponding drug pair to each drug
      synDFS[[samplename]][["drm"]][[drugname]] <- merge(synDFS[[samplename]][["drm"]][[drugname]], synDataset[[samplename]][["combinationData"]], 
                                                         by = c("Sample", "Destination.Plate.Barcode", "Plate.Number", "Combination.ID",
                                                                "Destination.Well", "Viability (%)"), suffixes = c(".1", ".2"))      
      # Remove same drug pairs and select specific columns
      synDFS[[samplename]][["drm"]][[drugname]] <- subset(synDFS[[samplename]][["drm"]][[drugname]], Drug.1 != Drug.2)
      synDFS[[samplename]][["drm"]][[drugname]] <- synDFS[[samplename]][["drm"]][[drugname]][c("Drug.1", "Drug.2", "Drug.Concentration.1", "Drug.Concentration.2", "Unit.1", "Unit.2", "Viability (%)")]
      
      # Rename selected columns
      names(synDFS[[samplename]][["drm"]][[drugname]])[names(synDFS[[samplename]][["drm"]][[drugname]]) %in%  c("Drug.1", "Drug.2", "Drug.Concentration.1", "Drug.Concentration.2", "Unit.1", "Unit.2", "Viability (%)")] <- c("drug_col", "drug_row", "conc_c", "conc_r", "conc_c_unit", "conc_r_unit", "response")
      
      # Reorder columns
      synDFS[[samplename]][["drm"]][[drugname]] <- synDFS[[samplename]][["drm"]][[drugname]][,c("response", "drug_row", "drug_col", "conc_r", "conc_c", "conc_r_unit", "conc_c_unit")]
      
      # Split the dose response matrix (drm) for each drugpair
      # synDFS[[samplename]][["drm"]][[drugname]] <- split(synDFS[[samplename]][["drm"]][[drugname]], synDFS[[samplename]][["drm"]][[drugname]]$drug_row)
      
      
      
      # Add to each dose belonging to the drugpair a corresponding dose zero to the reference drug. The addition of a dose zero is avoided to only drugpairs that are the same as the reference drug. 
      # Note: A zero dose is added to each of the 6 doses of the 64 drugs w/o the reference drug
      synDFS[[samplename]][["drm"]][[drugname]] <- rbind(synDFS[[samplename]][["drm"]][[drugname]], data.frame(response = synDataset[[samplename]][["singleDrugResponseData"]]$`Viability (%)`[which(synDataset[[samplename]][["singleDrugResponseData"]]$Drug != drugname)],
                                                                                                               drug_row = synDataset[[samplename]][["singleDrugResponseData"]]$Drug[which(synDataset[[samplename]][["singleDrugResponseData"]]$Drug != drugname)],
                                                                                                               drug_col = drugname,
                                                                                                               conc_r = synDataset[[samplename]][["singleDrugResponseData"]]$Drug.Concentration[which(synDataset[[samplename]][["singleDrugResponseData"]]$Drug != drugname)],
                                                                                                               conc_c = 0,
                                                                                                               conc_r_unit = synDataset[[samplename]][["singleDrugResponseData"]]$Unit[which(synDataset[[samplename]][["singleDrugResponseData"]]$Drug != drugname)],
                                                                                                               conc_c_unit = synDataset[[samplename]][["singleDrugResponseData"]]$Unit[which(synDataset[[samplename]][["singleDrugResponseData"]]$Drug == drugname)], check.names = FALSE, stringsAsFactors = FALSE))
      
      # As above, howevere now a dose zero is added to each dose fo the reference drug. Add to each dose belonging to the reference drug a corresponding dose zero to the durgpair.
      # Note: A zero dose is added to each of the 6 doses from the reference drug
      synDFS[[samplename]][["drm"]][[drugname]] <- rbind(synDFS[[samplename]][["drm"]][[drugname]], data.frame(response = synDataset[[samplename]][["singleDrugResponseData"]]$`Viability (%)`[which(synDataset[[samplename]][["singleDrugResponseData"]]$Drug == drugname)],
                                                                                                               # Important to sort or group the different drugs so that the response matches the order of each drug group
                                                                                                               drug_row = sort(synDataset[[samplename]][["singleDrugResponseData"]]$Drug[which(synDataset[[samplename]][["singleDrugResponseData"]]$Drug != drugname)]),
                                                                                                               drug_col = drugname,
                                                                                                               conc_r = 0,
                                                                                                               conc_c = synDataset[[samplename]][["singleDrugResponseData"]]$`Drug.Concentration`[which(synDataset[[samplename]][["singleDrugResponseData"]]$Drug == drugname)],
                                                                                                               conc_r_unit = synDataset[[samplename]][["singleDrugResponseData"]]$Unit[which(synDataset[[samplename]][["singleDrugResponseData"]]$Drug != drugname)],
                                                                                                               conc_c_unit = synDataset[[samplename]][["singleDrugResponseData"]]$Unit[which(synDataset[[samplename]][["singleDrugResponseData"]]$Drug == drugname)], check.names = FALSE, stringsAsFactors = FALSE))
      
      
      # Add a zero dose pair for each drugpair ~ reference drug combination.
      # Note: A single zero dose ~ zero dose pair is added 
      synDFS[[samplename]][["drm"]][[drugname]] <- rbind(synDFS[[samplename]][["drm"]][[drugname]], data.frame(response = 1, # Needs optimization: either randomly picking a DMSO, or using the median of all DMSO measurements from all plates
                                                                                                               drug_row = unique(synDataset[[samplename]][["singleDrugResponseData"]]$Drug[which(synDataset[[samplename]][["singleDrugResponseData"]]$Drug != drugname)]),
                                                                                                               drug_col = drugname,
                                                                                                               conc_r = 0,
                                                                                                               conc_c = 0,
                                                                                                               conc_r_unit = unique(synDataset[[samplename]][["singleDrugResponseData"]]$Unit[which(synDataset[[samplename]][["singleDrugResponseData"]]$Drug != drugname)]),
                                                                                                               conc_c_unit = unique(synDataset[[samplename]][["singleDrugResponseData"]]$Unit[which(synDataset[[samplename]][["singleDrugResponseData"]]$Drug == drugname)]), check.names = FALSE, stringsAsFactors = FALSE))
      
      
      
      # Select only doses that have been combined for each drug and drug pair, respectively
      synDFS[[samplename]][["drm"]][[drugname]] <- do.call(rbind, setNames( lapply(split(synDFS[[samplename]][["drm"]][[drugname]], synDFS[[samplename]][["drm"]][[drugname]]$drug_col), function (a) subset(a, (conc_c %in% c(unique(synDataset$MeWo$combinationData$Drug.Concentration[which(synDataset[[samplename]][["combinationData"]]$Drug == unique(a$drug_col) )]), 0)))), NULL))
      synDFS[[samplename]][["drm"]][[drugname]] <- do.call(rbind, setNames( lapply(split(synDFS[[samplename]][["drm"]][[drugname]], synDFS[[samplename]][["drm"]][[drugname]]$drug_row), function (a) subset(a, (conc_r %in% c(unique(synDataset$MeWo$combinationData$Drug.Concentration[which(synDataset[[samplename]][["combinationData"]]$Drug == unique(a$drug_row) )]), 0)))), NULL))
      
      # Re-arrange data set
      synDFS[[samplename]][["drm"]][[drugname]] <- synDFS[[samplename]][["drm"]][[drugname]][with(synDFS[[samplename]][["drm"]][[drugname]], order(drug_row, drug_col, -conc_r, -conc_c)),]
      
      
      # Create indices for drug, dose and replicate
      synDFS[[samplename]][["drm"]][[drugname]]$block_id <- with(synDFS[[samplename]][["drm"]][[drugname]], ave(drug_row, FUN = function(x) as.integer(factor(x, levels = unique(x)))))
      synDFS[[samplename]][["drm"]][[drugname]]$replicate <- with(synDFS[[samplename]][["drm"]][[drugname]], ave(drug_col, drug_row, conc_c, conc_r, FUN = length))
      
      # Calculate the median for the zero-dose response from the triplicates of the single drug treatments
      synDFS[[samplename]][["drm"]][[drugname]] <- within(synDFS[[samplename]][["drm"]][[drugname]], { response = ifelse(conc_c == 0, ave(response, list(drug_row, conc_r), FUN=median), response) })
      synDFS[[samplename]][["drm"]][[drugname]] <- within(synDFS[[samplename]][["drm"]][[drugname]], { response = ifelse(conc_r == 0, ave(response, list(drug_row, conc_c), FUN=median), response) })
      
      # Remove duplicate rows of triplicates
      synDFS[[samplename]][["drm"]][[drugname]] <- synDFS[[samplename]][["drm"]][[drugname]][!duplicated(synDFS[[samplename]][["drm"]][[drugname]]), ]
      
      
      # Reorder, and select columns
      synDFS[[samplename]][["drm"]][[drugname]] <- synDFS[[samplename]][["drm"]][[drugname]][,c("block_id", "drug_row", "drug_col", "conc_r", "conc_c", "response", "conc_r_unit", "conc_c_unit")]
      
      
      
      
      # Reshaping dose response matrix for synergfinder analysis
      cat('\r', "> Reshaping dose response data for: ", samplename, " : ",  drugname, strrep(" ", 50), sep = "")
      
      suppressMessages(
        synDFS[[samplename]][["ReshapeData"]][[drugname]] <- synergyfinder::ReshapeData(synDFS[[samplename]][["drm"]][[drugname]], data_type = "viability", impute = TRUE, noise = TRUE)
      )
      
      
      
    }
    
    cat('\r', "Finished reshaping data for: ", samplename, strrep(" ", 100), "\n", sep = "")
    
    
    rm(samplename, drugname)
    
  }
  
  
  
  
  # SECTION B: Calculating the synergy score fore each drug combination ##############################
  # Exporting a data set with synergy scores for each cell line and drug combination
  
  for(samplename in names(synDFS)){
    
    cat("Analyzing synergies for ", samplename, ".", "\n", sep = "")
    
    # Calculating individual synergy scores
    for (drugname in names(synDFS[[samplename]][["ReshapeData"]])){
      
      # Generating synergy scores for each drug combination
      cat('\r', "> Calculating synergy score using the ", synergymodel, " model for ", samplename, ": ", drugname, strrep(" ", 100), sep = "")
      
      suppressMessages(
        synDFS[[samplename]][["synergyfinder"]][[drugname]][[synergymodel]] <- synergyfinder::CalculateSynergy(data = synDFS[[samplename]][["ReshapeData"]][[drugname]], method = synergymodel, correct_baseline = "all")
      )
      
      
      # Save  output to file
      if(isTRUE(.saveoutput)){
        if(!file.exists(file.path(.saveto, "results", "data/synergyfinder", samplename))){ dir.create(file.path(.saveto, "results", "data/synergyfinder", samplename), showWarnings = FALSE, recursive = TRUE)}
        object <- synDFS[[samplename]][["synergyfinder"]][[drugname]][[synergymodel]]
        base::save(object, file = file.path(.saveto, "results", "data/synergyfinder", samplename, paste(gsub("/", " ", drugname), sep = "")))
      }
      
    }
    
    cat('\r', "Finished analyzing synergies for ", samplename, strrep(" ", 100), "\n", sep = "")
    

  }
  
  
  
  
  
  # SECTION C: Generating synergyfinder plots ##############################
  # Exporting graphs for each cell line and drug combination
  
  if(isTRUE(.plot)){
    
    for(samplename in names(synDFS)){
      
      cat("Plotting synergyfinder plots for ", samplename, ".", sep = "")
      
      
      # PLOTTING SYNERGIES
      for (drugname in names(synDFS[[samplename]][["ReshapeData"]])){
        
        
        cat('\r', "> Plotting synergyfinder plots for ", samplename, ": ", drugname, ".", strrep(" ", 50), sep = "")
        
        
        # Create subfolder if it does not exist
        if(!file.exists(file.path(.saveto, "results", "graphs/synergyfinder", samplename, gsub("/", " ", drugname)))){
          dir.create(file.path(.saveto, "results", "graphs/synergyfinder", samplename, gsub("/", " ", drugname)), showWarnings = FALSE, recursive = TRUE)
        } 
        
        
        for (block_id in 1:length(unique(synDFS[[samplename]][["drm"]][[drugname]]$drug_row))){
          
          drugpair = gsub("/", " ", unique(synDFS[[samplename]][["drm"]][[drugname]]$drug_row)[block_id])
          
          # Create subfolder if it does not exist
          if(!file.exists(file.path(.saveto, "results", "graphs/synergyfinder", samplename, gsub("/", " ", drugname), paste(gsub("/", " ", drugname), "+", drugpair, sep = " ")))){
            dir.create(file.path(.saveto, "results", "graphs/synergyfinder", samplename, gsub("/", " ", drugname), paste(gsub("/", " ", drugname), "+", drugpair, sep = " ")), showWarnings = FALSE, recursive = TRUE)
          }
          
          
          
          
          # PLOTTING DOSE RESPONSE CURVES
          # Generating dose response plots for each drug combination
          
          for (drugindex in 1:length(grep(colnames(synDFS[[samplename]][["drm"]][[drugname]]), pattern = "drug"))){
            
            # cat('\r', "> Plotting Dose-Response Curves: ", samplename, ": ", gsub("/", " ", unique(synDFS[[samplename]][["drm"]][[drugname]][which(synDFS[[samplename]][["drm"]][[drugname]]$block_id == block_id), grep(colnames(synDFS[[samplename]][["drm"]][[drugname]]), pattern = "drug")[drugindex]])), ".", sep = "")
            
            
            png(filename = file.path(.saveto, "results", "graphs/synergyfinder", samplename, gsub("/", " ", drugname), paste(gsub("/", " ", drugname), "+", drugpair, sep = " "), paste("Dose-Response Curve ", gsub("/", " ", unique(synDFS[[samplename]][["drm"]][[drugname]][which(synDFS[[samplename]][["drm"]][[drugname]]$block_id == block_id), grep(colnames(synDFS[[samplename]][["drm"]][[drugname]]), pattern = "drug")[drugindex]])), ".png", sep = "")),
                width = 1400, height = 900, units = "px", pointsize = 24)
            suppressWarnings(
              synergyfinder::PlotDoseResponseCurve(data = synDFS[[samplename]][["ReshapeData"]][[drugname]], plot_block = block_id, drug_index = drugindex, plot_new = TRUE, record_plot = FALSE)
            )
            dev.off()
            
            # cat('\r', "Finished plotting Dose-Response Curves: ", samplename, ": ", gsub("/", " ", unique(synDFS[[samplename]][["drm"]][[drugname]][which(synDFS[[samplename]][["drm"]][[drugname]]$block_id == block_id), grep(colnames(synDFS[[samplename]][["drm"]][[drugname]]), pattern = "drug")[drugindex]])), ".", sep = "")
            
          }
          
          
          
          # PLOTTING SYNERGYFINDER HEATMAPS
          # Generating dose-response matrix heatmap plots for each drug combination
          # cat('\r', "> Plotting Dose-Response Matrix Heatmap: ", samplename, ": ", paste(drugname, "+", drugpair, sep = " "), ".", sep = "")
          
          png(filename = file.path(.saveto, "results", "graphs/synergyfinder", samplename, gsub("/", " ", drugname), paste(gsub("/", " ", drugname), "+", drugpair, sep = " "), paste("Dose-Response Matrix Heatmap ", gsub("/", " ", drugname), " + ", drugpair, ".png", sep = "")), 
              width = 800, height = 800, units = "px", pointsize = 12)
          suppressWarnings(
            synergyfinder::Plot2DrugHeatmap(data = synDFS[[samplename]][["ReshapeData"]][[drugname]], plot_block = block_id, drugs = c(1, 2), plot_value = "response", dynamic = FALSE, summary_statistic = c("mean",  "median"))
          )
          dev.off()
          
          # cat('\r', "Finished plotting Dose-Response Matrix Heatmap: ", samplename, ": ", paste(gsub("/", " ", drugname), "+", drugpair, sep = " "), ".", sep = "")
          
          
          # Generating synergy score heatmap plots for each drug combination
          # cat('\r', "> Plotting ", synergymodel, " Synergy Score: ", samplename, ": ", paste(gsub("/", " ", drugname), "+", drugpair, sep = " "), ".", sep = "")
          
          png(filename = file.path(.saveto, "results", "graphs/synergyfinder", samplename, gsub("/", " ", drugname), paste(gsub("/", " ", drugname), "+", drugpair, sep = " "), paste(synergymodel, " Synergy Score Heatmap ", gsub("/", " ", drugname), " + ", drugpair, ".png", sep = "")), 
              width = 800, height = 800, units = "px", pointsize = 12)
          suppressWarnings(
            synergyfinder::Plot2DrugHeatmap(data = synDFS[[samplename]][["synergyfinder"]][[drugname]][[synergymodel]], plot_block = block_id, drugs = c(1, 2), plot_value = paste(synergymodel, "synergy", sep = "_"), dynamic = FALSE, summary_statistic = c("mean",  "median"))
          )
          dev.off()
          
          # cat('\r', "Finished plotting ", synergymodel, " Synergy Score: ", samplename, ": ", paste(drugname, "+", drugpair, sep = " "), ".", sep = "")
          
          
          
          # PLOTTING SYNERGYFINDER CONTOUR PLOTS
          # Generating dose-response matrix contour plots for each drug combination
          # cat('\r', "> Plotting Dose-Response Matrix: ", samplename, ": ", paste(drugname, "+", drugpair, sep = " "), ".", sep = "")
          
          png(filename = file.path(.saveto, "results", "graphs/synergyfinder", samplename, gsub("/", " ", drugname), paste(gsub("/", " ", drugname), "+", drugpair, sep = " "), paste("Dose-Response Matrix Contour ", gsub("/", " ", drugname), " + ", drugpair, ".png", sep = "")), 
              width = 800, height = 800, units = "px", pointsize = 12)
          suppressWarnings(
            synergyfinder::Plot2DrugContour(data = synDFS[[samplename]][["ReshapeData"]][[drugname]], plot_block = block_id, drugs = c(1, 2), plot_value = "response", dynamic = FALSE, summary_statistic = c("mean",  "median"))
          )
          dev.off()
          
          # cat('\r', "Finished plotting Dose-Response Matrix: ", samplename, ": ", paste(drugname, "+", drugpair, sep = " "), ".", sep = "")
          
          
          # Generating synergy score contour plots for each drug combination
          # cat('\r', "> Plotting ", synergymodel, " Synergy Score: ", samplename, ": ", paste(drugname, "+", drugpair, sep = " "), ".", sep = "")
          
          png(filename = file.path(.saveto, "results", "graphs/synergyfinder", samplename, gsub("/", " ", drugname), paste(gsub("/", " ", drugname), "+", drugpair, sep = " "), paste(synergymodel, " Synergy Score Contour ", gsub("/", " ", drugname), " + ", drugpair, ".png", sep = "")), 
              width = 800, height = 800, units = "px", pointsize = 12)
          suppressWarnings(
            synergyfinder::Plot2DrugContour(data = synDFS[[samplename]][["synergyfinder"]][[drugname]][[synergymodel]], plot_block = block_id, drugs = c(1, 2), plot_value = paste(synergymodel, "synergy", sep = "_"), dynamic = FALSE, summary_statistic = c("mean",  "median"))
          )
          dev.off()
          
          # cat('\r', "Finished plotting ", synergymodel, " Synergy Score: ", samplename, ": ", paste(drugname, "+", drugpair, sep = " "), ".", sep = "")
          
          
          
          # PLOTTING SYNERGYFINDER SURFACE PLOTS
          # Generating dose-response matrix surface plots for each drug combination
          # cat('\r', "> Plotting Dose-Response Matrix: ", samplename, ": ", paste(drugname, "+", drugpair, sep = " "), ".", sep = "")
          
          png(filename = file.path(.saveto, "results", "graphs/synergyfinder", samplename, gsub("/", " ", drugname), paste(gsub("/", " ", drugname), "+", drugpair, sep = " "), paste("Dose-Response Matrix Surface ", gsub("/", " ", drugname), " + ", drugpair, ".png", sep = "")), 
              width = 800, height = 800, units = "px", pointsize = 12)
          suppressWarnings(
            synergyfinder::Plot2DrugSurface(data = synDFS[[samplename]][["ReshapeData"]][[drugname]], plot_block = block_id, drugs = c(1, 2), plot_value = "response", dynamic = FALSE, summary_statistic = c("mean",  "median"))
          )
          dev.off()
          
          # cat('\r', "Finished plotting Dose-Response Matrix: ", samplename, ": ", paste(drugname, "+", drugpair, sep = " "), ".", sep = "")
          
          
          # Generating synergy score surface plots for each drug combination
          # cat('\r', "> Plotting ", synergymodel, " Synergy Score: ", samplename, ": ", paste(drugname, "+", drugpair, sep = " "), ".", sep = "")
          
          png(filename = file.path(.saveto, "results", "graphs/synergyfinder", samplename, gsub("/", " ", drugname), paste(gsub("/", " ", drugname), "+", drugpair, sep = " "), paste(synergymodel, " Synergy Score Surface ", gsub("/", " ", drugname), " + ", drugpair, ".png", sep = "")), 
              width = 800, height = 800, units = "px", pointsize = 12)
          suppressWarnings(
            synergyfinder::Plot2DrugSurface(data = synDFS[[samplename]][["synergyfinder"]][[drugname]][[synergymodel]], plot_block = block_id, drugs = c(1, 2), plot_value = paste(synergymodel, "synergy", sep = "_"), dynamic = FALSE, summary_statistic = c("mean",  "median"))
          )
          dev.off()
          
          # cat('\r', "Finished plotting ", synergymodel, " Synergy Score: ", samplename, ": ", paste(drugname, "+", drugpair, sep = " "), ".", sep = "")
          
          
          
          
          # PLOTTING SYNERGYFINDER DOSE RESPONSE COMPOSITE
          # Generating dose-response matrix surface plots for each drug combination
          # cat('\r', "> Plotting Dose-Response: ", samplename, ": ", paste(drugname, "+", drugpair, sep = " "), ".", sep = "")
          
          png(filename = file.path(.saveto, "results", "graphs/synergyfinder", samplename, gsub("/", " ", drugname), paste(gsub("/", " ", drugname), "+", drugpair, sep = " "), paste("Dose-Response ", gsub("/", " ", drugname), " + ", drugpair, ".png", sep = "")), 
              width = 1600, height = 869, units = "px", pointsize = 12)
          suppressWarnings(
            synergyfinder::PlotDoseResponse(data = synDFS[[samplename]][["ReshapeData"]][[drugname]], block_ids = block_id, drugs = c(1, 2), save_file = FALSE)
          )
          dev.off()
          
          # cat('\r', "Finished plotting Dose-Response: ", samplename, ": ", paste(drugname, "+", drugpair, sep = " "), ".", sep = "")
          
          
          
          
        }
        
      }
      
      cat('\r', "Finished plotting synergyfinder plots for ", samplename, ".", strrep(" ", 50), "\n", sep = "")
      
    }
    
  }
    
  
  .synergyscore = list()
  
  # Build output data object from synergyfinder data
  for (samplename in names(synDFS)){
    for (drugname in names(synDFS[[samplename]][["ReshapeData"]])){
      .synergyscore[[samplename]][[drugname]] <- synDFS[[samplename]][["synergyfinder"]][[drugname]][[synergymodel]][["drug_pairs"]]
      .synergyscore[[samplename]][[drugname]]$`Cell Type` <- samplename
      .synergyscore[[samplename]][[drugname]]$`package` <- "synergyfinder"
      .synergyscore[[samplename]][[drugname]]$`version` <- as.character(packageVersion("synergyfinder"))
      .synergyscore[[samplename]][[drugname]]$model <- synergymodel
    }
    .synergyscore <- do.call(rbind, setNames(.synergyscore[[samplename]], NULL))
  }
  

    
  # Remove duplicate drug pairs and sort by maximum score
  .synergyscore <- .synergyscore[!duplicated(t(apply(.synergyscore[c("drug1", "drug2")], 1, sort))),]
  
  # Sort the data set by synergy score and rank them accordingly
  # Select the column to sort by  dynamically based on the synergy model
  .synergyscore <- .synergyscore[do.call(order, -.synergyscore[paste(synergymodel, "synergy", sep = "_")]), ]
  .synergyscore$rank <- 1:nrow(.synergyscore)
  
  # Select only relevant columns and rename those columns
  .synergyscore <- .synergyscore[,c("Cell Type", "package", "version", "model", "rank", "drug1", "drug2", "input_type", grep("p_value", names(.synergyscore), value = TRUE), grep("synergy$", names(.synergyscore), value = TRUE))]
  names(.synergyscore) <- c("Cell Type", "Package", "Version", "Model", "Rank", "Drug A", "Drug B", "Response Type", "p-value", "Synergy Score")
  
  
  # Assign data to an object
  synDFS[[samplename]][["synergyscores"]] <- .synergyscore
  

  # Return object of class synergyfinderData
  class(synDFS) <- "synergyfinderData"
  return(synDFS)
    
  
}
