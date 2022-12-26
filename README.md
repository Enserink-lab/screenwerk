# screenwerk: a modular pipeline for drug sensitivity screens

An R package with a modular library for the design and analysis of high-throughput drug combination screens


code name: screenwerk   dev: alpha 0.2.0-1   author: roberthanes


<img src="https://github.com/Enserink-lab/screenwerk/blob/main/doc/screenwerk.png?raw=true" width="100%" align="left"></img><br />
<br>


The R-package consists of individual functions that can be used to (a) set-up a screen and generate a dispensing file, (b) to read raw measurements and consolidate different data sets, (c) perform a quality control, (d) analyze, (e) report and visualize the results of a drug combination screen.

#### INSTALLATION

```r
install.packages('devtools')
library(devtools)

# Install R-package from github
devtools::install_github('Enserink-lab/screenwerk', build = TRUE, build_opts = c("--no-resave-data", "--no-build-vignettes"))
# Load 'screenwerk' package
library(screenwerk)

```

A set of files can be found in inst/extdata/ and used as a reference for building data sets.



#### Workflow

The R-package consists of individual functions that can be used to (a) set-up a screen and generate a dispensing file, (b) to read raw measurements and consolidate different data sets, (c) perform a quality control, (d) analyze, (e) report and visualize the results of a drug screen.

#### ***(a) pre-experimental: setting-up a screen and generating a dispensing file***

A set of files can be found in inst/extdata/ and used as examples.

In order to generate a dispensing file for a drug combination screen, we need a list of combinations, among a few other files. All files are provided as .csv files. The list of combinations can be easily generated using one of the functions in this package.

First we need a list of drugs at different doses to be used in the screen, which will be combined with each other at selected doses.
This list can be provided in a wide format, as it would be done using Microsoft Excel, or any other spreadsheet tool.


This file can look as shown below:
<img src="https://github.com/Enserink-lab/screenwerk/blob/main/doc/figures/doses-wide.png?raw=true" align="left"></img><br>
<br>
<br>
<br>
<br>
<br>
<br>
<br>


Now we need to import this file into R using the code below:

```r
listofDoses <- read.csv(file=file.path("inst/extdata/library/listofDoses.csv"), check.names=FALSE, header=TRUE, stringsAsFactors=FALSE, na.strings="", sep=",", dec=".", skip=0)
```
and run the following function, which converts the list from a human readable wide-format to a long-format:
```r
listofDoses <- generateListofDoses(listofDoses, .dropCol = TRUE)
```
Note: *.dropCol = TRUE* drops columns that are not needed.

This will generate a data set that looks like the one below:
<br>
<img src="https://github.com/Enserink-lab/screenwerk/blob/main/doc/figures/doses-long.png?raw=true" align="left"></img><br>
<br>
<br>
<br>
<br>



Now that we have a list of doses, we can use it to generate a list of combinations:
```r
listofCombinations <- combineDrugs(listofDoses, .combineDoses=c(2:5), .noReplicates = 3, .drugRepAttrib = "single")
```
In the function above, we can provide the list of doses we just generated through the first argument. With *.combineDoses=c(2:5)* we specify, which of the doses we want to combine, in that case the doses 2 to 5 are combined, leaving out the lowest and highest dose out. With *.noReplicates = 3* we specify how many replicates we want, and with *.drugRepAttrib = "single"* we can specify that only single drug treatments should be replicated.

This will generate the following data set, with a drug, dose and unit column for each drug pair:
<img src="https://github.com/Enserink-lab/screenwerk/blob/main/doc/figures/combinations.png?raw=true" align="left"></img><br />
<br>
<br>
<br>
<br>
<br>

Now we are ready to use the list of combinations we just generated, to build a dispensing data set. However, in order to do that, we need to provide a few additional files. First, let's generate a list of wells we want to exclude from being dispensed into. This can be achieved using the function below:
```r
listofExWells <- excludeWells(1536, outer.wells=TRUE)
```

This will generate the following list of wells:
```r
> listofExWells
[1] "A1"   "A2"   "A3"   "A4"   "A5"   "A6"   "A7"   "A8"   "A9"   "A10"  "A11"  "A12"  "A13"  "A14"  "A15"  "A16"  "A17"  "A18"  "A19"  "A20"  "A21"  "A22"  "A23"  "A24"  "A25"
[26] "A26"  "A27"  "A28"  "A29"  "A30"  "A31"  "A32"  "A33"  "A34"  "A35"  "A36"  "A37"  "A38"  "A39"  "A40"  "A41"  "A42"  "A43"  "A44"  "A45"  "A46"  "A47"  "A48"  "B1"   "B48"
[51] "C1"   "C48"  "D1"   "D48"  "E1"   "E48"  "F1"   "F48"  "G1"   "G48"  "H1"   "H48"  "I1"   "I48"  "J1"   "J48"  "K1"   "K48"  "L1"   "L48"  "M1"   "M48"  "N1"   "N48"  "O1"  
[76] "O48"  "P1"   "P48"  "Q1"   "Q48"  "R1"   "R48"  "S1"   "S48"  "T1"   "T48"  "U1"   "U48"  "V1"   "V48"  "W1"   "W48"  "X1"   "X48"  "Y1"   "Y48"  "Z1"   "Z48"  "AA1"  "AA48"
[101] "AB1"  "AB48" "AC1"  "AC48" "AD1"  "AD48" "AE1"  "AE48" "AF1"  "AF2"  "AF3"  "AF4"  "AF5"  "AF6"  "AF7"  "AF8"  "AF9"  "AF10" "AF11" "AF12" "AF13" "AF14" "AF15" "AF16" "AF17"
[126] "AF18" "AF19" "AF20" "AF21" "AF22" "AF23" "AF24" "AF25" "AF26" "AF27" "AF28" "AF29" "AF30" "AF31" "AF32" "AF33" "AF34" "AF35" "AF36" "AF37" "AF38" "AF39" "AF40" "AF41" "AF42"
[151] "AF43" "AF44" "AF45" "AF46" "AF47" "AF48"
```
In the first argument, we specify the plate format we plan to use. With *outer.wells=TRUE* we simply state that we want to exclude all outer wells. The function will then automatically generate the list. If additional or  a custom set of wells need to be excluded, have a look at the R documentation for the package 'screenwerk'.
<br>
<br>


Now that this is done, we need to provide the following files, in addition to the ones already generated:

 - a list of drugs
 - a list of volumes
 - a list of controls
 - a list of stock concentrations
 - one or more source plates
 <br>

**The list of drugs can look as shown below:**<br>
<img src="https://github.com/Enserink-lab/screenwerk/blob/main/doc/figures/drugs.png?raw=true" align="left"></img><br />
<br>
<br>
<br>
<br>



**The list of volumes, can be provided as the list of doses in a wide-format:**<br>
<img src="https://github.com/Enserink-lab/screenwerk/blob/main/doc/figures/volumes.png?raw=true" align="left"></img><br />
<br>
<br>
<br>
<br>
<br>



**The list of controls:**<br>
<img src="https://github.com/Enserink-lab/screenwerk/blob/main/doc/figures/controls.png?raw=true" align="left"></img><br />
<br>
<br>
<br>


**The list of stock concentrations:**<br>
<img src="https://github.com/Enserink-lab/screenwerk/blob/main/doc/figures/concentrations.png?raw=true" align="left"></img><br />
<br>
<br>
<br>
<br>


Let's import all those files:
```r
listofDrugs <- read.csv(file=file.path("inst/extdata/library/listofdrugs.csv"), check.names=FALSE, header=TRUE, stringsAsFactors=FALSE, na.strings="", sep=",", dec=".", skip=0)
listofVolumes <- read.csv(file=file.path("inst/extdata/library/listofVolumes.csv"), check.names=FALSE, header=TRUE, stringsAsFactors=FALSE, na.strings="", sep=",", dec=".", skip=0)
listofCtrls <- read.csv(file=file.path("inst/extdata/library/listofctrls.csv"), check.names=FALSE, header=TRUE, stringsAsFactors=FALSE, na.strings="", sep=",", dec=".", skip=0)
listofStockConcentrations <- read.csv2(file=file.path("inst/extdata/library/listofStockConcentrations.csv"), check.names=FALSE, header=TRUE, stringsAsFactors=FALSE, na.strings="", sep=",", dec=".", skip=0)
```
<br>

In order to import a source plate from a .PlateMap file, we can use the following function:
```r
sourcePlate <- importPlateMap("inst/extdata/library/", .sourcePlateConv = TRUE)
```
The first argument sets the location of one or more .PlateMap files, which to import.\
Note: The function will import all .PlateMap files it finds in that folder.
<br>
<br>
This will generate a data set as shown below:<br>
<img src="https://github.com/Enserink-lab/screenwerk/blob/main/doc/figures/sourceplate.png?raw=true" align="left"></img><br />
<br>
<br>
<br>
<br>
<br>
<br>


Now we are ready to generate the dispensing data set:
```r
dispensingData <- generateDispensingData(listofCombinations, listofDrugs, listofDoses, listofVolumes, listofCtrls, listofStockConcentrations, sourcePlate,
                       listofExWells, .ctrlReplicates = 8, .addUntreated = list(name = "Untreated", replicates = 8),
                       .finalWellVolume = 5, .plateFormat = 1536, .destinationPlateID = "0521",
                       .randomizeDispensing = TRUE, .probeDispensing = FALSE)
```
A short explanation of the arguments:\
**ctrlReplicates**, sets the number of replicates for the controls provided in the list of controls *listofCtrls*\
**addUntreated**, includes a set of untreated controls, by specifying the name and number of replicates.\
**finalWellVolume**, the final volume in all the wells. This is important for mapping the right dose on the source plate to the dose on the destination plate.\
**plateFormat**, the plate format used for the destination plates.\
**destinationPlateID**, a unique destination plate ID, which will be used to generate the destination plate barcode.\
**randomizeDispensing**, specifying that the drug treatments should be randomized across all plates.\
**probeDispensing**, specifying, if the dispensing should be probed before generating the data set, which will provide a short summary with the number of plates, drug treatments and much more.

For a more detailed explanation, have a look at the R documentation for the package 'screenwerk'.

Now we should have an R object of class S3, with the dispensing data, as well as the original input data and a number of lists included.
<br>
<br>
<br>


In order to *save* the dispensing data, *print* a summary, *plot* or *export* the dispensing files, we can use a number of different functions:\
<br>
**print**, to print and extract the dispensing data set\
**summary**, to show a dispensing summary, with the number of plates, drug treatments and much more.\
**plot**, to plot the dispensing layout for each plate.\
**save**, to save and export the dispensing data to individual files for each set.
<br>

*Extracting the dispensing data:*
```r
print(dispensingData)
```
<br>

*Showing a summary for the dispensing:*
```r
summary(dispensingData)
```
<br>

*Plotting the dispensing layout for each plate:*
```r
plot(dispensingData, .saveto = "results\myDispeningPlots")
```
<br>

*Generating a set of dispensing files:*
```r
save(dispensingData, .saveto = "results\myDispeningFiles", .sets = 12, .labels = "alphabetic", .split = FALSE, .format = "Echo")
```

**.sets**, specifies how many sets should be generated.\
**.labels**, provides the labels for each set. The labels can be predefined, or provided as custom labels. See package documentation for more information. In the example above, we use a list of alphabetic letters, A-L, for each of the 12 sets.\
**.split**, if the dispensing file should be split into individual files for each set, otherwise, a single file will be generated.\
**.format**, specifies the final format of the dispensing file. The format is dependent on the machine and robot used to carry out the dispensing. For more information, see package documentation.
<br>
<br>

#### ***(b) post-experimental: reading raw data and consolidating data sets***

Once an experiment has been carried out and the assay read, we can use the function below to read all raw data files, in order to be able to build a master data set by consolidating the dispensing data with the raw measurements.


Reading raw data:
```r
rfs <- readRAWData(.readfrom = "inst/extdata/raw", .fileformat = c(".csv", ".txt"), .format = "EnVision")
```
**.readfrom**, set the folder location to the raw data files.\
**.fileformat**, set the file format of those files. The function supports the most common text-based file types (.txt, .csv) in various formats.\
**.format**, select the format of the read-out files. This argument is machine specific, since each manufacturer and plate reader will use a different, sometimes proprietary, format. Alternatively, the data can be provided in a raw data format, independent of the machine with the raw measurements based on the well plate layout. The function is capable of detecting the export-format and identifying the used text delimiter, such as comma-, semicolon- or tab-separated, as well as the plate format and find the raw data of interest. For more information, see package documentation.

<br>
<br>

With the raw measurements read, we can now build the final reference data set, used for downstream analysis of the drug sensitivity screen. This can be achieved by consolidating the raw measurements with the dispensing data.

Before we do that, we need to import a barcode reference list with the names of the samples used in the drug screen and by associating them to the corresponding plate id and set.

This file can look as shown below:<br>
<img src="https://github.com/Enserink-lab/screenwerk/blob/main/doc/figures/barcode.png?raw=true" align="left"></img><br>
<br>
<br>
<br>
<br>

and imported with:
```r
.barcodeReference <- read.csv(file=file.path(libDirectory, "platebarcode.csv"), check.names=FALSE, header=TRUE, stringsAsFactors=FALSE, colClasses=c("PlateID"="character"), comment.char = "#", blank.lines.skip	= TRUE, na.strings="", sep=",", dec=".", skip=0)
```
<br>
<br>

Now we can consolidate the raw measurements with the dispensing data, using the function below:
```r
clData <- consolidateData(dispensingData = dispensingData, rawMeasurements = rfs, .barcodeReference)
```
**.dispensingData**, an object of class 'dispensingData'.\
**.rawMeasurements**, an object of class 'rawMeasurements'.\
**.barcodeReference**, a list of samples for each set of plates.
<br>
<br>
Once consolidated, the data can be processed for downstream analysis with screenwerk, or exported for the use with other tools. This can be achieved by using again *save* function, as shown below:
<br>

*Export and save data set to file:*
```r
save(clData, .saveto = "export\", .fileformat = ".csv", .sep = ";")
```
<br>
