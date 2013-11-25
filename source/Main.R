
# Load the TclTk package
require(tcltk)

# Check for required libraries, install if not installed.
# (Taken care of in install file).

# Load libraries
library(tkrplot)
#library(drc)
library(gplots)

# Parallel processing
# (need to include in installation script)
show ("Loading parallel processing library")
library(foreach)
#library(doMC)
#registerDoMC()



# Source files used in the program:
source("Analysis3D.R")
source("Analysis2D.R")
source("DataPrep.R")
source("FinalAnalysis.R")
source("HillFitFunctions.R")
source("Colorbar.R")
source("BinPartition.R")
source("MultiFileEntry.R")
source("BatchPrep.R")
source("Preferences.R")
source("Colorspec.R")
source("LogPlotting.R")





# Set directories appropriately

#StorePrefs()
#SetPrefs()
#setwd(Preferences$home.directory)
programDir		<- getwd()




################################################################################
## Create variables#############################################################


## VARS ##########################################

# Apply across titrations titrations equally
Vars 	<- list()
Vars[["num.doses"]] 		<- 6	# Number of doses per titration
Vars[["num.titr"]]			<- 2	# Number of titrations to analyze in parallel

Vars[["programDir"]]		<- programDir
Vars[["PlotOptions"]]		<- list()
Vars[["cores"]]		<- getDoParWorkers()
Vars[["islog"]]		<- list()
Vars$islog$x			<- T
Vars$islog$y			<- T
Vars$islog$z			<- T
Vars[["stats"]]		<- list()
Vars$stats$minx			<- 0
Vars$stats$miny			<- 0
Vars$stats$minz			<- 0
Vars$stats$maxx			<- 4
Vars$stats$maxy			<- 4
Vars$stats$maxz			<- 10000
Vars$stats$samplesize	<- 10000
Vars[["key"]]<- list()
Vars$key[["short"]] 	<- c("ec","amp","hill")
Vars$key[["long"]] 	<- c("EC50","Amplitude","Hillslope")

Filefolder	<- "~"

# Variables - per titration ######
# (each will be a list with as many elements as there are titrations)

# Set in MultiFileEntry
Vars[["Datafiles"]] 	<- list()	# Where headers for datafiles will be kept temporarily
Vars[["Datas"]] 		<- list()	# After files are entered, contains only selected x, y and z columns
Vars[["units"]]			<- list()
Vars[["Filenames"]] 	<- list()
Vars[["Concs"]] 		<- list()
Vars[["Control.File"]]	<- list()
Vars$units[[1]]		<- "pM"
Vars[["prefix"]]	<- list()
Vars$prefix[[1]]	<- "Titration 1"

# For DataPrep
Vars[["Channels.orig"]]	<- list()	# Names of all channels on exported files
Vars[["Channels"]]		<- list()	# Names of all channels entered by user (in DataPrep)
Vars[["Chan.track"]]	<- list()	# List of channels, used to track log, lin, x,y,z, etc.
Vars[["xyzLabels"]] 	<- list()	# Stores the names of only the x, y, and z channels

jet.colors 	<- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan",
                     "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
Jet <- jet.colors(32)





# Tracking progress through the program, used for deciding whether to use default
# or user specified values.
Track	<- list()
Track[["FirstDataprep"]] 	<- T
Track[["FirstEntry"]] 		<- T
Track[["FirstInitialPlot"]]	<- T





### Options for plotting etc.  #####################################
Options	<- list()

Options[["parallelize"]] <- T
Options[["thresh"]]	<- 10

Options[["fit"]]	<- "means"
# determines whether hill fits will be performed on all data ("all"), means of
# data ("means") or a percentage of positive events, basd on a threshold ("pos_pct")

Options[["plots"]]	<- list()
Options$plots$min$x	<- 0
Options$plots$min$y	<- 0
Options$plots$minz	<- 0
Options$plots$max$x	<- 4
Options$plots$max$y	<- 4
Options$plots$maxz	<- 10^4

Options$plots$initialZmin	<- 1
Options$plots$initialZmax	<- 10^5

Options$plots[["islog"]]	<- list()
Options$plots$islog$ec		<- T
Options$plots$islog$amp		<- T
Options$plots$islog$hill	<- F

Options$plots$showz		<- F
Options$plots$zcutoff	<- 10^4

Options$plots$pcutoff		<- 1
Options$plots$ecREcutoff.idx	<- 1e5
Options$plots$ampREcutoff.idx	<- 1e5
Options$plots$p.AndOr.RE	<- "&" #alternate value is "|"

#2D plot options
Options$plots$hist2d	<- "off" #alternates = "bw" or "col"

Options[["bins"]]	<- list()
Options$bins$min$x	<- 0
Options$bins$min$y	<- 0
Options$bins$max$x	<- 4
Options$bins$max$y	<- 4
Options$bins$n$x		<- 20
Options$bins$n$y		<- 20

Options[["analysis"]]	<- list()
Options$analysis$zchoice	<- "mfi" #vs. "pct"
Options$analysis$pctgate	<- 2
Options$analysis$hillchoice	<- "fixed" # vs. "free": value for hill coefficient will be fixed at 1.
Options$analysis$method		<- "hillfit" # or "subtract": will be hillfit unless too few doses are entered (3 or fewer)






### HOLDS BINNING RESULTS ######################################################
Bins	<- list()
# Parameters applying to all titrations
Bins$nx	<- 20	# Number of bins (n)
Bins$ny	<- 20
Bins$seqx	<- c()	# Sequence of bins, length=n+1
Bins$seqy	<- c()
Bins$centersx	<- c() # calculated centers of bins, length=n
Bins$centersy	<- c()

# Used to determine whether to re-bin or not
Bins$params	<- c("0") 	# Concatenated current binning parameters for each titration
						# Used to determine if binning has been done or if rebinning is necessary.

# One list per titration                                   
Bins$sizes	<- list()	# Number of events per bin, Titration level list
Bins$sizes[[1]]	<- list()	# Dose
Bins$sizes[[1]][[1]] <- matrix(0)	# Matrix for bins

Bins$means	<- list()	# Z value means Titration level
Bins$means[[1]]	<- list()	#Dose level
Bins$means[[1]][[1]] <- matrix(0)	#

Bins$sigmas	<- list()	# Z value SD; per Titration
Bins$sigmas[[1]]	<- list()	#Dose level
Bins$sigmas[[1]][[1]] <- matrix(0)	#

Bins$pos.cut.mat	<- list()	# Threshold for determining "positive" cells, set by control file
Bins$pos.cut.mat[[1]]	<- matrix(0)

Bins$pospct	<- list()	# Percent positive, determined from control file -> pos.cut.mat
Bins$pospct[[1]] <-list()
Bins$pospct[[1]][[1]] <- matrix(0)

Bins$Full	<- list()	# Contains z values for final binning (large list)
Bins$Full[[1]]	<- list() # Each titration
Bins$Full[[1]][[1]]	<- list() # One per dose, named with X and Y coordinates

# 2D binning variables
Bins[["x"]] <- list()
Bins$x[["Full"]]	<- list()
Bins$x[["means"]]	<- list()
Bins$x[["sizes"]]	<- list()
Bins$x[["pospct"]]	<- list()
Bins$x[["sigmas"]]	<- list()

Bins[["y"]] <- Bins$x

# Bins$x and $y will have a sub-list for each titration, with Full, means, sizes, and pospct
# each of the above will have a vector 

# Run program
startBox(Vars,Track,Options,Bins)





