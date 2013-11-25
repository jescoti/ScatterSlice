BinPartition3d <- function(Bins,Vars, Options, titr){
	show("BinPartition3d")
	
	
		
	#show(sprintf("Binning titration %i", titr))
	#show(sprintf("length of Bins$sizes=%i", length(Bins$sizes)))
	
	Partition <- function(filenum, titr, Bins, datamatrix, Options, controlfile=F) {
		#show(sprintf("Binning and analyzing file %d of titration %d", filenum, titr))
		ScsBins	<<- Bins
		
		
		
		totalbins 		<- Bins$nx*Bins$ny	# Counter for matrix sizing
		Binz		<- list()
		Binz$Full	<- list()
		mat			<- matrix(rep(0,times=totalbins), ncol=Bins$ny)
		sizemat 	<- mat  # number of cells in bin
		meanmat		<- mat	# mean of Z variable
		sigmamat	<- mat	# sd of a bin
		logmeanmat	<- mat	# geometric mean
		medmat		<- mat	# median
		
		pos.pctmat	<- matrix(rep(NA,times=totalbins), ncol=Bins$ny)	# percentage positive cells
		pos.cut.mat	<- mat	# 
		
		if (controlfile){
			pos.pctmat <- matrix(rep(Options$analysis$pctgate,times=totalbins), ncol=Bins$ny)
		}
		
		# d is a temporary matrix to hold the values
		# of Data that will be processed for binning
		
		
		# Splitting seems like a fast way to bin, but I don't yet have it working.
		#d<-as.data.frame(datamatrix)
		#names(d)	<- c("x", "y", "z")
		#SPLITTEST <<- split(d, findInterval(d$x,Bins$seqx))
		#lowbin	<- which(names(SPLITTEST)=="1")
		#SPLITTEST	<<- SPLITTEST[lowbin:(lowbin+Bins$nx-1)]
		#show(sprintf("number of bins = %d, number of split bins = %d", Bins$nx, length(SPLITTEST)))
		#show(sprintf("median of bottom bin = %.3g", log10(median(SPLITTEST[[1]]$z))))
		#show(sprintf("lowest bin edges = %.3g, %.3g", Bins$seqx[1], Bins$seqx[2]))
		
		
		
		d	<- datamatrix
		
		##### Binning function
		#show("Dividing Bins")
		for (y.idx in 1:Bins$ny){
			temp.low.y 	<- matrix(d[d[,2] > Bins$seqy[y.idx],], ncol=3)
			temp.high.y	<- matrix(temp.low.y[temp.low.y[,2]<= Bins$seqy[y.idx+1],], ncol=3)
			for (x.idx in 1:Bins$nx){
				temp.low.x 	<- matrix(temp.high.y[temp.high.y[,1] > Bins$seqx[x.idx],], ncol=3)
				final.bin	<- matrix(temp.low.x[temp.low.x[,1]<= Bins$seqx[x.idx+1],], ncol=3)
				FINALBIN	<<- final.bin
				
				
				#meanmat[x.idx, y.idx]	<- mean(final.bin[,3])
				
				finalpos	<- final.bin[,3][final.bin[,3]>0]	# Only take positive numbers
				themean		<- 10^mean(log10(finalpos)) # Geometric mean of positive numbers 
				meanmat[x.idx,y.idx]	<- themean
				
				themedian	<- median(final.bin[,3])	# Planning to switch to median
				medmat[x.idx,y.idx]		<- themedian
				
				### Planning to switch out mean for median
				#if(!is.na(themean) && !is.na(themedian)){
				#if(themean != 0 && themedian != 0){
				#	
				#	show(sprintf("bin %d, %d  |  mean=%.6g, median = %.6g, ratio= %.6g", x.idx, y.idx, themean, themedian, themean/themedian))
				#}
				#}
				
				binsize		<- length(final.bin)
				sizemat[x.idx, y.idx] 	<- binsize
				
				if(binsize>2){
					sigmamat[x.idx,y.idx]	<- sd(log10(finalpos)) 
				}
				else {sigmamat[x.idx,y.idx]	<- 0}
				
				
				final.z <- final.bin[,3]
				bin.idx			<- sprintf("%i,%i", x.idx, y.idx)
				##show(bin.idx)
				##show(head(final.bin))
				#show(sprintf("cells in bin = %d, threshold = %d", length(final.z), Options$thresh))
				Binz$Full[[bin.idx]]	<- final.bin[,3]
				if (length(final.z)>Options$thresh){
					if(controlfile){
						# Sets cutoff for positive cells based on control
						pos.cut.mat[x.idx, y.idx] <- quantile(final.z, (100-Options$analysis$pctgate)/100, na.rm=T)
					}else{
						# Evaluates percentage of positive cells
						if (Bins$pos.cut.mat[[titr]][x.idx, y.idx] > 0){
							pos.pctmat[x.idx, y.idx]	<- length(final.z[final.z > Bins$pos.cut.mat[[titr]][x.idx, y.idx]])/length(final.z)*100
							#show("control had enough cells")
						}else{
							pos.pctmat[x.idx, y.idx] <- NA
							#show("control didn't have enough cells")
						}
					}
				}
			}
		}
		#show("Division completed")
		# Store statistics of bins in appropriate entry in Bins$sizes or Bins$means
		Binz$means 		<- meanmat
		Binz$sizes 		<- sizemat
		Binz$pospct		<- pos.pctmat
		Binz$pos.cut.mat	<- pos.cut.mat
		Binz$filenum	<- filenum
		Binz$sigmas		<- sigmamat
		Binz$logmeans	<- logmeanmat
		Binz$medians	<- medmat
		#show("Binning completed")
		#Binz$Full<-NA
		
		return(Binz)
		
		#show("Sizes:")
		#show(sizemat)
		
	}
	
	
	current.params	<- paste(	Options$bins$min$x,Options$bins$min$y,Options$bins$max$x,
									Options$bins$max$y,Options$bins$n$x,Options$bins$n$y, 
									Vars$islog$x, Vars$islog$y, Vars$islog$z,
									sep="_")
	
	if(Bins$params[titr] != current.params){
		Bins$params[titr]	<- current.params
		#show("New binning parameters -> binning")
		#show(sprintf("%d X, %d Y bins, %d total", Bins$nx,Bins$ny,Bins$nx*Bins$ny))	
		starttime <- proc.time()[1] #for timing
		
		
		## Bin "Control" file
		ctl			<-Vars$Control.File[[titr]]
		ctl.dmtx	<- as.matrix(Vars$Datas[[titr]][[ctl]][,1:3])
		tempBins	<- Partition(filenum=ctl, titr, Bins, datamatrix=ctl.dmtx, Options, controlfile=T)
		Bins$means[[titr]][[ctl]] 	<- tempBins[["means"]]
		Bins$sizes[[titr]][[ctl]] 	<- tempBins[["sizes"]]
		Bins$pospct[[titr]][[ctl]]	<- tempBins[["pospct"]]
		Bins$sigmas[[titr]][[ctl]]	<- tempBins[["sigmas"]]
		Bins$Full[[titr]][[ctl]]	<- tempBins[["Full"]]
		Bins$pos.cut.mat[[titr]]	<- tempBins[["pos.cut.mat"]]
		
		
		
		means 	<- Bins$means[[titr]][[ctl]]
		sizes	<- Bins$sizes[[titr]][[ctl]]
		
		sizecheck.tf	<- sizes >= Options$thresh
		meancheck.tf	<- !is.nan(means)
		binsmin <- min(means[sizecheck.tf & meancheck.tf])
		binsmax	<- max(means[sizecheck.tf & meancheck.tf])
		
		ScsBins	<<- Bins
		#show("Control file binned, sizes: ")
		#show(Bins$y$sizes[[titr]][[ctl]])
		
		## Bin all other files
		filevec	<- c(1:Vars$num.doses)
		treated.filenums	<- filevec[-Vars$Control.File[[titr]]]
		if (Options$parallelize==T) {
			# Use foreach package to partition bins in parallel
			
			parallelBins <- foreach(filenum=treated.filenums) %do% {  ### SHOULD BE ABLE TO PARALLELIZE THIS, BUT IT'S REALLY SLOW...
				#show(sprintf("Binning file %d",filenum))
				Partition(filenum, titr, Bins, datamatrix=as.matrix(Vars$Datas[[titr]][[filenum]][,1:3]), Options, controlfile=F)
				
			}
			
			PARBINS<<-parallelBins
			
			non.control.files	<- filevec[-Vars$Control.File[[titr]]]
			
			for (file.idx in 1:length(parallelBins)) {
				
				filenum	<- parallelBins[[file.idx]][["filenum"]]
				
				#show(sprintf("moving data for file %i",filenum))
				
				Bins$means[[titr]][[filenum]] 	<- parallelBins[[file.idx]][["means"]]
				Bins$sizes[[titr]][[filenum]] 	<- parallelBins[[file.idx]][["sizes"]]
				Bins$pospct[[titr]][[filenum]]	<- parallelBins[[file.idx]][["pospct"]]
				Bins$Full[[titr]][[filenum]]	<- parallelBins[[file.idx]][["Full"]]
				Bins$sigmas[[titr]][[filenum]]	<- parallelBins[[file.idx]][["sigmas"]]
				
				#show("moved data")
				means <- Bins$means[[titr]][[filenum]]
				sizes <- Bins$sizes[[titr]][[filenum]]
				sizecheck.tf	<- sizes >= Options$thresh
				meancheck.tf	<- !is.nan(means)
				##show(means)
				##show(sizecheck.tf)
				##show(meancheck.tf)
				#AND <- sizecheck.tf & meancheck.tf
				##show(binsmin)
				
				binsmin <- min(means[sizecheck.tf & meancheck.tf], binsmin)
				binsmax	<- max(means[sizecheck.tf & meancheck.tf], binsmax)
				#show("finished check")
			}
			
			BINSPAR<<-Bins
			#show("finished binning")

			
		}else{
			
			for (filenum in filevec[-Vars$Control.File[[titr]]]) {
				Bins <- Partition(filenum, titr, Bins, datamatrix=as.matrix(Vars$Datas[[filenum]][,1:3]), Options, controlfile=F)
				
				means <- Bins$means[[titr]][[filenum]]
				sizes <- Bins$sizes[[titr]][[filenum]]
				sizecheck.tf	<- sizes >= Options$thresh
				meancheck.tf	<- !is.nan(means)
				binsmin <- min(means[sizecheck.tf & meancheck.tf], binsmin)
				binsmax	<- max(means[sizecheck.tf & meancheck.tf], binsmax)
			}
		}
		
		#### INCLUDED AS AN OPTION TO COMPARE PARALLEL AND SEQUENTIAL BINNING TIMES
		#compareParallel<-F
		#if(compareParallel==T){
		#	sequentialtime	<- system.time(
		#	foreach(filenum=filevec[-Vars$Control.File[[titr]]]) %do% {
		#		Partition(filenum, titr, Bins, datamatrix=as.matrix(Vars$Datas[[titr]][[filenum]][,1:3]), Options, controlfile=F)
		#	})[3]
		#	show(sprintf("Parallel binning took %.3f seconds", paralleltime))
		#	show(sprintf("Sequential binning takes %.3f seconds", sequentialtime))
		#}
		
		BINS<<-Bins
		
		
		Bins$binsmin <- binsmin
		Bins$binsmax <- binsmax
		
		ScsDatas <<- Vars$Datas
		
		#show("Finished Binning")
		##show(sprintf("Filled %i out of %i total bins.", length(Bins$sizes[Bins$sizes>0]),as.integer(totalbins)))
	}
	#show(Bins$y$sizes[[1]][[1]])
	return(Bins)
}

#########################################################################################################
#########################################################################################################

BinPartition2d <- function(Bins, Vars, Options, titr){
	#show("BinPartition2d")
	show(sprintf("Binning titration %i",titr))
	starttime <- proc.time()[1]
	
	# Set up matrices to receive values
	# sizeMat holds number of values in each bin
	# meanMat holds Z variable mean in each bin
	
	
	Partition <- function(filenum, titr, Bins, Vars, Options, xy) {
		#show("names(Bins)=")
		#show(names(Bins))
		nbins <- Bins[[paste("n", xy, sep="")]]
		n<-nbins
		vect		<- rep(0,times=nbins)
		Binz	<- list()
		Binz$Full	<- list()
		sizevect 	<- vect
		meanvect	<- vect
		idxvect		<- vect
		pos.pctvect	<- vect
		pos.cutoff	<- vect
		sigmavec	<- vect
		
		seqxy		<- paste("seq",xy, sep="")
		
		if (filenum == Vars$Control.File[[titr]]){
			pos.pctvect	<- rep(Options$analysis$pctgate,times=n)
		}
		
		
		# d is a temporary matrix to hold the values
		# of Data that will be processed for binning
		
		d <- as.data.frame(Vars$Datas[[titr]][[filenum]][,1:3])
		
		
		##### Binning function
		for (idx in 1:nbins){
			
			temp.low 	<- data.frame(d[[xy]][d[[xy]] > Bins[[seqxy]][idx]], d$z[d[[xy]] > Bins[[seqxy]][idx]])
			names(temp.low)	<- c(xy, "z")
			final.bin		<- data.frame(temp.low[[xy]][temp.low[[xy]]<= Bins[[seqxy]][idx+1]], temp.low$z[temp.low[[xy]]<= Bins[[seqxy]][idx+1]])
			names(final.bin)<- c(xy, "z")
			final.z 		<- final.bin$z
			finalpos		<- final.z[final.z>0]
			
			sizevect[idx] 	<- length(finalpos)
			meanvect[idx]	<- 10^mean(log10(finalpos))
			sigmavec[idx] 	<- ifelse((sizevect[idx]>2), sd(log10(finalpos)),0)
			                                  		
			Finalbin 		<<-finalpos		
			if (length(final.z)>Options$thresh){
				bin.idx			<- paste(xy, idx, ".z", filenum, sep="")
				idxvect[idx]	<- bin.idx
				Binz[["Full"]][[bin.idx]]	<- final.z
				
				if (Options$analysis$zchoice =="pct"){
					if(filenum==Vars$Control.File[[titr]]){
						# Sets cutoff for positive cells basd on control
						#show("Setting cutoffs")
						pos.cutoff <- sort(final.z)[floor(length(final.z)*(1-(Options$analysis$pctgate/100)))]
					}
					else {
						# Evaluates fraction of positive cells
						if (Bins[[xy]][["pos.cut"]][idx] > 0){
							pos.pctvect[idx]	<- length(final.z[final.z > Bins[[xy]][["pos.cut"]][idx]]) / length(final.z)*100
						}
						else {
							pos.pctvect[idx] <- 0
			}
		}}}}
		
		# Store statistics of bins in appropriate entry in Bins$sizes or Bins$means
		Binz[["means"]] 	<- meanvect
		Binz[["sizes"]] 	<- sizevect
		Binz[["pospct"]]	<- pos.pctvect
		Binz[["pos.cutoff"]]<- pos.cutoff
		Binz[["filenum"]]	<- filenum
		Binz[["sigmas"]]	<- sigmavec

		
		return(Binz)
	}
	
	moveData	<- function(xy,titr,dose,Bins,tempBins, file.idx){
		Bins[[xy]]$means[[titr]][[dose]] 	<- tempBins[[file.idx]][["means"]]
		Bins[[xy]]$sizes[[titr]][[dose]] 	<- tempBins[[file.idx]][["sizes"]]
		Bins[[xy]]$pospct[[titr]][[dose]]	<- tempBins[[file.idx]][["pospct"]]
		Bins[[xy]]$pos.cut[[titr]]			<- tempBins[[file.idx]][["pos.cutoff"]]
		Bins[[xy]]$sigmas[[titr]][[dose]]	<- tempBins[[file.idx]][["sigmas"]]
		return(Bins)
	}
	
	### Bin control files for x and y
	show("Binning Control file")
	ctl	<- Vars$Control.File[[titr]]
	tempBins	<- list()
	tempBins[[1]] <- Partition(ctl, titr, Bins, Vars, Options, "x")
	Bins <- moveData(xy="x",titr,dose=ctl,Bins,tempBins,1)   
	
	tempBins[[1]] <- Partition(ctl, titr, Bins, Vars, Options, "y")
	Bins <- moveData(xy="y",titr,dose=ctl,Bins,tempBins,1)	                                                          
	 
	ScsBins	<<- Bins
	
	filevec	<- c(1:Vars$num.doses)
	treated.filenums	<- filevec[-Vars$Control.File[[titr]]]

	#### Bin x
	parallelBins	<- foreach (filenum=treated.filenums) %do% { 
		Partition(filenum, titr, Bins, Vars, Options, "x")
	}
	#show(length(parallelBins))
	for (file.idx in 1:length(parallelBins)) {
		#show(file.idx)
		filenum	<- parallelBins[[file.idx]][["filenum"]]
		Bins <- moveData(xy="x", titr, dose=filenum, Bins, tempBins=parallelBins, file.idx)
	}                       
	
	#### Bin Y:
	parallelBins	<- foreach (filenum=treated.filenums) %do% { 
		Partition(filenum, titr, Bins, Vars, Options, "y")
	}
	for (file.idx in 1:length(parallelBins)) {
		filenum	<- parallelBins[[file.idx]][["filenum"]]
		Bins <- moveData(xy="y", titr, dose=filenum, Bins, tempBins=parallelBins, file.idx)
	}
	ScsBins <<- Bins	     
	stoptime <- proc.time()[1]
	show(sprintf("Binning time was %.2f seconds",stoptime-starttime))
	return(Bins)
}
