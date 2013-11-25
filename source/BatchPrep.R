BatchPrep	<- function(Vars,Track,Options,Bins){
	
	show("BatchPrep")
	t4 <- tktoplevel()
	#show(Bins)
		
	num.doses 	<- Vars$num.doses
	num.titr	<- Vars$num.titr
	Datas		<- Vars$Datas
	
	# Scaling for window
	myhscale <- .85
	myvscale <- 1
	
	smallfile <- F
	
	current.params	<- paste(	Options$bins$min$x,Options$bins$min$y,Options$bins$max$x,
									Options$bins$max$y,Options$bins$n$x,Options$bins$n$y, 
									Vars$islog$x, Vars$islog$y, Vars$islog$z,
									
									sep="_")
	
									
	if ((Track$FirstInitialPlot == T) | any(Bins$params != current.params)) {
		# Set the binning and display parameters according to data
		
		### Set binning variables to correct size
		length(Bins$sizes) 		<- Vars$num.titr
		length(Bins$means) 		<- Vars$num.titr
		
		### Set up lists to receive binning data.
		### Determine maximum and minimum of X Y Z across all Datas.
		
		# Used for binning and plotting.
		# The values for x and y are log10 transformed if the "log" box was checked.
		
		maxx	<- 0
		maxy	<- 0
		maxz	<- 0                                     
		minx	<- 1e6
		miny	<- 1e6
		minz	<- 1e6        
		
		statmatrix	<- matrix(rep(0, num.doses*num.titr),nrow=num.doses, ncol=num.titr)
		z1ths	<- statmatrix
		z99ths	<- statmatrix
		x1ths	<- statmatrix
		x99ths	<- statmatrix
		y1ths	<- statmatrix
		y99ths	<- statmatrix
		samplesizes	<- statmatrix
		
		for (titr.idx in 1:num.titr){
			Bins$sizes[[titr.idx]] 	<- list()
			Bins$means[[titr.idx]] 	<- list()
			Bins$pospct[[titr.idx]] <- list()
			Bins$sigmas[[titr.idx]] <- list()
			Bins$Full[[titr.idx]] 			<- list()
			Bins$pos.cut.mat[[titr.idx]] 	<- matrix()
			
			
			#2D variables
			Bins$x[[titr.idx]]	<- list()
			Bins$x[["Full"]][[titr.idx]]	<- list()
			Bins$x[["means"]][[titr.idx]]	<- list()
			Bins$x[["sizes"]][[titr.idx]]	<- list()
			Bins$x[["pospct"]][[titr.idx]]	<- list()
			Bins$x$sigmas[[titr.idx]]		<- list()
			
			
			for (dose.idx in 1:num.doses){
				
				#3D binning variables sizing
				Bins$sizes[[titr.idx]][[dose.idx]]		<- matrix(0)
				Bins$means[[titr.idx]][[dose.idx]]		<- matrix(0)
				Bins$pospct[[titr.idx]][[dose.idx]]		<- matrix(0)
				Bins$Full[[titr.idx]][[dose.idx]]		<- list()
				Bins$pos.cut.mat[[titr.idx]][[dose.idx]]<- matrix(0)
				Bins$sigmas[[titr.idx]][[dose.idx]]		<- matrix(0)
				
				#2D binning
				Bins$x[["Full"]][[titr.idx]][[dose.idx]]	<- c(0)
				Bins$x[["means"]][[titr.idx]][[dose.idx]]	<- c(0)
				Bins$x[["sizes"]][[titr.idx]][[dose.idx]]	<- c(0)
				Bins$x[["pospct"]][[titr.idx]][[dose.idx]]	<- c(0)
				Bins$x[["sigmas"]][[titr.idx]][[dose.idx]]		<- c(0)
				# Stats of data files
				
				maxx	<- max(c(Datas[[titr.idx]][[dose.idx]][,1],maxx))
				maxy	<- max(c(Datas[[titr.idx]][[dose.idx]][,2],maxy))
				maxz	<- max(c(Datas[[titr.idx]][[dose.idx]][,3],maxz))
				
				minx	<- max(0,min(c(Datas[[titr.idx]][[dose.idx]][,1],maxx)))
				miny	<- max(0,min(c(Datas[[titr.idx]][[dose.idx]][,2],miny)))
				minz	<- max(0,min(c(Datas[[titr.idx]][[dose.idx]][,3],minz)))
				
				# Determine 99th and 1st percentile of data to set binning and display ranges
				
				x1ths[dose.idx,titr.idx]	<- quantile(Datas[[titr.idx]][[dose.idx]][,1], 0.001)
				y1ths[dose.idx,titr.idx]	<- quantile(Datas[[titr.idx]][[dose.idx]][,2], 0.001)
				z1ths[dose.idx,titr.idx]	<- quantile(Datas[[titr.idx]][[dose.idx]][,3], 0.001)
				
				x99ths[dose.idx,titr.idx]   <- quantile(Datas[[titr.idx]][[dose.idx]][,1], 0.999)
				y99ths[dose.idx,titr.idx]   <- quantile(Datas[[titr.idx]][[dose.idx]][,2], 0.999)
				z99ths[dose.idx,titr.idx]   <- quantile(Datas[[titr.idx]][[dose.idx]][,3], 0.999)
				
				# evaluate number of cells for binning
				samplesizes[dose.idx,titr.idx] <- length(Datas[[titr.idx]][[dose.idx]][,1])
			}
		}
		Bins[["y"]] <- Bins$x
		
		Vars$stats$maxx <- maxx
		Vars$stats$maxy <- maxy                     
		Vars$stats$maxz <- maxz
		Vars$stats$minx	<- minx
		Vars$stats$miny	<- miny
		Vars$stats$minz	<- minz
		
		## Z values for plotting
		zcutoff	<- max(z99ths)
		Options$plots$zcutoff <- zcutoff
		
		#Options$plots$initialZmin	<- Bins$binsmin
		#Options$plots$initialZmax	<- Bins$binsmax
		
		Vars$stats$zrange		<- Vars$stats$maxz - Vars$stats$minz
		Vars$stats$samplesize	<- median(samplesizes)
		
		ScsVars <<- Vars
		
	
		
		if(Vars$islog$x){
			# Set x ranges for log scale
			Options$bins$min$x 	<- round(max(min(x1ths)-.2, minx),2)
			Options$bins$max$x 	<- round(min(max(x99ths)+0.2,maxx),2)   
			
			xrange	<- Options$bins$max$x - Options$bins$min$x    
			Options$bins$n$x	<- floor(xrange*5)
			
			max.xy.x	<- Options$bins$max$x
		
		}else{ 
			# Set x ranges for linear scale
			Options$bins$min$x 	<- round(max(1, floor(min(x1ths)*0.8))) # a little less minimum 0.1%ile, being sure it's positive  
			Options$bins$max$x 	<- round(min(max(x99ths)*1.2, maxx))   # a little more than the 99.9%ile, but less than the max      
			
			xrange	<- Options$bins$max$x - Options$bins$min$x
			Options$bins$n$x	<- floor(log10(xrange)*5)
			max.xy.x	<- log10(Options$bins$max$x)
		
		}
		
		if(Vars$islog$y){          
			# Set y range in log scale
			Options$bins$min$y 	<- round(max(min(y1ths)-.2,  miny),2)              
			Options$bins$max$y 	<- round(min(max(y99ths)+0.2,maxy),2)
			
			yrange				<- Options$bins$max$y - Options$bins$min$y
			Options$bins$n$y	<- floor(yrange*5)
			max.xy.y			<- Options$bins$max$y
		
		}else{
			# Set y ranges for linear scale
			Options$bins$min$y 	<- round(max(1, floor(min(y1ths)*0.8))) # a little less minimum 0.1%ile, being sure it's positive             
			Options$bins$max$y 	<- round(min(max(y99ths)*1.2, maxy))# a little more than the 99.9%ile, but less than the max     
			
			yrange				<- Options$bins$max$y - Options$bins$min$y
			Options$bins$n$y	<- floor(log10(yrange)*5)
			max.xy.y			<- log10(Options$bins$max$y)
		}
		
		
		# Square up the axes
		max.xy.square <- max(max.xy.x, max.xy.y)
		if (max.xy.square > 4) {
			Options$plots$max$x <- ifelse(Vars$islog$x, max.xy.square, 10^max.xy.square)
			Options$plots$max$y <- ifelse(Vars$islog$y, max.xy.square, 10^max.xy.square)  
		}
		
		ScsOptions <<- Options
		
		Track$FirstInitialPlot <- F
	}
	
	
	
	
	
	#### BIN SPECIFICATION #####################################################
	# This function takes as input a tcl value of bin number  
	# and creates variables for number and size of bins, 
	# with edges and centers of all bins as vectors.  Also contains lists
	# to hold stats about bins and bin data.
	BinSpec <- function(Options, Vars){
		#show("Specifying bins (BinSpec)")
		#show(sprintf("length of Bins$sizes=%i", length(Bins$sizes)))
		Bins$nx		<<- Options$bins$n$x
		Bins$ny		<<- Options$bins$n$y
		Bins$seqx 	<<- seq(Options$bins$min$x, Options$bins$max$x, length.out=(Bins$nx+1))
		Bins$sizex 	<<- (Bins$seqx[2]-Bins$seqx[1])
		Bins$centersx 	<<- Bins$seqx[1:Bins$nx]+0.5*Bins$sizex
		Bins$seqy 		<<- seq(Options$bins$min$y, Options$bins$max$y, length.out=(Bins$ny+1))
		Bins$sizey 		<<- (Bins$seqy[2]-Bins$seqy[1])
		Bins$centersy 	<<- Bins$seqy[1:Bins$ny]+0.5*Bins$sizey
		ScsBins 	<<- Bins 
	}
	Bins <- BinSpec(Options, Vars)
	#show(sprintf("%d X, %d Y bins, %d total", Bins$nx, Bins$ny, Bins$nx*Bins$ny))
	
	
	Bins <- BinPartition3d(Bins, Vars, Options, titr=1)
	show("Binning complete")
	ScsBins <<- Bins
	
	#Use bin values to set upper and lower display limits for bin coloring
	tempMFI		<- unlist(Bins$means)
	tempSIZE	<- unlist(Bins$sizes)
	tempMFI <- tempMFI[tempSIZE>5]     
	tempMFI	<- tempMFI[tempMFI>0]
	density.cutoff <- 0.05
	Options$plots$initialZmin 	<- quantile(tempMFI, density.cutoff, na.rm=T)
	Options$plots$initialZmax	<- quantile(tempMFI, (1-density.cutoff), na.rm=T)
	
	ScsOptions <<- Options
	
	
	#### NAVIGATION ############################################################
	# Store data to plot first file.  Stores
	# x and y values in temporary variables, to be overwritten
	# as the user moves through files.
	dose.current 	<- 1
	titr.current	<- 1
	ptstemp <- Datas[[titr.current]][[dose.current]][,1:3]
		
	# Move to next or previous file in list
	NextFile <- function(){
		if (dose.current < num.doses) {
			dose.current <<- dose.current+1
		}else{
			dose.current <<- 1
			if(titr.current < num.titr){
				titr.current <<- titr.current+1
				Bins <<- BinPartition3d(Bins, Vars, Options, titr.current)
				
			}else{
				titr.current <<- 1
				Bins <<- BinPartition3d(Bins, Vars, Options, titr.current)
			}
		}                        
		ptstemp 	<<- Datas[[titr.current]][[dose.current]][,1:3]
		
		tkrreplot(initialplot.img)
		tkbind(initialplot.img, "<Button-1>",OnLeftClickInit)
	}
	
	PrevFile <- function(){
		if (dose.current > 1){
			dose.current <<- dose.current-1
		}else{
			dose.current <<- num.doses
			if(titr.current > 1){
				titr.current<- titr.current-1
				Bins <<- BinPartition3d(Bins, Vars, Options, titr.current)
			}else{
				titr.current <- num.titr
				Bins <<- BinPartition3d(Bins, Vars, Options, titr.current)
			}
		}
		ptstemp 	<<- Datas[[titr.current]][[dose.current]][,1:3]
		
		tkrreplot(initialplot.img)
		tkbind(initialplot.img, "<Button-1>",OnLeftClickInit)
	}
	
	
	# Goes back to previous step (Choosing and labeling X Y and Z)
	OnBack <- function(){
		tkdestroy(t4)
		DataPrep(Vars,Track,Options,Bins,titr.current=1)
	}
	
	
	
	
	
	########################################################################################
	#### PLOTTING #####
	
	
	parPlotSize <- c()
	usrCoords 	<- c()
	
	Plot <- function(){
		#show(c("Plotting"))
		starttime <- proc.time()[1]
		
		## Set graphical parameters for a square plot.
		
		oldpar <- par(no.readonly = T)
		parplt1 <- par("plt")[1]
		parplt2	<- 0.99
		parplt3 <- par("plt")[3]
		parplt4 <- par("plt")[4]
		par("plt"=c(parplt1, parplt2, parplt3, parplt4))
		ysize 	<- parplt4-parplt3
		xsize 	<- parplt2-parplt1
		myhscale <- (ysize/xsize)
		

		
		### PLOTTING
		# Set axis limits
		xlim	<- c(Options$plots$min$x, Options$plots$max$x)
		ylim	<- c(Options$plots$min$y, Options$plots$max$y)
		
		params <- par(bg="white")
		
		plot.title <- paste(Vars$prefix[[titr.current]], ":", Vars$Concs[[titr.current]][[dose.current]], Vars$units[[titr.current]])
		
		plot(xlim, ylim, pch=".", col="white", 
				main = plot.title, xlim=xlim, ylim=ylim,
				xlab=NA, ylab=NA, type="n", xaxt="n", yaxt="n"
				)
		
		niceaxes(xy="x", limits=xlim, axtitle=Vars$xyzLabels[[titr.current]]$x, islog=Vars$islog$x)
		niceaxes(xy="y", limits=ylim, axtitle=Vars$xyzLabels[[titr.current]]$y, islog=Vars$islog$y)
		
		
		# Plot the points.
		if(smallfile == F){
			points(ptstemp[,1], ptstemp[,2], pch=".")
		}
		
		
		# Draw outlines of bins        
		for (x in Bins$seqx){lines(c(x, x), c(Options$bins$min$y, Options$bins$max$y), col=colors()[45])}
		for (y in Bins$seqy){lines(c(Options$bins$min$x, Options$bins$max$x), c(y, y), col=colors()[45])}
		
		#show(sprintf("Drawing %i x and %i y bins", Bins$nx, Bins$ny))
		# Analyze bin occupancy and draw a rectangle.  Rect edge is red if bin 
		# has above the specified number of cells.  Can be colored by Z MFI
		# or number of positive gated events (set in options).
		if(Options$plots$initialZmin==0){
			minval <- 0
		}else{
			minval		<- sign(Options$plots$initialZmin)*log10(abs(Options$plots$initialZmin))
		}
		maxval 		<- sign(Options$plots$initialZmax)*log10(abs(Options$plots$initialZmax))
		binrange 	<- maxval - minval                                      
	
		
		for (x.idx in 1:Bins$nx){
			for (y.idx in 1:Bins$ny){
				binfill <- Bins$sizes[[titr.current]][[dose.current]][x.idx,y.idx]
				bin.fill.color <- NA
				#show(sprintf("bin %i,%i has %i cells",x.idx, y.idx, binfill))
				if (binfill >= Options$thresh){
					if (Options$plots$showz) {
						if (Options$analysis$zchoice == "mfi"){
							zval	<- Bins$means[[titr.current]][[dose.current]][x.idx,y.idx]
							if(!is.nan(zval)){
								zval <- sign(zval)*log10(abs(zval))
								if(zval > minval & zval < maxval){
									#show(sprintf("zval = %.2f", zval))
									colr.ratio	<- (zval-minval) / binrange
									#show(sprintf("colr.ratio = %.2f", colr.ratio))
									colr.jet	<- colr.ratio*32
									colr.jet	<- ceiling(colr.jet)
									bin.fill.color	<- Jet[colr.jet]
								}
								else{
									if(zval>=maxval){
										bin.fill.color <- Jet[32]
										#show("above range")
									}
									if(zval<=minval){
										bin.fill.color <- Jet[1]
										#show("below range")
									}
								}
							}
						}
						if (Options$analysis$zchoice == "pct") {
							pos.percent	<- Bins$pospct[[titr.current]][[dose.current]][x.idx, y.idx]
							if(!is.na(pos.percent)){
								if(pos.percent==0){
									bin.fill.color	<- Jet[1]	
								}
								if(pos.percent>0) {
									colr.jet	<- pos.percent/100*32
									colr.jet	<- ceiling(colr.jet)
									bin.fill.color	<- Jet[colr.jet]
								}
							}
						}
						
					}
					#show("plotting rectangle")
					rect(Bins$seqx[x.idx], Bins$seqy[y.idx], Bins$seqx[x.idx+1], Bins$seqy[y.idx+1], col=bin.fill.color, border="red")
				}
			}
		}
		
		parPlotSize <<- par("plt")
		usrCoords   <<- par("usr")
		par(oldpar)
		
		stoptime <- proc.time()[1]
		#show(sprintf("plotting time was %.2f seconds",stoptime-starttime))
		
	}
	
	
	#########################################################################################
	### COLORBAR ###
	PlotColorbar	<- function() {
		if (Options$plots$showz) {
			
			if (Options$analysis$zchoice == "mfi") {
				Colorbar(bottom=Options$plots$initialZmin, top=Options$plots$initialZmax, titlestring="MFI", islog=Vars$islog$z)
			}
			if (Options$analysis$zchoice == "pct") {
				Colorbar(bottom=0, top=100, titlestring="% positive", islog=F)
			}
		}
		else {  #if showz=F
			plot(c(0,1),c(0,1), pch="", xlab=NA, ylab=NA, xaxt="n", yaxt="n") 
		}
	}
	
	
	
	
	##################################################################################################################
	#### HISTOGRAM PLOTTING ##########################################################################################	
	OnLeftClickInit <- function(x,y){
		# Will find the point in a list closest to the user's click.
		# Mostly copied directly from example by Wettenhall
		
		xClick <- x
		yClick <- y
		width  <- as.numeric(tclvalue(tkwinfo("reqwidth",initialplot.img)))
		height <- as.numeric(tclvalue(tkwinfo("reqheight",initialplot.img)))
		xMin <- parPlotSize[1] * width
		xMax <- parPlotSize[2] * width
		yMin <- parPlotSize[3] * height
		yMax <- parPlotSize[4] * height
		rangeX <- usrCoords[2] - usrCoords[1]
		rangeY <- usrCoords[4] - usrCoords[3]
		
		xClick <- as.numeric(xClick)+0.5
		yClick <- as.numeric(yClick)+0.5
		yClick <- height - yClick
		
		##X-Y values of click
		xPlotCoord <- usrCoords[1]+(xClick-xMin)*rangeX/(xMax-xMin)
		yPlotCoord <- usrCoords[3]+(yClick-yMin)*rangeY/(yMax-yMin)
		#show(sprintf("clicked at (%.2f,%.2f)",xPlotCoord,yPlotCoord))
		
		if (xPlotCoord>=Options$bins$min$x && xPlotCoord<=Options$bins$max$x && yPlotCoord>=Options$bins$min$y && yPlotCoord<=Options$bins$max$y){
			x.idx <- ceiling((xPlotCoord-Options$bins$min$x)/(Options$bins$max$x-Options$bins$min$x) * Bins$nx)
			y.idx <- ceiling((yPlotCoord-Options$bins$min$y)/(Options$bins$max$y-Options$bins$min$y) * Bins$ny)
			#show(sprintf("x.idx = %i; y.idx=%i",x.idx,y.idx))
			points.idx	<<- paste(x.idx,",",y.idx, sep="")
			x.idx <<- x.idx
			y.idx <<- y.idx
		}
		else { points.idx <<- "outside fit area" }
		HistPlot(Bins,x.idx,y.idx)
		
	}
	
	HistPlot<- function(Bins, x.idx,y.idx){
		#show("Plotting histogram")
		hist.idx	<- sprintf("%i,%i", x.idx, y.idx)
		#show(sprintf("hist.idx is %s, which contains %i cells", hist.idx, length(Bins$Full[[titr.current]][[dose.current]][[hist.idx]])))
		
		if(length(Bins$Full[[titr.current]][[dose.current]][[hist.idx]]) > 0){
			histdata 	<- Bins$Full[[titr.current]][[dose.current]][[hist.idx]]
			t4.hist <- tktoplevel()
			tkwm.title(t4.hist, paste("Histogram of", Vars$xyzLabels[[titr.current]]$z, "values"))
			
			
			gatemin <- Bins[["pos.cut.mat"]][[titr.current]][x.idx,y.idx]
			if (Vars$islog$z){
				histdata	<- log10(histdata[histdata >= 1])
				zcutoff		<- log10(Vars$stats$maxz)
				gatemin		<- log10(gatemin)
				
			}
			else {
				histdata	<- histdata[histdata > Options$plots$minz]
				histdata	<- histdata[histdata < zcutoff]  #zcutoff defined above as 99.9th percentile
			}
			
			
			PlotHist<- function(){
				par(bg="white")
				histo <- hist(histdata, breaks=seq(0, zcutoff, length.out=30))
				#xticks	<- pretty(c(0,zcutoff))
				tkclipboard.clear()
				plot(histo, xlab=paste("log10", Vars$xyzLabels[[titr.current]][["z"]]), ylab="number of cells", main=paste("X bin =", x.idx, "Y bin =", y.idx))
				if (Options$analysis$zchoice == "pct"){
					# Draw gate for positive cells
					lines(c(gatemin, zcutoff), c(max(histo$counts)*.75, max(histo$counts)*.75), col="red")
					lines(c(gatemin, gatemin), c(max(histo$counts)*0.8, max(histo$counts)*0.7), col="red")
					lines(c(zcutoff, zcutoff), c(max(histo$counts)*0.8, max(histo$counts)*0.7), col="red")
					text(x=mean(c(gatemin, zcutoff)), y=max(histo$counts)*.75, label="+", pos=3, cex=1.5, col="red")
				}
				
			}
			
			#Histogram options
			
			OnClose		<- function(){tkdestroy(t4.hist)}
			# Window objects
			hist.img	<- tkrplot(t4.hist, fun=PlotHist, hscale=0.7, vscale=0.7)
			close.but	<- tkbutton(t4.hist, text=" Close ", command=OnClose)
			tkgrid(hist.img, columnspan=2)
			tkgrid(tklabel(t4.hist, text=" "))
			tkgrid(close.but)
		}
	}
	####
	
	
	
	
	########################################################################################################################
	#### OPTIONS ###########################################################################################################
	
	OnOptions <- function(){
		t4.options <- tktoplevel()
		tkwm.title(t4.options, "Initial Plot Options")
		
		Set.t4.options <- function(){
			
			# Ensures that limits of binning aren't outside of limits of plotting
			for (xy in c("x", "y")){
				for (minmax in c("min", "max")) {
					idx	<- paste(minmax, xy, sep="")
					Options$bins[[minmax]][[xy]]	<<- as.numeric(tclvalue(binminmax.tk[[idx]]))
					Options$plots[[minmax]][[xy]]	<<- as.numeric(tclvalue(plotminmax.tk[[idx]]))
				}
				if (Options$bins$min$x < Options$plots$min$x) {Options$plots$min$x <<- Options$bins$min$x}
				if (Options$bins$min$y < Options$plots$min$y) {Options$plots$min$y <<- Options$bins$min$y}
				if (Options$bins$max$x > Options$plots$max$x) {Options$plots$max$x <<- Options$bins$max$x}
				if (Options$bins$max$y > Options$plots$max$y) {Options$plots$max$y <<- Options$bins$max$y}
				
			}
			
			Options$thresh				<<- as.numeric(tclvalue(thresh.tcl))
			Options$bins$n$x			<<- as.numeric(tclvalue(nxbins.tcl))
			Options$bins$n$y			<<- as.numeric(tclvalue(nybins.tcl))
			Options$analysis$zchoice 	<<- as.character(tclvalue(zchoice.tk))
			Options$analysis$pctgate 	<<- as.numeric(tclvalue(pctgate.tk))
			Options$analysis$hillchoice	<<- as.character(tclvalue(hillchoice.tk))
			Options$plots$showz			<<- as.logical(as.numeric(tclvalue(showz.tk)))
			Options$plots$initialZmin	<<- as.numeric(tclvalue(initialZmin.tcl))
			Options$plots$initialZmax	<<- as.numeric(tclvalue(initialZmax.tcl))
			Bins$maxx	<<- Options$bins$max$x
			Bins$maxy	<<- Options$bins$max$y
			
			Bins <<- BinSpec(Options, Vars)
			Bins <<- BinPartition3d(Bins, Vars, Options, titr.current)
			ScsBins <<- Bins
			ScsOptions <<- Options

			tkrreplot(initialplot.img)
			tkrreplot(colorbar.img)
		}
		
		OnCloseOptions	<- function(){
			tkdestroy(t4.options)
		}
		
		OnOkOptions		<- function(){
			Set.t4.options()
			tkdestroy(t4.options)
		}                                                           
		#show("Options functions specified")
		
		
		
		### GET VALUES to fill in entries, etc. ################################
		xmax.tk 	<- tclVar(Options$plots$max$x)
		ymax.tk 	<- tclVar(Options$plots$max$y)
		zdec.tk 	<- tclVar(log10(Options$plots$maxz))
		
		binminmax.tk	<- list()
		binminmax.entry	<- list()
		plotminmax.tk	<- list()
		plotminmax.entry<- list()
		
		
		for (xy in c("x", "y")){            
			for (minmax in c("min", "max")) {
				idx	<- paste(minmax, xy, sep="")
				plotminmax.tk[[idx]]	<- tclVar(Options$plots[[minmax]][[xy]])
				plotminmax.entry[[idx]]	<- tkentry(t4.options, width="10", textvariable=plotminmax.tk[[idx]])
				binminmax.tk[[idx]]		<- tclVar(Options$bins[[minmax]][[xy]])     
				binminmax.entry[[idx]]	<- tkentry(t4.options, width="10", textvariable=binminmax.tk[[idx]])
			}
		}
		
		# Value of cell number cutoff for bins to analyze by hill fit
		thresh.tcl <- tclVar(Options$thresh)
		
		#Set initial number of bins
		nxbins.tcl <- tclVar(Options$bins$n$x)
		nybins.tcl <- tclVar(Options$bins$n$y)
		
		# Specifies initial choice for whether to #show z values                      
		showz.tk 	<- tclVar(Options$plots$showz)
		zchoice.tk	<- tclVar(Options$analysis$zchoice)
		pctgate.tk	<- tclVar(Options$analysis$pctgate)  
		
		# Values for colorbars for plotting
		#show(Options$plots)
		initialZmin.tcl	<- tclVar(floor(Options$plots$initialZmin))
		initialZmax.tcl	<- tclVar(ceiling(Options$plots$initialZmax)) 
		
		
		### LAYOUT #############################################################
		
		xminmax.label	<- tklabel(t4.options, text=" x-axis limits ")
		yminmax.label	<- tklabel(t4.options, text=" y-axis limits ")
		
		label.head	<- tklabel(t4.options, text="                   ")
		min.head	<- tklabel(t4.options, text=" Min ")
		max.head	<- tklabel(t4.options, text=" Max ")
		
		zdisplay.frm	<- tkframe(t4.options)
		zdisplayhead	<- tklabel(zdisplay.frm, text=" Colorbar range ")
		zmin.head	<- tklabel(zdisplay.frm, text=" Min ")
		zmax.head	<- tklabel(zdisplay.frm, text=" Max ")
		zmin.entry	<- tkentry(zdisplay.frm, width="10", textvariable=initialZmin.tcl)
		zmax.entry	<- tkentry(zdisplay.frm, width="10", textvariable=initialZmax.tcl)
		tkgrid(zdisplayhead)                                      
		tkgrid(zmin.head, zmax.head)
		tkgrid(zmin.entry, zmax.entry)
		
		
		## Plotting options layout
		tkgrid(tklabel(t4.options, text="Plotting Options ")) #Row1
		tkgrid(label.head, min.head, max.head) #2
		tkgrid(xminmax.label, plotminmax.entry$minx, plotminmax.entry$maxx) #3
		tkgrid(yminmax.label, plotminmax.entry$miny, plotminmax.entry$maxy) #4
		tkgrid(tklabel(t4.options, text="  "))
		
		## Binning options layout
		tkgrid(tklabel(t4.options, text="Binning Options"))
		tkgrid(tklabel(t4.options, text="Min"), row=5, column=1)
		tkgrid(tklabel(t4.options, text="Max"), row=5, column=2)
		tkgrid(tklabel(t4.options, text="X Bin Limits "), binminmax.entry$minx, binminmax.entry$maxx)
		tkgrid(tklabel(t4.options, text="Y Bin Limits "), binminmax.entry$miny, binminmax.entry$maxy)
		
		nxbins.entry	<- tkentry(t4.options, width="6", textvariable=nxbins.tcl)
		nybins.entry	<- tkentry(t4.options, width="6", textvariable=nybins.tcl)
		thresh.entry	<- tkentry(t4.options, width="6", textvariable=thresh.tcl)
		label.nbins	<- tklabel(t4.options, text="Number of bins")
		
		label.thresh	<- tklabel(t4.options, text="Min. cells per bin")
		x.label <- tklabel(t4.options, text="  X  ")
		y.label <- tklabel(t4.options, text="  Y  ")
		tkgrid(tklabel(t4.options, text=" "), x.label, y.label)
		tkgrid.configure(y.label, sticky = "w")
		tkgrid(label.nbins, nxbins.entry, nybins.entry)
		tkgrid.configure(label.nbins, sticky="e")
		tkgrid.configure(nybins.entry, sticky="w")
		tkgrid(tklabel(t4.options, text="               "))
		tkgrid(label.thresh, thresh.entry)
		tkgrid.configure(label.thresh, sticky="e")    
		tkgrid(tklabel(t4.options, text=" "))
		
		
		tkgrid(tklabel(t4.options, text="Display Options"))
		display.fr		<- tkframe(t4.options)
		occupied.label 	<- tklabel(display.fr, text="Occupancy")
		withzvals.label	<- tklabel(display.fr, text="with Z values")
		occupied.rb		<- tkradiobutton(display.fr)
		zval.rb			<- tkradiobutton(display.fr) 
		
		tkconfigure(occupied.rb, variable=showz.tk, value=0)
		tkconfigure(zval.rb, 	 variable=showz.tk, value=1)
		tkgrid(occupied.rb, occupied.label, tklabel(t4.options, text="   "), zval.rb, withzvals.label)
		tkgrid.configure(occupied.label, sticky="w")
		tkgrid.configure(occupied.rb, sticky="e")
		tkgrid.configure(withzvals.label, zval.rb, sticky="w")
		tkgrid.configure(zval.rb, sticky="e")
		tkgrid(display.fr, columnspan=4, sticky="w")
		tkgrid(tklabel(t4.options, text=" "))
		tkgrid(zdisplay.frm, columnspan=3, sticky="w")
		tkgrid(tklabel(t4.options, text=" "))
		                                        
		## Analysis Options
		# Options to choose the numbers that will be used for fitting of
		# hill functions.  MFI will use mean fluorescence intensity, whereas
		# pct will set a gate  that includes only the top  
		# (user specified)% of cells in the control, then analyzes how many cells
		# are in that gate in other samples.  
		                                                        
		tkgrid(tklabel(t4.options, text="Analysis Options"))
		mfipct.fr	<- tkframe(t4.options)
		mfi.rb		<- tkradiobutton(mfipct.fr)
		pct.rb		<- tkradiobutton(mfipct.fr)
		tkconfigure(mfi.rb, variable=zchoice.tk, value="mfi")
		tkconfigure(pct.rb, variable=zchoice.tk, value="pct")
		pct.entry	<- tkentry(mfipct.fr, width="3", textvariable=pctgate.tk)
		mfi.label	<- tklabel(mfipct.fr, text="MFI")                      
		pct1.label	<- tklabel(mfipct.fr, text="Gate top ")
		pct2.label	<- tklabel(mfipct.fr, text="% of control")
		spacer	<- tklabel(mfipct.fr, text="   ")
		tkgrid(mfi.rb, mfi.label, spacer, pct.rb, pct1.label,pct.entry, pct2.label)
		tkgrid(mfipct.fr, columnspan=4, sticky="w")
		                            
		# Hill coefficient options
		hillchoice.tk	<- tclVar(Options$analysis$hillchoice)
		hill.frm	<- tkframe(t4.options)
		hill.label	<- tklabel(hill.frm, text="Hill coeff: ")
		hill.fixd.label	<- tklabel(hill.frm, text="Fixed (=1)   ")
		hill.free.label	<- tklabel(hill.frm, text="Free")
		hill.free.rb	<- tkradiobutton(hill.frm)
		hill.fixd.rb	<- tkradiobutton(hill.frm)
		tkconfigure(hill.free.rb, variable=hillchoice.tk, value="free")
		tkconfigure(hill.fixd.rb, variable=hillchoice.tk, value="fixed")
		tkgrid(hill.label, hill.fixd.rb, hill.fixd.label, hill.free.rb, hill.free.label)
		tkgrid(hill.frm, columnspan=4, sticky="w")
		
		                                             
		## Navigation ###
		apply.options.but	<- tkbutton(t4.options, text="  Update ", command=Set.t4.options)
		OK.options.but		<- tkbutton(t4.options, text="   OK    ", command=OnOkOptions)
		close.options.but	<- tkbutton(t4.options, text="  Close  ", command=OnCloseOptions)
		                                                                         
		tkgrid(tklabel(t4.options, text=" "))
		tkgrid(close.options.but, apply.options.but, OK.options.but)
		tkgrid(tklabel(t4.options, text=" "))
		
	}
	###############################
	
	# Goes to next step, doing Hill fits of all data
	OnProcess <- function() {
		
		for(titr.idx in 1:Vars$num.titr){
			Bins	<- BinPartition3d(Bins,Vars,Options, titr.idx)
		}
		
		# These lists are stored as global variables for debugging and are
		# updated at each transition.
		ScsOptions 	<<- Options
		ScsVars		<<- Vars
		ScsTrack	<<- Track
		ScsBins		<<- Bins
		
		Bins$Full	<- NULL
		gc()
		
		
		
		Answers <- FinalAnalyze(Vars,Options,Bins)
		Analysis3D(Answers,Vars,Options,Bins,Track)
	}
	
	
	OnProcess2D	<- function() {
		
		# Store global variables for debugging
		ScsOptions 	<<- Options
		ScsVars		<<- Vars
		ScsTrack	<<- Track
		ScsBins		<<- Bins
		
		# Bin in 2D and analyze
		Plots2d(Vars,Options,Bins)
	}
	
	
	##### EXPORT ###################
	
	OnPreExport <- function(){
		initialplots.export <- tktoplevel()
		tkwm.title(initialplots.export, "Export Options")
		exportprefix.tk	<- tclVar("Experiment 1")
		prefix.label 	<- tklabel(initialplots.export, text=" File Name: ")
		prefix.entry	<- tkentry(initialplots.export, width=20, textvariable=exportprefix.tk)
		example.label	<- tklabel(initialplots.export, text=" Example: ")
		example.text	<- tklabel(initialplots.export, textvariable=tclVar(sprintf("%s.pdf",tclvalue(exportprefix.tk))))
		
		OnExport <- function(){
			exportfolder <- tkchooseDirectory(initialdir=Filefolder)
			exportfolder <- tclvalue(exportfolder)
			smallfile	<<- T
			if(!nchar(exportfolder)){
			tkmessageBox(message="No folder selected")}
			else{
				filename <- paste(exportfolder,"/", tclvalue(exportprefix.tk), ".pdf", sep="")
				pdf(filename)
				
				for(titr.idx in 1:num.titr){
					for(file.idx in 1:num.doses){
						dose.current <<- file.idx
						ptstemp 	<<- Datas[[titr.current]][[dose.current]][,1:3]
						Plot()
					}
				}
				PlotColorbar()
				dev.off()
			}
			smallfile	<<- F
			tkdestroy(initialplots.export)
		}	
		
		OnCancelExport		<- function(){tkdestroy(initialplots.export)}
		
		export.cancel.but	<-tkbutton(initialplots.export, text=" Cancel ", command=OnCancelExport)
		final.export.but	<- tkbutton(initialplots.export, text=" Export ", command=OnExport)
		
		tkgrid(prefix.label, prefix.entry)
		tkgrid.configure(prefix.entry, sticky="ew")
		tkgrid.configure(prefix.label, sticky="e")
		tkgrid(example.label, example.text)
		tkgrid.configure(example.label, sticky="e")
		tkgrid.configure(example.text, sticky="w")
		tkgrid(tklabel(initialplots.export, text="  "))
		tkgrid(export.cancel.but, final.export.but, tklabel(initialplots.export, text="  "))
		tkgrid.columnconfigure(initialplots.export, 1, weight=1)
		tkbind(initialplots.export, "<Return>", OnExport)
	}
	
	
	
	####  LAYOUT  ########################################################
	
	tkwm.title(t4, "Initial Plot for Binning Analysis")
	
	plotwindow.frm	<- tkframe(t4)
	initialplot.img <- tkrplot(plotwindow.frm, fun=Plot, hscale=myhscale)
	colorbar.img	<- tkrplot(plotwindow.frm, fun=PlotColorbar, hscale=.25)
	tkgrid(initialplot.img, colorbar.img)
	tkgrid(plotwindow.frm)
	tkbind(initialplot.img, "<Button-1>", OnLeftClickInit)
	tkconfigure(initialplot.img, cursor="hand2")
	
	fwdrev.frm	<- tkframe(t4)
	fwd.but 	<- tkbutton(fwdrev.frm, text=">>", command=NextFile)
	rev.but 	<- tkbutton(fwdrev.frm, text="<<", command=PrevFile)
	fwdrev.label<- tklabel(fwdrev.frm, text=" File ")
	tkgrid(rev.but,fwdrev.label,fwd.but)
	tkgrid.configure(rev.but, sticky="e")
	tkgrid.configure(fwd.but, sticky="w")
	tkgrid(fwdrev.frm)
	
	sep <- ttkseparator(t4, orient="horizontal")
	tkgrid(sep, columnspan=1, sticky="ew")
	
	settings.frm	<- tkframe(t4)
	options.t4.but	<- tkbutton(settings.frm, text=" Options ", command=OnOptions)
	tkgrid(options.t4.but)
	tkgrid(settings.frm)
	
	sep <- ttkseparator(t4, orient="horizontal")
	tkgrid(sep, columnspan=1, sticky="ew")
	
	nav.frm		<- tkframe(t4)
	back.but		<- tkbutton(nav.frm, text="  Back  ", command=OnBack)
	export.but		<- tkbutton(nav.frm, text=" Export ", command=OnPreExport)
	analyze2d.but	<- tkbutton(nav.frm, text=" Analyze 2D  -> ", command=OnProcess2D)
	analyze3d.but	<- tkbutton(nav.frm, text=" Analyze 3D  -> ", command=OnProcess)
	tkgrid(analyze2d.but, column=2)
	tkgrid(back.but, export.but, analyze3d.but)
	tkgrid.configure(back.but, sticky="w")
	tkbind(t4, "<Return>", OnProcess)
	tkgrid(nav.frm)
}


