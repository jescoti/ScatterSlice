Analysis3D <- function(Answers,Vars,Options,Bins,Track){
	show("Analysis3D")
	
	num.doses 	<- Vars$num.doses
	num.titr	<- Vars$num.titr
	titr.current	<- 1
	if (!("hillfit" %in% names(Options$Analysis))){
		Options$analysis$method	<- ifelse((Vars$num.doses <= 4), "subtract", "hillfit")
	}
	
	# Initial values for upper and lower limits of colorbar values will come from input values or
	# maxima and minima of fit values, but are customizable through an options
	# window (t5.options)
	
	ec		<- list()
	amp		<- list()
	hill	<- list()
	bas		<- list()
	top		<- list()
	ecRE	<- list()
	ampRE	<- list()
	hillRE	<- list()
	
	limits			<- list() # list for setting limits across all titrations
	limits[["ea"]]	<- list() # sub-list for setting limits for each titration

	density.cutoff	<- 0.05
	
	for(titr.idx in 1:num.titr) {
		ec[[titr.idx]]		<- Answers[[titr.idx]]$ec
		amp[[titr.idx]]		<- Answers[[titr.idx]]$amp
		bas[[titr.idx]]		<- Answers[[titr.idx]]$bas
		top[[titr.idx]]		<- Answers[[titr.idx]]$bas+Answers[[titr.idx]]$amp
		ecRE[[titr.idx]]	<- Answers[[titr.idx]]$ecRE
		ampRE[[titr.idx]]	<- Answers[[titr.idx]]$ampRE
		if(Options$analysis$hillchoice =="free"){
			hill[[titr.idx]]	<- Answers[[titr.idx]]$hill
			hillRE[[titr.idx]]	<- Answers[[titr.idx]]$hillRE
		}else{
			hill[[titr.idx]]	<- rep(0,length(Answers[[titr.idx]]$hill))
			hillRE[[titr.idx]]	<- rep(0,length(Answers[[titr.idx]]$hill))
		}
	}
	
	#show(amp)
	
	allanswers	<- data.frame(	ec=  as.numeric(as.vector(sapply(ec, c))), # will produce NAs if negative, but will be taken care of later
								amp= as.numeric(as.vector(sapply(amp,c))),
								hill=as.numeric(as.vector(sapply(hill,c))),
								bas= as.numeric(as.vector(sapply(bas,c))),
								top= as.numeric(as.vector(sapply(top,c)))
							 )
	#allanswers	<- log10(allanswers)
	allanswers	<- na.omit(allanswers)
	ALLANS		<<- allanswers

	# Set Amp scaling regardless of method. Needed for future ability to display amplitude estimates.
	bottom	<- quantile(allanswers[["amp"]], density.cutoff  , na.rm=T)
	top		<- quantile(allanswers[["amp"]], 1-density.cutoff, na.rm=T)
	limits[["amp.low.tcl"]]		<- tclVar(signif(bottom,3))
	limits[["amp.high.tcl"]]	<- tclVar(signif(top,3))
	
	#show(c(bottom, top))
	#show("finished")
	
	# Set other scalings only if hill fitting was used
	if(Options$analysis$method =="hillfit"){
		#show("Setting plot scaling")
		
		if (sum(!is.na(allanswers$ec))>2){
			bottom	<- max(quantile(allanswers$ec, density.cutoff  , na.rm=T), min(allanswers$ec[allanswers$ec>0]))
			top		<- quantile(allanswers$ec, 1-density.cutoff, na.rm=T)
			limits[["ec.low.tcl"]]	<- tclVar(signif(bottom,3))
			limits[["ec.high.tcl"]]	<- tclVar(signif(top,3))
		}
		
		if (sum(!is.na(allanswers$hill))>2){
			bottom	<- quantile(allanswers$hill, 0.05  , na.rm=T)
			top		<- max(5, quantile(allanswers$hill, 0.9, na.rm=T))
			limits[["hill.low.tcl" ]]	<- tclVar(signif(bottom,3))
			limits[["hill.high.tcl"]]	<- tclVar(signif(top,3))
		}

		
		loglins <- list()					
		loglins$ec.islog.tkrb	<- tclVar(as.numeric(Options$plots$islog$ec))
		loglins$amp.islog.tkrb  <- tclVar(as.numeric(Options$plots$islog$amp))
		loglins$hill.islog.tkrb <- tclVar(as.numeric(Options$plots$islog$hill))
		
		# If z values were chosen as percent positive, set appropriate defaults
		# Used for plotting fits within bins upon click
		if (Options$analysis$zchoice == "pct") {
			limits[["Z.min"]] 			<- tclVar(0)
			limits[["Z.max"]]			<- tclVar(100)
			limits[["amp.low.tcl"]]		<- tclVar(0)
			limits[["amp.high.tcl"]]	<- tclVar(100)
			Options$plots$islog$amp		<- F
		}else{
			# Otherwise, set values using output from fits
			ztop		<- 10^quantile(allanswers$top, (1-density.cutoff), na.rm=T)
			limits[["Z.max"]]			<- tclVar(signif(ztop,3))
			zbot		<- 10^quantile(allanswers$bas, density.cutoff, na.rm=T)
			if (Vars$islog$z == T) {
				limits[["Z.min"]] 		<- tclVar(max(signif(zbot,3),1))
			}else{
				limits[["Z.min"]] 		<- tclVar(signif(zbot,3))
			}
			#show(sprintf("z limits are %.3g and %.3g", zbot, ztop))
		}
		
		if(Options$analysis$hillchoice=="free"){
			allstats	<- data.frame(ecRE=as.vector(sapply(ecRE,c)), ampRE=as.vector(sapply(ampRE,c)), hillRE=as.vector(sapply(hillRE,c)))
		}else{
			allstats	<- data.frame(ecRE=as.vector(sapply(ecRE,c)), ampRE=as.vector(sapply(ampRE,c)))
		}
		allstats	<- na.omit(allstats)
		allstats	<- allstats[allstats$ecRE>0,]
		allstats	<- allstats[allstats$ampRE>0,]
		allstats$ecRE	<- sort(allstats$ecRE)
		allstats$ampRE	<- sort(allstats$ampRE)
		if(Options$analysis$hillchoice=="free"){
			allstats	<- allstats[allstats$hillRE>0,]
			allstats$hillRE	<- sort(allstats$hillRE)
			Options$plots$hillREcutoff.idx	<- length(allstats$hillRE)
		}
		Options$plots$ecREcutoff.idx	<- length(allstats$ecRE)
		Options$plots$ampREcutoff.idx	<- length(allstats$ampRE)
		
		ALLSTATS	<<- allstats
		
		
	}
	
	############################################################################
	##  Variables to locate clicks for plotting of individual fits
	parPlotSize 	<- c()
	usrCoords 		<- c()
	
	# scaling variables for plots
	ysize 		<- 1
	xsize 		<- 1
	myhscale	<- 0.82
	# eventually, these will be set dynamically by window size
	#show("Finished setting initial variable values")
	
	############################################################################
	
	
	# Function to plot rectangles.  Essentially a custom "image" function.
	RectPlot <- function(graphtitle, flipjet=F, inname, Options, Bins, Answer.curr){
		#show("Plotting")
		inmatrix	<-Answer.curr[[inname]]
		inREmatrix	<- Answer.curr[[paste(inname,"RE", sep="")]]
		
		xlim	<- c(Options$plots$min$x, Options$plots$max$x)
		ylim	<- c(Options$plots$min$y, Options$plots$max$y)
		
		plot(xlim, ylim, pch=".", col="white", xlim=xlim, ylim=ylim,
				xlab=NA, ylab=NA, type="n", xaxt="n", yaxt="n"
				)
				
		title(paste(Vars$xyzLabels[[titr.current]]$z,graphtitle))
		
		niceaxes(xy="x", limits=xlim, axtitle=Vars$xyzLabels[[titr.current]]$x, islog=Vars$islog$x)
		niceaxes(xy="y", limits=ylim, axtitle=Vars$xyzLabels[[titr.current]]$y, islog=Vars$islog$y)
		
		# Determine scaling for bin color and colorbar
		minval	<- as.numeric(tclvalue(limits[[paste(inname,".low.tcl", sep="")]]))
		maxval	<- as.numeric(tclvalue(limits[[paste(inname,".high.tcl", sep="")]]))
		
		# Determine color for each bin (using Colorspec) and plot 
		for (binynum in 1:Bins$ny){
			for (binxnum in 1:Bins$nx){
				
				fitTF		<- F
				hillfitTF	<- F
				colr.rgb <- rgb(1,1,1)
				
				#show(c(binxnum, binynum))
				if(Options$analysis$method=="hillfit"){
					# This checks the stats to see if the point should be plotted
					method	<- Answer.curr$method[binxnum, binynum]
					hillfitTF	<- grepl("Hill",method)
					
					pval	<- Answer.curr$pvals[binxnum, binynum]
					reName	<- paste(inname,"RE", sep="")
					RE		<- Answer.curr[[reName]][binxnum, binynum]
					#show(sprintf("pvalue = %.2g, RE = %.2g", pval, RE))
					if (is.na(pval) | is.na(RE)) {
						fit.TF	<- F
						# if either pval or RE is NA, don't plot
					}else{
						pTF			<- pval < Options$plots$pcutoff
						reTF		<- RE < allstats[[reName]][Options$plots[[paste(reName,"cutoff.idx",sep="")]]]
						andOr		<- Options$plots$p.AndOr.RE #will either be "&" or "|"
						if (is.na(pTF)|is.na(reTF)){
							fitTF	<- F
						}else{
							fitTF	<- eval(parse(text=paste(pTF,andOr,reTF)))
						}
					}
					if(Options$analysis$zchoice=="pct"){
						fitTF	<- T
						# The statistics for deciding whether to plot values fit by % positive don't work
						# so we're just going to plot all of them.
					}
				}
				
				if ((fitTF & hillfitTF) | (Options$analysis$method == "subtract") ) {
					zval	<- inmatrix[binxnum,binynum]
					if(is.nan(zval) | is.na(zval)){ #ignore
					}else{
						if(zval != 0){
							neg.ok		<- (inname == "amp")
							#show(sprintf("Rectplot, Before Colorspec: minval=%.2g, maxval=%.2g, zval=%.2g", minval, maxval, zval))
							colr.rgb 	<- Colorspec(zval, minval, maxval, flipjet, logdisp=Options$plots$islog[[inname]], negvals.ok=neg.ok)
						
						}
					}
				}
				rect(Bins$seqx[binxnum],Bins$seqy[binynum],Bins$seqx[binxnum+1], Bins$seqy[binynum+1],col=colr.rgb, border=NA)
			}
		}
		parPlotSize <<- par("plt")
		usrCoords   <<- par("usr")
	}
	
	
	whichfun <- "PlotEC" # Default value
	if (num.doses <=3) { 
		whichfun = "PlotAmp" 
	}
	
	PlotGraph <- function(){
		# For resetting variables at the end of the plot
		#show("PlotGraph")
		oldpar <- par(no.readonly = T)
		par(bg="white")
		parplt1 <- par("plt")[1]
		parplt2	<- 0.99
		parplt3 <- par("plt")[3]
		parplt4 <- par("plt")[4]
		par("plt"=c(parplt1, parplt2, parplt3, parplt4))
		
		ysize 	<- parplt4-parplt3
		xsize 	<-	parplt2-parplt1
		myhscale <- (ysize/xsize)
		##show(myhscale)
		
		if (whichfun == "PlotEC"){
			if(Options$analysis$method=="subtract"){
				plot(c(0,4),c(0,4),pch=".", col="white", xlab=NA, ylab=NA, xaxt="n", yaxt="n")
				text(2,2, "No EC50s calculated \n (Too few doses in titrations.)")
			}else{
				RectPlot(paste(Vars$prefix[[titr.current]],"EC50 values"), flipjet=T, inname="ec", Options, Bins, Answers[[titr.current]])
			}
		}
		if (whichfun == "PlotAmp"){
			RectPlot(paste(Vars$prefix[[titr.current]],"Amplitudes") , flipjet=F, inname="amp",Options,Bins, Answers[[titr.current]])
		}
		if(whichfun == "PlotHill") {
			if (Options$analysis$hillchoice == "fixed"){
				plot(c(0,4),c(0,4),pch=".", col="white", xlab=NA, ylab=NA, xaxt="n", yaxt="n")
				text(2,2, "Fits performed with \n hill coefficients fixed \n at 1")
			}else {
				RectPlot(paste(Vars$prefix[[titr.current]],"Hill Coeffs"), flipjet=F, inname="hill", Options, Bins, Answers[[titr.current]])
			}
		}
		par(oldpar)
	}
	
	PlotColorbar <- function(){
		#show("PlotColorbar")
		if (whichfun == "PlotEC"){
			if (Options$analysis$method == "subtract"){
				plot(c(0,1,0,1), col="white", xlab=NA, ylab=NA, type="n", xaxt="n", yaxt="n")
			}else{
				bottom 	<- as.numeric(tclvalue(limits$ec.low.tcl))
				top		<- as.numeric(tclvalue(limits$ec.high.tcl))
				Colorbar(bottom, top, Vars$units[[titr.current]], islog=Options$plots$islog$ec, flip=T)
			}
		}
		if (whichfun == "PlotAmp"){
			bottom 	<- as.numeric(tclvalue(limits$amp.low.tcl))
			top		<- as.numeric(tclvalue(limits$amp.high.tcl))
			if (Options$analysis$zchoice =="pct") {
				Colorbar(bottom, top, "% positive", islog=Options$plots$islog$amp)
			}
			else {
				Colorbar(bottom, top, "A.U.", islog=Options$plots$islog$amp)
			}
		}
		if(whichfun == "PlotHill") {
			if (Options$analysis$hillchoice == "fixed"){
				plot(c(0,1,0,1), col="white", xlab=NA, ylab=NA, type="n", xaxt="n", yaxt="n")
			}
			else {
				bottom 	<- as.numeric(tclvalue(limits$hill.low.tcl))
				top		<- as.numeric(tclvalue(limits$hill.high.tcl))
				Colorbar(bottom, top, "", islog=Options$plots$islog$hill)
			}
		}
	}
	
	# Functions for buttons
	ECPlot	<- function(){
		whichfun <<- "PlotEC"
		tkrreplot(plot3d.img)
		tkrreplot(colorbar.img)
	}
	AmpPlot	<- function(){
		whichfun <<- "PlotAmp"
		tkrreplot(plot3d.img)
		tkrreplot(colorbar.img)
	}
	HillPlot	<- function(){
		whichfun <<- "PlotHill"
		tkrreplot(plot3d.img)
		tkrreplot(colorbar.img)
	}
	
	NextTitr	<- function(){
		if(titr.current < num.titr){
				titr.current <<- titr.current+1
			}else{
				titr.current <<- 1
			}
		tkrreplot(plot3d.img)
		tkrreplot(colorbar.img)
	}
	
	PrevTitr	<- function(){
		if(titr.current > 1){
				titr.current <<- titr.current-1
			}else{
				titr.current <<- num.titr
			}
		tkrreplot(plot3d.img)
		tkrreplot(colorbar.img)
	}
	
	####
	OnBack <- function(){
		tkdestroy(t5)
		BatchPrep(Vars, Track, Options,Bins)
	}
	####
	
	
	
	
	
	#############################################################################
	##### OPTIONS ###############################################################
	pop.rbval	<- tclVar("means")
	Set3dOptions.f <- function(){
		
		entries 	<- list()
		loglins		<- list()
		showest		<- list()
		
		t5.options 	<- tktoplevel()
		tkwm.title(t5.options, " 3D Analysis Options ")
		tkgrid(tklabel(t5.options, text="Set Options for main plots"))
		tkgrid(ttkseparator(t5.options, orient="horizontal"), sticky="ew")
		
		
		t5.options.f1 <- tkframe(t5.options)
		tkgrid(tklabel(t5.options.f1, text="   "), tklabel(t5.options.f1, text="Min"),tklabel(t5.options.f1, text="Max"))
		
		## Min/Max of displayed values #########################################
		##show(Vars$key)
		for (i in 1:length(Vars$key$short)){
			idx <- Vars$key$short[i]
			#show(tclvalue(limits[[paste(idx,".low.tcl", sep="")]]))
			entries[[paste(idx,".low.ent",sep="")]]		<- tkentry(t5.options.f1, width="11", textvariable=limits[[paste(idx,".low.tcl", sep="")]])
			entries[[paste(idx,".high.ent",sep="")]]	<- tkentry(t5.options.f1, width="11", textvariable=limits[[paste(idx,".high.tcl", sep="")]])
			rowlabel <- tklabel(t5.options.f1, text=Vars$key$long[i])
			tkgrid(rowlabel, entries[[paste(idx,".low.ent",sep="")]], entries[[paste(idx,".high.ent",sep="")]])
			tkgrid.configure(rowlabel, sticky ="e")
			
		}
		
		tkgrid(tklabel(t5.options.f1, text="  "))
		tkgrid(ttkseparator(t5.options.f1, orient="horizontal"), sticky="ew", columnspan=3)
		tkgrid(tklabel(t5.options.f1, text=" Scale: "), tklabel(t5.options.f1, text="Log"),tklabel(t5.options.f1, text="Linear"))
		
		## Log/Lin #############################################################
		for (i in 1:length(Vars$key$short)){
			idx <- Vars$key$short[i]
			loglins[[paste(idx,".islog.tkrb", sep="")]] <- tclVar(as.numeric(Options$plots$islog[[idx]]))
			loglins[[paste(idx,".log.rb",sep="")]]	<- tkradiobutton(t5.options.f1, variable=loglins[[paste(idx,".islog.tkrb", sep="")]], value=1)
			loglins[[paste(idx,".lin.rb",sep="")]]	<- tkradiobutton(t5.options.f1, variable=loglins[[paste(idx,".islog.tkrb", sep="")]], value=0)
			rowlabel <- tklabel(t5.options.f1, text=Vars$key$long[i])
			tkgrid(rowlabel, loglins[[paste(idx,".log.rb",sep="")]], loglins[[paste(idx,".lin.rb",sep="")]])
			tkgrid.configure(rowlabel, sticky ="e")
			
		}
		tkgrid(t5.options.f1)
		tkgrid(ttkseparator(t5.options, orient="horizontal"), sticky="ew")
		###################################################################################################################################################################
		##########################################################################################
		## Data cutoffs ########################################################
		p.options.fr 	<- tkframe(t5.options)
		pcutoff.tk 		<- tclVar(Options$plots$pcutoff)
		pcutoff.entry	<- tkentry(p.options.fr, width="6", textvariable=pcutoff.tk)
		pcutoff.label	<- tklabel(p.options.fr, text=" p-value cutoff for plots " )
		tkgrid(pcutoff.label, pcutoff.entry)
		tkgrid(p.options.fr)
		
		andor.fr	<- tkframe(t5.options)
		p.AndOr.RE	<- tclVar(Options$plots$p.AndOr.RE)
		pAndRE.rb	<- tkradiobutton(andor.fr, variable=p.AndOr.RE, value="&")
		pOrRE.rb	<- tkradiobutton(andor.fr, variable=p.AndOr.RE, value="|")
		pAndRE.label	<- tklabel(andor.fr, text="AND")
		pOrRE.label		<- tklabel(andor.fr, text="OR")
		tkgrid(pAndRE.label, pAndRE.rb, pOrRE.label, pOrRE.rb)
		tkgrid.configure(pAndRE.label, sticky="e")
		tkgrid.configure(pOrRE.label, sticky="e")
		tkgrid.configure(pAndRE.rb, sticky="w")
		tkgrid.configure(pOrRE.rb, sticky="w")
		tkgrid(andor.fr)

		## Relative error sliders ###########
		sliders.fr	<- tkframe(t5.options)
		
		#ec50 ###
		ecREcutoff.idx		<- Options$plots$ecREcutoff.idx
		ecREs	<- allstats$ecRE
		ecREs	<- sort(ecREs)
		ecREslider.idx	<- tclVar(ecREcutoff.idx)
		ecREslider.val	<- tclVar(as.character(ecREs[as.numeric(tclvalue(ecREslider.idx))]))
		ecSliderVal.label	<- tklabel(sliders.fr, text=ecREslider.val)
		
		ecUpdateLabel	<- function(...){
			ecREslider.val <<- tclVar(as.character(signif(ecREs[as.numeric(tclvalue(ecREslider.idx))],3)))
			tkconfigure(ecSliderVal.label, textvariable=ecREslider.val)
			tkrreplot(plot3d.img)
			Options$plots$ecREcutoff.idx	<<- as.numeric(tclvalue(ecREslider.idx))
			Options$plots$ampREcutoff.idx	<<- as.numeric(tclvalue(ampREslider.idx))
		}
		
		tkgrid(tklabel(sliders.fr, text="EC50 RE cutoff: "),ecSliderVal.label)
		tkconfigure(ecSliderVal.label, textvariable=ecREslider.val)
		ecSlider.sc	<- tkscale(sliders.fr, from=1, to=length(allstats$ecRE), showvalue=F,
					variable=ecREslider.idx, resolution=1, orient="horizontal",
					command=ecUpdateLabel)
		tkgrid(ecSlider.sc, columnspan=2)
		
		#amplitude ###
		ampREcutoff.idx		<- Options$plots$ampREcutoff.idx
		ampREs	<- allstats$ampRE
		ampREs	<- sort(ampREs)
		ampREslider.idx	<- tclVar(ampREcutoff.idx)
		ampREslider.val	<- tclVar(as.character(ampREs[as.numeric(tclvalue(ampREslider.idx))]))
		ampSliderVal.label	<- tklabel(sliders.fr, text=ampREslider.val)
		ampUpdateLabel	<- function(...){
			ampREslider.val <<- tclVar(as.character(signif(ampREs[as.numeric(tclvalue(ampREslider.idx))],3)))
			tkconfigure(ampSliderVal.label, textvariable=ampREslider.val)
			tkrreplot(plot3d.img)
			Options$plots$ecREcutoff.idx	<<- as.numeric(tclvalue(ecREslider.idx))
			Options$plots$ampREcutoff.idx	<<- as.numeric(tclvalue(ampREslider.idx))
		}
		tkgrid(tklabel(sliders.fr, text="Amp RE cutoff: "),ampSliderVal.label)
		tkconfigure(ampSliderVal.label, textvariable=ampREslider.val)
		ampSlider.sc	<- tkscale(sliders.fr, from=1, to=length(allstats$ampRE), showvalue=F,
					variable=ampREslider.idx, resolution=1, orient="horizontal",
					command=ampUpdateLabel)
		tkgrid(ampSlider.sc, columnspan=2)
		
		tkgrid(sliders.fr)
		tkgrid(ttkseparator(t5.options, orient="horizontal"), sticky="ew")
		
		
		
		# Buttons ##########################
		ApplyChg	<- function(){
			
			Options$plots$islog$ec		<<- as.logical(as.numeric(tclvalue(loglins$ec.islog.tkrb)))
			Options$plots$islog$amp		<<- as.logical(as.numeric(tclvalue(loglins$amp.islog.tkrb)))
			Options$plots$islog$hill	<<- as.logical(as.numeric(tclvalue(loglins$hill.islog.tkrb)))
			Options$plots$pcutoff		<<- as.numeric(tclvalue(pcutoff.tk))
			Options$plots$ecREcutoff.idx	<<- as.numeric(tclvalue(ecREslider.idx))
			Options$plots$ampREcutoff.idx	<<- as.numeric(tclvalue(ampREslider.idx))
			Options$plots$p.AndOr.RE		<<- tclvalue(p.AndOr.RE)
			
			tkrreplot(plot3d.img)
			tkrreplot(colorbar.img)
			ScsOptions <<- Options
		}
		OnCloseOpt	<- function(){
			tkdestroy(t5.options)
		}
		
		tkbind(t5.options, "<Return>", ApplyChg)
		t5.options.b.fr <- tkframe(t5.options)
		Apply.but		<- tkbutton(t5.options.b.fr, text="  Apply ", command=ApplyChg)
		Close.t5.o.but	<- tkbutton(t5.options.b.fr, text="  Close ", command=OnCloseOpt)
		tkgrid(Close.t5.o.but, Apply.but)
		tkgrid(t5.options.b.fr)
		
		
	}
	
	############################################################################
	##### PLOTTING BINS ########################################################
	
	#Function to plot individual fits using left mouseclick as input
	OnLeftClick3D <- function(x,y){	
		
		#######################
		# Will find the point in a list closest to the user's click.
		# Mostly copied directly from example by Wettenhall
		
		xClick <- x
		yClick <- y
		width  <- as.numeric(tclvalue(tkwinfo("reqwidth",plot3d.img)))
		height <- as.numeric(tclvalue(tkwinfo("reqheight",plot3d.img)))
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
			points.idx	<- paste(x.idx, ",", y.idx, sep="")
			x.idx <- x.idx
			y.idx <- y.idx
		}
		else { points.idx <- "outside binned area" }
		#############################
		
		
		
		## Plotting
		PlotFit	<- function(){
			par(bg="white")
			plotdata	<- Answers[[titr.current]]$fitdata[[points.idx]]
			pdata		<- plotdata$indata
			PDATA		<<- pdata
			if(length(pdata)>0){
				
				
				# Set curve to NA, in case fitting failed
				fitconcs	<- NA
				fit			<- NA
				
				axislog 	<- "x"
				#if (Vars$islog$z && Options$analysis$zchoice =="mfi"){ 
				#	axislog <- "xy"
				#}
				allconcs	<- Vars$Concs[[titr.current]]
				zeropoint	<- min(allconcs[allconcs != 0])/100 #this is where 0 will be shown in log scale
				xlimits 	<- c(zeropoint,max(allconcs))
				
				concs 	<- pdata$concs
				zero.idx				<- which(concs==0)
				concs[zero.idx]			<- zeropoint
						
				logrange	<- log10( c(min(allconcs[allconcs != 0]), max(allconcs)))
				conc.pts	<- c(zeropoint, 10^pretty(logrange, n=4))
				toprange	<- pretty(logrange, n=4)
				show(conc.pts)
				conc.labels	<- expression()
				conc.labels[1]	<- "0"
				for (i in 1:length(toprange)){
					conc.labels[i+1]	<- as.expression(bquote(.(10)^.(toprange[i])))
				}
				show(conc.labels)
				
				cat("\n")			
				show(sprintf("Point index is %s", points.idx))
				show(pdata)
				
				if(points.idx != "outside binned area"){
					cat(sprintf('Data shown: ScsAnswers[[%i]]$fitdata[["%s"]]$indata\n', titr.current, points.idx))
					show(Answers[[titr.current]]$fitdata[[points.idx]]$indata)
		    	
				}
				# Determine y limits based on error bars
				errbar.tops	<- log10(pdata$zvalues)+pdata$sigmas
				errbar.bots	<- log10(pdata$zvalues)-pdata$sigmas
				ymin	<- min(errbar.bots)
				ymax	<- max(errbar.tops)
				
				# Calculate fit curve if fitting was successful
				#show("Calculating fit curve...")
				if(grepl("Hill",Answers[[titr.current]]$method[x.idx, y.idx])){
					
					bas		<- Answers[[titr.current]]$bas[x.idx, y.idx]
					amp		<- Answers[[titr.current]]$amp[x.idx,y.idx]
					ec		<- Answers[[titr.current]]$ec[x.idx,y.idx]
					hill	<- 1
					if (Options$analysis$hillchoice == "free") {
						if (!is.na(Answers[[titr.current]]$hill[x.idx, y.idx])){
							hill <- Answers[[titr.current]]$hill[x.idx, y.idx]
						}
					}
					fitconcs	<- 10^seq(log10(xlimits[1]), log10(xlimits[2]), length.out=100)
					
					if(Options$analysis$zchoice=="pct"){
						## Alternate values to be used if data was fit to percent positive
						yvals	<- pdata$zvalues
						ymax	<- min(max(yvals)+10, 100)
						ymin	<- max(min(yvals)-10, 0)
						fit <- bas+amp*fitconcs^hill/(fitconcs^hill + ec^hill)
						subtitle	<- "curve fit to percent positive"
						ylabel	<- paste("%", Vars$xyzLabels[[titr.current]]$z, "+")
						
					}else{
						
						fit			<- log10(bas+amp*fitconcs^hill/(fitconcs^hill + ec^hill))
						FITCONCS	<<- fitconcs
						FIT			<<- fit
						yvals	<- log10(pdata$zvalues)
						# Adjust y limits to include fits if necessary
						ymin	<- min(c(errbar.bots,fit))
						ymax	<- max(c(errbar.tops,fit))
						subtitle	<- sprintf("Fit p=%.2f, EC50 RE=%.3f, amp RE= %.3f", 
											Answers[[titr.current]]$pvals[x.idx, y.idx],
											Answers[[titr.current]]$ecRE[x.idx, y.idx], 
											Answers[[titr.current]]$ampRE[x.idx, y.idx]
											)
						ylabel	<- paste(Vars$xyzLabels[[titr.current]]$z, "MFI")
					}
					
					show("Fitting output for this index:")
					show(sprintf("EC50 = %.2g; Amplitude = %.2g; Hillslope = %.2g", ec, amp, hill))
				}else{
					show("Fitting unsuccessful")
				}
				
				
				ylimits	<- c(ymin,ymax)
				graphtitle 	<- sprintf("Hill fit curve for bin (%i, %i)", x.idx, y.idx)
				
				## Plotting points ± error bars
				if(Options$analysis$zchoice=="pct"){
					# Plot points only if fit to percent positive
					plot(concs, yvals, log=axislog,
						type="p", pch=20, col="black", cex=1.2,
						xlim=xlimits, ylim=ylimits, xaxt="n", ylab=NA, xlab=NA
					)
					
					
				}else{
					# Plot points with error bars for fit to means
					plotCI(concs, yvals, (pdata$sigmas)/sqrt(pdata$N), log=axislog,
						type="p", pch=20, col="black", gap=0,
						xlim=xlimits, ylim=ylimits, xaxt="n", ann=F, ylab=NA, xlab=NA
					)
					
				}
				
				
				## Add information about fit
				if(grepl("Hill",Answers[[titr.current]]$method[x.idx, y.idx])){
					lines(fitconcs, fit, col="red", lwd=2)
				}else{
					subtitle <- "Hill fitting unsuccessful"
				}
				
				
				## Add titles
				title(xlab=sprintf("concentration (%s)", Vars$units[[titr.current]]), ylab=ylabel, cex.axis=0.7, line=2.5,
					main=graphtitle, cex.main=0.9)	
				title(sub=subtitle, cex.sub=0.75, line=3.5)
				axis(side=1, at=conc.pts, labels=conc.labels)
			
			}
		}
		
		
		
		## Options for fit plots ##
		OnOptions	<- function(){
			t6.options	<- tktoplevel()
			means.rb	<- tkradiobutton(t6.options)
			tkconfigure(means.rb, variable=pop.rbval, value="means")
			alldata.rb	<- tkradiobutton(t6.options)
			tkconfigure(alldata.rb, variable=pop.rbval, value="alldata")
			
			tkgrid(tklabel(t6.options, text="   Set options for pop-out plots   "), columnspan=3)
			tkgrid(tklabel(t6.options, text=" Means only "), means.rb)
			tkgrid(tklabel(t6.options, text=" All Data "), alldata.rb)
			
			tkgrid(tklabel(t6.options, text="   "))
			zmin.label	<- tklabel(t6.options, text=paste(" min",Vars$xyzLabels[[titr.current]]$z))
			zmax.label	<- tklabel(t6.options, text=paste(" max",Vars$xyzLabels[[titr.current]]$z))
			zmin.entry	<- tkentry(t6.options, width="11", textvariable=limits[["Z.min"]])
			zmax.entry	<- tkentry(t6.options, width="11", textvariable=limits[["Z.max"]])
			tkgrid(zmin.label, zmin.entry)
			tkgrid(zmax.label, zmax.entry)
			tkgrid.configure(zmin.label, zmax.label, sticky="e")
			
			tkgrid(tklabel(t6.options, text="   "))
			ApplyChg	<- function(){
				tkrreplot(fit.img)
				Options$plots$maxz	<- tclvalue(limits$Z.max)
				Options$plots$minz	<- tclvalue(limits$Z.min)
			}
			OnCloseOpt	<- function(){
				tkdestroy(t6.options)
			}
			
			Apply.but		<- tkbutton(t6.options, text="  Apply ", command=ApplyChg)
			Close.t6.o.but	<- tkbutton(t6.options, text="  Close ", command=OnCloseOpt)
			tkgrid(tklabel(t6.options, text="   "))
			tkgrid(Close.t6.o.but, Apply.but)
		}
		
		
		## Layout ##
		if (points.idx != "outside binned area"){
			t6 		<- tktoplevel()
			tkwm.title(t6, "Bin Fitting")
			OnClose		<- function(){tkdestroy(t6)}
			fit.img		<- tkrplot(t6, fun=PlotFit, hscale=0.8, vscale=0.8)
			close.but	<- tkbutton(t6, text="  Close  ", command=OnClose)
			options.but	<- tkbutton(t6, text=" Options ", command=OnOptions)
			tkgrid(fit.img, columnspan=3)
			tkgrid(tklabel(t6, text=" "))
			tkgrid(close.but, options.but)
		}
	}
	                                     
	##############################################################################################
	####### SAVING DATA #########################    ###############################################
	
	
	OnSave	<- function(){
		t5.save	<- tktoplevel()
		tkwm.title(t5.save, "Save")
		tkgrid(tklabel(t5.save, text=" File Name "), columnspan=2)
		
		filename.tk	<- tclVar("")
		filename.enter	<- tkentry(t5.save, width=40, textvariable=filename.tk)
		tkgrid(filename.enter, columnspan=2)
		
		
		
		onOK	<- function(){
			saveData	<- list(Vars=Vars, Answers=Answers, Bins=Bins, Track=Track, Options=Options)
			savefolder	<- tkchooseDirectory(initialdir=Filefolder)
			fullpath	<- paste(tclvalue(savefolder),"/", tclvalue(filename.tk), ".scs", sep="")
			tkdestroy(t5.save)
			save(saveData,file=fullpath, compress=T)
			
		}
		
		onCancel	<- function(){
			tkdestroy(t5.save)
		}
		
		ok.but		<- tkbutton(t5.save, text="   OK   ", command=onOK)
		cancel.but	<- tkbutton(t5.save, text=" Cancel ", command=onCancel)
		tkgrid(cancel.but, ok.but)
		
	}
	
	
	##############################################################################################
	######## EXPORT ##############################################################################
	OnExportOpt <- function() {
		t5.export <- tktoplevel()
		tkwm.title(t5.export, "Export Options")
		
		## Choosing options ##
		if(Options$analysis$hillchoice=="free"){
			exports 	<- c("ec","amp","hill")
		}else{
			exports	<- c("ec","amp")
		}
		
		expt.labels <- list()
		expt.vars.tables	<- list()
		expt.buts.tables	<- list()
		expt.vars.figs		<- list()
		expt.buts.figs		<- list()
		
		t5.export.top <- tkframe(t5.export)
		# Labels
		expt.labels[["ec"]]		<- tklabel(t5.export.top, text=" EC50 values ")
		expt.labels[["amp"]]	<- tklabel(t5.export.top, text=" Amplitudes ")
		expt.labels[["hill"]]	<- tklabel(t5.export.top, text=" Hill Coefficients ")
		tables.heading			<- tklabel(t5.export.top, text="Table")
		figs.heading			<- tklabel(t5.export.top, text="Figure")
		spacer1					<- tklabel(t5.export.top, text="  ")
		
		# Checkbox configuration and layout # Default is to export all
		tkgrid(spacer1, tables.heading, figs.heading)
		for (export in exports){ 
			expt.vars.tables[[export]]<- tclVar("1")
			expt.buts.tables[[export]]<- tkcheckbutton(t5.export.top)
			tkconfigure(expt.buts.tables[[export]], variable=expt.vars.tables[[export]])
			expt.vars.figs[[export]] <- tclVar("1")
			expt.buts.figs[[export]]<- tkcheckbutton(t5.export.top)
			tkconfigure(expt.buts.figs[[export]], variable=expt.vars.figs[[export]])
			tkgrid(expt.labels[[export]], expt.buts.tables[[export]], expt.buts.figs[[export]])
			tkgrid.configure(expt.labels[[export]], sticky="w")
		}
		if(Options$analysis$hillchoice=="fixed"){
			expt.vars.tables$hill	<- tclVar("0")
			expt.vars.figs$hill		<- tclVar("0")
		}			
		tkgrid(t5.export.top, columnspan=2)
		sep <- ttkseparator(t5.export, orient="horizontal")
		tkgrid(sep, sticky="ew", columnspan=2)
		
		
		## Exporting ###
		
		#Make bin names
		Xnames <- c()
		Ynames <- c()
		for (xbin.idx in 1:Bins$nx){
			Xnames[xbin.idx]<- paste("X", xbin.idx, sep="")
		}
		for (ybin.idx in 1:Bins$ny){
			Ynames[ybin.idx]<- paste("Y", ybin.idx, sep="")
		}
		
		# Export functions
		OnExport <- function(){
			#show(sapply(expt.vars.figs,tclvalue))
			#show(sapply(expt.vars.tables,tclvalue))
			exportfolder <- tkchooseDirectory(initialdir=Filefolder)
			show(exportfolder)
			exportfolder <- tclvalue(exportfolder)
			if(!nchar(exportfolder)){
				tkmessageBox(message="No folder selected")
			}else{
				exportheader <- date()
				exportheader <- paste(exportheader,"Values exported from ScatterSlice", sep="\n")
				exportheader <- paste(exportheader,"Rows are X bins, columns are Y bins", sep="\n")
				if(Options$analysis$zchoice=="pct"){
					zline	<- paste("Z values were calculated as percent with ",Vars$xyzLabels[[titr.current]]$z," higher than ", 100-Options$analysis$pctgate,"% of control.", sep="")
				}else{
					zline	<- paste("Z values calculated from ", Vars$xyzLabels[[titr.current]]$z, "MFI.", sep="")
				}
				exportheader <- paste(exportheader, zline, sep="\n")
				exportheader <- paste(exportheader, paste("Fit with", Options$hillchoice, "hill coefficient."), sep="\n")
				
				titr.now	<- titr.current
				export.time	<- system.time(
				foreach(titr.idx=1:num.titr) %do% {
					titr.current <<- titr.idx
					for (export in exports){
						filename.base	<- paste(exportfolder,"/", Vars$prefix[[titr.idx]], " " , paste(Vars$xyzLabels[[titr.idx]]$x, Vars$xyzLabels[[titr.idx]]$y, sep=","), sep="")
						if (tclvalue(expt.vars.tables[[export]]) == "1"){
							filename <- paste(filename.base, "_", export,".txt", sep="")
							write(export, filename)
							write(exportheader, filename, append=T)
							write("\n", filename, append=T)
							write.table(Bins$centersx, filename, append=T, row.names=F, col.names="X_Bin_Centers")
							write("\n", filename, append=T)
							write.table(Bins$centersy, filename, append=T, row.names=F, col.names="Y_Bin_Centers")
							write("\n", filename, append=T)
							write.table(Answers[[titr.idx]][[export]], filename, sep="\t", append=T, row.names=Xnames, col.names=Ynames)
						}
						#show(tclvalue(expt.vars.figs[[export]]))
						if (tclvalue(expt.vars.figs[[export]]) == "1") {
							#show("Outputting figure")
							filename <- paste(filename.base, "_", export,".pdf", sep="")
							pdf(filename)
							if (export == "ec") {
								whichfun <<- "PlotEC"
							}
							if (export == "amp") {
								whichfun <<- "PlotAmp"
							}
							if (export == "hill") {
								whichfun <<- "PlotHill"
							}
							
							PlotGraph()
							dev.off()
							filename <- paste(filename.base, "_", export,"_Colorbar.pdf", sep="")
							pdf(filename, width=1.75)
							PlotColorbar()
							dev.off()
						}
					}
				}
				)[3]
				titr.current	<- titr.now
			}
			#show(sprintf("export time was %.3f",export.time))
			tkdestroy(t5.export)
		}
		tkbind(t5.export, "<Return>", OnExport)
		final.export.but	<- tkbutton(t5.export, text="  Export...", command=OnExport)
		OnCancel.export <- function(){tkdestroy(t5.export)}
		export.cancel.but	<- tkbutton(t5.export, text="  Cancel   ", command=OnCancel.export)
		tkgrid(export.cancel.but, final.export.but, tklabel(t5.export, text="  "))
		tkgrid.columnconfigure(t5.export, 1, weight=1)
	}
	
	
	#### LAYOUT ####
	t5 <- tktoplevel()
	T5layout <- function(){
		tkwm.title(t5, "3D Analysis")
		graphs.frm	<- tkframe(t5)
		colorbar.img 	<<- tkrplot(graphs.frm, fun=PlotColorbar, hscale=0.25, vscale=1)
		plot3d.img 		<<- tkrplot(graphs.frm, fun=PlotGraph, hscale=myhscale)
		plotEC.but		<- tkbutton(t5, text="   EC50 Values   ", command=ECPlot)
		plotAmp.but		<- tkbutton(t5, text="    Amplitudes   ", command=AmpPlot)
		plotHill.but	<- tkbutton(t5, text="   Hill Coeffs   ", command=HillPlot)
		
		
		next.titr.but	<- tkbutton(t5, text=" >> ", command=NextTitr)
		prev.titr.but	<- tkbutton(t5, text=" << ", command=PrevTitr)
		tkgrid(plot3d.img, colorbar.img)
		tkgrid(graphs.frm, columnspan=3)
		titr.nav.frm	<- tkframe(t5)
		titr.label		<- tklabel(t5, text="Titrations")
		tkgrid(prev.titr.but, titr.label, next.titr.but)
		tkgrid(tklabel(t5, text="  "), tklabel(t5, text="Click a bin to display fit"))
		tkgrid(tklabel(t5, text="  "))
		tkgrid(plotEC.but, plotAmp.but, plotHill.but)
		tkgrid.configure(plotEC.but, sticky="ew")
		tkgrid.configure(plotAmp.but, sticky="ew")
		tkgrid.configure(plotHill.but, sticky="ew")
		sep <- ttkseparator(t5, orient="horizontal")
		tkgrid(tklabel(t5, text="  "))
		tkgrid(sep, columnspan=3, sticky="ew")
		
		nav.fr		<- tkframe(t5)
		back.but 	<- tkbutton(nav.fr, text=" << Back    ", command=OnBack)
		export.but	<- tkbutton(nav.fr, text=" Export...  ", command=OnExportOpt)
		options.but	<- tkbutton(nav.fr, text=" Options... ", command=Set3dOptions.f)
		save.but	<- tkbutton(nav.fr, text="   Save...  ", command=OnSave)
		tkgrid(back.but, options.but, save.but, export.but)
		tkgrid(nav.fr, columnspan=3)
		
		tkbind(plot3d.img, "<Button-1>", OnLeftClick3D)
		tkconfigure(plot3d.img,cursor="hand2")
	}
	
	T5layout()
}

