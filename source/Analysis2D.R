#######  PLOTTING X-Z AND Y-Z 2D PLOTS  ########
Plots2d <- function(Vars, Options, Bins){
	num.doses 	<- Vars$num.doses
	num.titr	<- Vars$num.titr
	#Concs		<- Vars$Concs
	thresh		<- Options$thresh
	
	histopt.tk	<- tclVar(Options$plots$hist2d)
	for(titr.idx in 1:Vars$num.titr){
		Bins <- BinPartition2d(Bins, Vars, Options,titr.idx)
	}
	TEMPindata <<- list()
	ScsBins <<- Bins
	
	#Fit hill equations
	FitHill2d <- function(Bins, xy, Vars, Options, titr.idx){
		#show(sprintf("Analyzing %s bins", xy))
		Concs <- Vars$Concs[[titr.idx]]
		
		binlist	<- Bins[[xy]]
		nbins	<- Bins[[paste("n",xy, sep="")]]
		#show(binlist)

		#show(sprintf("Fitting %i %s bins.", nbins, xy))
		if (Vars$num.doses < 4){
			# method if fewer than four doses
			method	<- vec2d
			Amps 	<- Bins[[xy]][["means"]][[which.max(Concs)]] - Bins[[xy]][["means"]][[which.min(Concs)]]
			Amps	<- Amps*((Bins[[xy]][["sizes"]][[which.max(Concs)]] >= Options$thresh)*(Bins[[xy]][["sizes"]][[which.min(Concs)]] >= Options$thresh))
			method[!is.nan(Amps)]	<- "subtract"
			method[is.nan(Amps)]	<- "Fail"
			method[((Bins[[xy]][["sizes"]][[which.max(Concs)]] < Options$thresh)|(Bins[[xy]][["sizes"]][[which.min(Concs)]] < Options$thresh))]	<- "Fail"
			Amps[is.nan(Amps)]		<- 0
			xyzfits[[paste(xy,"zamp", sep="")]]		<- Amps
			xyzfits[[paste(xy,"method", sep="")]]	<- method
		}else{
			
			#### PARALLEL FITTING ########################################################################
			fitting2D <- foreach(bin.idx=1:nbins) %do% { #################################################
				
				#fitting2D <- list()
				#for (bin.idx in 1:nbins) {
				# General method for most fitting
				
				xyzfits	<- list()
				  
				
				#show(sprintf("Analyzing %s bin number %i", xy, bin.idx))
				answer_lookup	<- paste(xy,bin.idx, sep="")
				sizevector	<- rep(0,num.doses)
				meanvector	<- rep(0,num.doses)
				posvector	<- rep(0,num.doses)
				sigvector	<- rep(0,num.doses)
				
				for(z.idx in 1:num.doses){
					sizevector[z.idx]	<- binlist$sizes[[titr.idx]][[z.idx]][bin.idx]
					meanvector[z.idx]	<- binlist$means[[titr.idx]][[z.idx]][bin.idx]
					posvector[z.idx]	<- binlist$pospct[[titr.idx]][[z.idx]][bin.idx]
					sigvector[z.idx]	<- binlist$sigmas[[titr.idx]][[z.idx]][bin.idx]
				}
				
				#show(c("sizevector", sizevector))
				okbins		<- which(sizevector >= thresh)
				n.okbins	<- length(okbins)
				
				if (n.okbins > 3) {
					# Choose indata based on whether analysis is percent or MFI
					if (Options$analysis$zchoice =="pct"){
						#show("collecting percentages for indata")
						indata 	<- na.omit(data.frame(concs=Concs[okbins], zvalues=posvector[okbins]))
					}
					else {
						#show("collecting MFIs for indata")
						indata 	<- na.omit(data.frame(concs=Concs[okbins], zvalues=meanvector[okbins], sigmas=sigvector[okbins], N=sizevector[okbins]))
					}
					
					answer_lookup <- paste(xy,bin.idx,sep="")
					
					
					TEMPindata[[answer_lookup]] <- list()
					TEMPindata[[answer_lookup]][["indata"]] <<- indata
					TEMPindata[[answer_lookup]][["concs"]] 	<<- Concs[okbins]
					TEMPindata[[answer_lookup]][["means"]]	<<- meanvector[okbins]
					
					
					#show(sprintf("indata has %i values with %i unique concs", length(indata$zvalues), length(unique(indata$concs))))
					
					outdata <- HillFit(indata, Options)  #<<<<<<<<<<<<<<<<-------------------HILL FITTING HERE 
					#show(indata)
					#show(outdata)
					if (outdata$method == "fail"){
						#show("Hill fitting failed completely")
						xyzfits[["models"]]	<- NA
						xyzfits[["method"]]	<- "Fail"
						xyzfits[["ec"]]		<- NA
						xyzfits[["ecRE"]]	<- NA
						xyzfits[["amp"]]	<- NA
						xyzfits[["ampRE"]]	<- NA
						xyzfits[["pvals"]]	<- NA
						xyzfits[["hill"]]	<- NA
						xyzfits[["indata"]]	<- indata
					}else{ 
						out	<- outdata$hilleq
						xyzfits[["indata"]]	<- indata
						xyzfits[["models"]]		<- out
						xyzfits[["method"]]		<- outdata$method
						
						bas                 <-  out$estimate[1]
						ec	                <-	out$estimate[3]
						amp					<-	out$estimate[2]
						xyzfits[["ec"]]		<- ec	
						xyzfits[["amp"]]	<- amp	
						
						if(length(out$estimate)==4){
							hill <- out$estimate[4]
							xyzfits$hill <- hill
						}else{
							hill <- 1
						}
						
						#show("Calculating residuals")
						fitVals	<- bas + amp*indata$concs^hill/(indata$concs^hill + ec^hill)
						
						res	<- log10(indata$zvalues)-log10(fitVals)
						
						RES	<<-res
						SIGVECTOR <<-sigvector
						
						if (Options$analysis$zchoice == "mfi"){
							N			<- length(indata$zvalues)
							chi2		<- sum(res^2*N/indata$sigmas^2)
							deg.free	<- N-length(out$estimate)
							pval		<- pgamma(0.5*chi2, 0.5*deg.free, 1, lower.tail = TRUE, log.p = FALSE)
						}else{
							pval	<- 1
						}
                        xyzfits[["pvals"]]		<- pval
						
						#show(sprintf("p-value = %.4f",pval))
						
						RE	<- try(sqrt( diag( 2*out$minimum/(deg.free)*solve(out$hessian)) )/out$estimate)
						if(is.numeric(RE)){
							xyzfits[["ecRE"]]	<- RE[3]
							xyzfits[["ampRE"]]	<- RE[1]
						}
					}
					
				}else{ # if n.okbins <= 3
					#show(sprintf("%s slice number %i has too few full bins for hill fit ", xy, bin.idx))
					xyzfits[["amp"]]	<- NA
					#Fix -> Need method to estimate amp	
					xyzfits$models 	<- 0
					xyzfits[["method"]]	<- "None"
					xyzfits$models		<- 0
					xyzfits[["ec"]]		<- NA
					xyzfits[["ecRE"]]	<- NA
					xyzfits[["ampRE"]]	<- NA
					#xyzfits[["indata"]]	<- indata
				}
				return(xyzfits)
			}
			#show("Fitting bins complete")
		}
		#show(fitting2D)
		return(fitting2D)
	}
	fits2D	<- list()
	fits2D[["x"]]	<- list()
	fits2D[["y"]]	<- list()
	
	fittingtime	<- system.time(
	for(titr.idx in 1:num.titr){
		fits2D[["x"]][[titr.idx]] <- FitHill2d(Bins, "x", Vars, Options, titr.idx)
		fits2D[["y"]][[titr.idx]] <- FitHill2d(Bins, "y", Vars, Options, titr.idx)
	})[3]
	
	
	
	#show(sprintf("Fitting time was %.3f",fittingtime))
	
	TEMPfits2D <<- fits2D
	#show(fits2D)
	
	#show(fits2D)
	t5.2dplots 		<- tktoplevel()
	tkwm.title(t5.2dplots, "2D Analysis")
	plots2d.list 	<- list()
	
	## To be fixed when storage names are changed to match the rest of the program
	bincenters		<- list()
	columns			<- list()
	bincenters[["x"]] <- Bins$centersx
	bincenters[["y"]] <- Bins$centersy
	columns[["x"]]	<- Vars$xyzLabels[[titr.idx]][["x"]]

	
	## Collect data into nice vectors and lists
	collected	<- list()
	collected[["x"]]	<- list()
	collected$x[["ec"]] 	<- list()
	collected$x[["ecRE"]] 	<- list()
	collected$x[["amp"]] 	<- list()
	collected$x[["ampRE"]] 	<- list()
	if(Options$analysis$hillchoice=="free"){
		collected$x[["hill"]] 	<- list()
		collected$x[["hillRE"]]	<- list()
	}
	for(xy in c("x","y")){
		nbins	<- Bins[[paste("n",xy, sep="")]]
		for(output.idx in names(collected$x)){
			for(titr.idx in 1:num.titr){
				colname <- paste(Vars$prefix[[titr.idx]],Vars$xyzLabels[[titr.idx]][[xy]], sep=":")
				vec	<- c()
				for(bin.idx in 1:nbins){
					replacement	<- fits2D[[xy]][[titr.idx]][[bin.idx]][[output.idx]]
					vec[bin.idx] <- ifelse(length(replacement)>0,replacement,NA)
				}
				collected[[xy]][[output.idx]][[colname]] <- vec
	}}}
	
	COLLECTED <<- collected
	
	
	whichplot	<- "plotEC"
	
	ecPlot <- function(){
		whichplot	<<- "plotEC"
		tkrreplot(plot2ds.img)
	}
	
	ampPlot <- function(){
		whichplot	<<- "plotAmp"
		tkrreplot(plot2ds.img)
	}
	
	# Colors and pch to be used for plotting
	pointchars	<- c(15,16,17,18,19,20,1,2,5,6,7,9,10)
	pointchars	<- rep(pointchars, times=ceiling(num.titr/length(pointchars)))
	colorchoices<- c(33,128,448,90,642,29,575,143,120,564,497,142)
	colorchoices<- rep(colorchoices, times=ceiling(num.titr/length(colorchoices)))
	
	plot2ds.func <- function(){
		#show("Plotting")
		
		par(mfrow=c(1,3))	
		for (xy in c("x","y")){
			
			
			
			if(whichplot=="plotEC"){
				ylabel <- sprintf("EC50 (%s)",unique(Vars$units))
				output.idx	<- "ec"
			}
			else{
				ylabel <- "Amplitude"
				output.idx	<- "amp"
			}
			
			method.list	<- collected[[xy]][[output.idx]]$method
			RE.list	<- collected[[xy]][[paste(output.idx,"RE", sep="")]]
			
			if(num.doses < 4){
				# Need to fix this, there isn't really a good method for this at the moment 	
				zvalues <- zvalues[method=="subtract"]
				xvalues	<- bincenters[[xy]][method=="subtract"]
			}
			else{
				xvalues	<- bincenters[[xy]]			
			}
			
			#show("xvalues")
			#show(xvalues)
			
			zvalues	<- collected[[xy]][[output.idx]]
			#show("zvalues")
			#show(zvalues)
			totalz	<- unlist(zvalues)
			totalz	<- totalz[!is.na(totalz)]
			totalz	<- totalz[totalz>0]
			
			#show("totalz")
			#show(totalz)
			
			label.parts <- c()
			for(titr.idx in 1:length(Vars$xyzLabels)){
				label.parts[titr.idx]<- Vars$xyzLabels[[titr.idx]][[xy]]
			}
			xlabel <- paste(xy,"bins:", paste(unique(label.parts), sep=" / "))
			plot(range(xvalues[!is.na(xvalues)]), range(totalz), log="y", xlab=xlabel, ylab=ylabel, pch="", asp=1)

			pch.idx 	<- 1
			color.idx		<- 1
			for(output.idx in names(zvalues)){
				lines(xvalues,zvalues[[output.idx]], col=colors()[colorchoices[color.idx]], lwd=2)
				points(xvalues,zvalues[[output.idx]], col=colors()[colorchoices[color.idx]], pch=pointchars[pch.idx], cex=2)
				pch.idx 	<- pch.idx+1
				color.idx		<- color.idx+1
			}
			
			
		}
		
		#legend
		pch.idx 	<- 1
		color.idx	<- 1
		plot(c(0,1),c(0,1), , pch="", xlab=NA, ylab=NA, xaxt="n", yaxt="n", xlim=c(0,1), ylim=c(0,1))
		for(output.idx in 1:length(zvalues)){
			#show(output.idx)
			yval <- 1-(output.idx)/10
			#show(yval)
			points(.05,yval, lty=1, pch=pointchars[pch.idx], col=colors()[colorchoices[color.idx]], cex=1.5)
			lines(c(0,.1), c(yval,yval), col=colors()[colorchoices[color.idx]])
			text(0.15,yval,ScsVars$prefix[[output.idx]], pos=4, col=colors()[colorchoices[color.idx]], offset=0)
			pch.idx 	<- pch.idx+1
			color.idx	<- color.idx+1
		}
		
	}
	
	options2ds.func	<- function(){
		options2d.win	<- tktoplevel()
		
		histoptions.frm	<- tkframe(options2d.win)
		hist.label		<- tklabel(histoptions.frm, text="Histogram:  ")
		off.label		<- tklabel(histoptions.frm, text="Off")
		bw.label		<- tklabel(histoptions.frm, text="B/W")
		col.label		<- tklabel(histoptions.frm, text="Color")
		spacer			<- tklabel(histoptions.frm, text="  ")
		hist.off.rb		<- tkradiobutton(histoptions.frm, variable=histopt.tk, value="off")
		hist.bw.rb		<- tkradiobutton(histoptions.frm, variable=histopt.tk, value="bw")
		hist.col.rb		<- tkradiobutton(histoptions.frm, variable=histopt.tk, value="col")
		tkgrid(spacer, off.label, bw.label, col.label)
		tkgrid(hist.label, hist.off.rb, hist.bw.rb, hist.col.rb)
		
		tkgrid(histoptions.frm)
	}
	
	export2ds.func <- function(){
		exportfolder <- tkchooseDirectory(initialdir=Filefolder)
		exportfolder <- tclvalue(exportfolder)
		if(!nchar(exportfolder)){
		tkmessageBox(message="No folder selected")}
		else{
			filename <- paste(c(exportfolder,"/",Vars$prefix,"Analysis_in_2D.txt"), collapse="")
			#show(filename)
			exportheader <- "Values exported from ScatterSlice"
			exportheader <- paste(exportheader,date(), sep="\n")
			write(exportheader, filename, append=T)
			write("\n", filename, append=T)
			write("MFI values represent bin centers", filename, append=T)
			write("\n", filename, append=T)
			
			
			#Output 2D fits
			for(xy in c("x","y")){
				for(whichout in names(collected$x)){
			
					write(whichout, filename, append=T)
					outdata <- data.frame(bincenters[[xy]],collected[[xy]][[whichout]])
					write.table(outdata, filename, append=T, row.names=F, col.names=c(xy, names(collected[[xy]][[whichout]])), quote=F, sep="\t")
					write("\n", filename, append=T)
			}}
		}
	}
	
	
	plot2ds.img		<<- tkrplot(t5.2dplots, fun=plot2ds.func, hscale=1.5, vscale=0.75)
	tkgrid(plot2ds.img)
	
	ecOrAmp.fr	<- tkframe(t5.2dplots)
	ec.but		<- tkbutton(ecOrAmp.fr, text="   EC50    ", command=ecPlot)
	amp.but		<- tkbutton(ecOrAmp.fr, text=" Amplitude ", command=ampPlot)
	tkgrid(ec.but, tklabel(ecOrAmp.fr, text= "              "), amp.but)
	tkgrid(ecOrAmp.fr)
	
	sep <- ttkseparator(t5.2dplots, orient="horizontal")
	tkgrid(sep, columnspan=1, sticky="ew")
	sep <- ttkseparator(t5.2dplots, orient="horizontal")
	tkgrid(sep, columnspan=1, sticky="ew")
	
	nav.fr <- tkframe(t5.2dplots)
	Close2ds.func	<- function(){tkdestroy(t5.2dplots)}
	close2ds.but	<- tkbutton(nav.fr, text="   Close     ", command=Close2ds.func)
	options.but		<- tkbutton(nav.fr, text="   Options   ", command=options2ds.func)
	export2ds.but	<- tkbutton(nav.fr, text="   Export    ", command=export2ds.func)
	tkgrid(close2ds.but)
	tkgrid(export2ds.but)
	tkgrid(nav.fr)
}
