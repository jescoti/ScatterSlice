FinalAnalyze <- function(Vars, Options, Bins){
	show("FinalAnalyze")
	Answers <- list()
	for(titr.idx in 1:Vars$num.titr) {
		hillwin	<- 0
		hilltot	<- 0
		
		show(sprintf("@@@@@@@@@@@@@   TITRATION %i @@@@@@@@@@@@@@@",titr.idx))
		thresh 		<- Options$thresh
		models		<- list()
		drc.output	<- list()
		fits		<- list()
		
		outmatx		<- matrix(rep(NA,times=Bins$nx*Bins$ny), ncol=Bins$ny)
		EC50s		<- outmatx
		Amps		<- outmatx
		Hill		<- outmatx
		Method		<- outmatx
		Bas			<- outmatx
		pvals		<- outmatx
		EC.RE		<- outmatx
		Amp.RE		<- outmatx
		Hill.RE		<- outmatx
		
		FinalZData	<<- list()
		
		
		Concs 			<- Vars$Concs[[titr.idx]]
		num.doses		<- Vars$num.doses
		
		minpoints	<- 4
		if(Options$analysis$hillchoice=="free"){ minpoints <- 4 }
		
		
		#### BEGIN TIMING  #####
		starttime <- proc.time()[1]
		
		
		#If fewer than 4 files, just subtract min conc z value means from max conc
		# to calculate Amplitudes
		if(Options$analysis$method=="subtract") {
			show("Subtracting background for amplitude estimation")
			Amps 	<- Bins$means[[titr.idx]][[which.max(Concs)]] - Bins$means[[titr.idx]][[which.min(Concs)]]
			Amps	<- Amps*((Bins$sizes[[titr.idx]][[which.max(Concs)]] >= Options$thresh)*(Bins$sizes[[titr.idx]][[which.min(Concs)]] >= Options$thresh))
			Method[!is.nan(Amps)]	<- "subtract"
			Method[is.nan(Amps)]	<- "Fail"
			Method[((Bins$sizes[[titr.idx]][[which.max(Concs)]] < Options$thresh)|(Bins$sizes[[titr.idx]][[which.min(Concs)]] < Options$thresh))]	<- "Fail"
			Amps[is.nan(Amps)]		<- NA
			Amps[Amps==0]	<- NA
		
		}else{
			### Perform Hill fits on binned data ####
			show("Hill Fitting")
			for (x.idx in 1:Bins$nx){
				for(y.idx in 1:Bins$ny){
					sizevector	<- vector("numeric", num.doses)
					meanvector	<- vector("numeric", num.doses)
					posvector	<- vector("numeric", num.doses)
					sigvector	<- vector("numeric", num.doses)
					
					
					#cat("\n\n@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n")
					show(sprintf("x = %i, y = %i",x.idx, y.idx))
					answer.lookup 	<- paste(x.idx,",",y.idx,sep="")
					
					# Collect summary data
					for(z.idx in 1:num.doses){
						show(sprintf("Collecting data for dose %i", z.idx))
						sizevector[z.idx]<- Bins$sizes[[titr.idx]][[z.idx]][x.idx, y.idx]
						meanvector[z.idx]<- Bins$means[[titr.idx]][[z.idx]][x.idx, y.idx]
						posvector[z.idx] <- Bins$pospct[[titr.idx]][[z.idx]][x.idx, y.idx]
						sigvector[z.idx] <- Bins$sigmas[[titr.idx]][[z.idx]][x.idx, y.idx]
					}
					show("Collected")
					# Determine which bins are sufficiently full
					okbins		<- which(sizevector >= thresh)
					badbins		<- is.na(meanvector)
					n.okbins	<- length(okbins)
					
					# Decide whether or not to attempt hill fitting
					hillfit.tf	<- (n.okbins >= minpoints)
					if(Options$analysis$zchoice =="pct"){
							control.size.ok	<- Bins$sizes[[titr.idx]][[Vars$Control.File[[titr.idx]]]][x.idx, y.idx] > thresh
							notNAbins	<- !is.na(posvector)
							pcts.ok	<- (sum(notNAbins) > thresh)
							hillfit.tf <- (control.size.ok & pcts.ok & hillfit.tf)
							show(c(control.size.ok, pcts.ok, hillfit.tf))
					}
					
					show(posvector)
					show(badbins)
					
					if (hillfit.tf) {
						show(sprintf("%i OK bins", n.okbins))
						#mean.zdata	<- rep(0, n.okbins)
						#means.idx	<- 1
						#zvalues		<-c(0)
						#concs  		<-c(0)
						#length(zvalues)	<- 0
						#length(concs)	<- 0
						
						if (Options$fit =="means") {
							indata <- data.frame(concs=Concs[okbins], zvalues=meanvector[okbins], sigmas=sigvector[okbins], N=sizevector[okbins])
						}
						
						if (Options$analysis$zchoice =="pct"){
							indata 		<- na.omit(data.frame(concs=Concs[okbins], zvalues=posvector[okbins], sigmas=sigvector[okbins], N=sizevector[okbins]))
						}
											
						FinalZData[[answer.lookup]][["indata"]] <<- indata
						FinalZData[[answer.lookup]][["concs"]] 	<<- Concs[okbins]
						FinalZData[[answer.lookup]][["means"]]	<<- meanvector[okbins]
						
						##### HILL FITTING #########################
						show(indata)
						show("Attempting Hill fit")
						outdata	<- HillFit(indata, Options)
						hilltot	<- hilltot+1
						OUTDATA	<<- outdata
						show(sprintf("hill fit was a %s",outdata$method))
						
						if (outdata$method != "fail"){
							show("Hill fit successful")
							hillwin	<- hillwin+1
							out	<- outdata$hilleq
							
							models[[answer.lookup]]		<- out
							
							Method[x.idx,y.idx]			<- outdata$method
							
							bas                 <-  out$estimate[1]
							ec	                <-	out$estimate[3]
							amp					<-	out$estimate[2]	
							
							Bas[x.idx,y.idx]	<-  out$estimate[1]
							EC50s[x.idx,y.idx]	<-	out$estimate[3]
							Amps[x.idx,y.idx]	<-	out$estimate[2] 
							
							if(length(out$estimate)==4){
								hill <- out$estimate[4]
								Hill[x.idx, y.idx] <- hill
							}else{
								hill <- 1
							}
							
							show("Calculating residuals")
							fitVals	<- bas + amp*indata$concs^hill/(indata$concs^hill + ec^hill)
							
							show(sprintf("fit ec=%.2g; amp=%.2g; bas=%.2g", ec, amp, bas))
							
							res	<- log10(indata$zvalues)-log10(fitVals)
							fits[[answer.lookup]] <- data.frame(concs=indata$concs, fitZvalues=fitVals, expZvals=indata$zvalues, residuals=res)
							
							show(fits)
							
							RES	<<-res
							SIGVECTOR <<-sigvector
							
							N			<- length(indata$zvalues)
							chi2		<- sum(res^2*N/indata$sigmas^2)
							deg.free	<- N-length(out$estimate)
                        
							 	
							pval	<- pgamma(0.5*chi2, 0.5*deg.free, 1, lower.tail = TRUE, log.p = FALSE)
							pvals[x.idx,y.idx]	<- pval
							
							show(sprintf("p-value = %.2g",pval))
							show(sprintf("Chi2 = %.2g", chi2))
							
							RE	<- try(sqrt( diag( 2*out$minimum/(deg.free)*solve(out$hessian)) )/out$estimate)
							if(is.numeric(RE)){
								EC.RE[x.idx,y.idx]	<- RE[3]
								Amp.RE[x.idx,y.idx]	<- RE[1]
								if(length(out$estimate==4)){
									Hill.RE[x.idx,y.idx]	<- RE[4]
								}
							}
							show("RE")
							show(RE)
							show("Values stored")
						}
						
					}else{
						# Matches check for bin fullness, need method to estimate amplitude for fewer than 4 files
						
						show(c("x,y index -", x.idx, y.idx))
						show("Too many insufficiently full bins, fitting not attempted")
						Amps[x.idx,y.idx]		<- NA
						models[[answer.lookup]]	<- NA
						Method[x.idx,y.idx]		<- "None"
						
					}
				}
			}
		}
		show(sprintf("hill fitting succeeded %i out of %i attempts", hillwin, hilltot))
		Answers[[titr.idx]] <- list(bas=Bas, models=models, ec=EC50s, amp=Amps,hill=Hill, pvals=pvals,
			fitdata=FinalZData, method=Method, drc.output=drc.output, ecRE=EC.RE, ampRE=Amp.RE, hillRE=Hill.RE)
		Answers[[titr.idx]]$models[["outside fit area"]] <- 0
	}
	
	show("FinalAnalysis Complete")
	
	ScsAnswers <<- Answers
	return(Answers)
	
}

