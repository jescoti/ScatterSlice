# Fit hill equations, using concentration values entered previously
# and calculated mean Z values

# indata will be a data frame of four columns, $zvalues, $concs, $sigmas, $N
# Options used are Options$analysis$hillchoice == "fixed" or "free" (default is "free", 4-parameter hill fit)

HillFit <- function(indata, Options=list(analysis=list(hillchoice="free"))){
	#show("Running HillFit")
	
	if(max(indata$zvalues, na.rm=T) == min(indata$zvalues, na.rm=T)){
			outdata$method <- "fail"
			
	}else{
	
		#############################################################################
		## Set initial values for nlm function
		p0=c(0,0,0) # (base, amp, ec50)
		p0[1]=indata$zvalues[which.min(indata$concs)] # base
		p0[2]=indata$zvalues[which.max(indata$concs)]-indata$zvalues[which.min(indata$concs)]# amplitude
		p0[3]=indata$concs[which.min(abs((indata$zvalues-min(indata$zvalues, na.rm=T))/(max(indata$zvalues, na.rm=T)-min(indata$zvalues, na.rm=T))-0.5))]
		#show(indata)
		#show(sprintf("Estimates: EC50=%.2g, amp=%.2g, base=%.2g", p0[3], p0[2], p0[1]))
		
		############################################################################
		hillfitfunc	<- function(npar){
			outdata <- list()
			outdata$method	<- "fail"
			
			guesses	<- p0
			if(npar==3){
				# Fit with hill coeff=1
				fn <- function(p) sum((log10(indata$zvalues)-log10(p[1]+p[2]*indata$concs/(indata$concs+p[3])))^2*indata$N/indata$sigmas^2)
			}else{
				#show("Attempting fit with free hill coefficient")
				fn <- function(p) sum((log10(indata$zvalues)-log10(p[1]+p[2]*indata$concs^p[4]/(indata$concs^p[4]+p[3]^p[4])))^2*indata$N/indata$sigmas^2)
				guesses[4]<-1
			}
			
			## Fit
			hilleq 	<- try(nlm(fn, p = guesses, hessian = TRUE, steptol=1e-18, gradtol=1e-8, iterlim=1000))
			
			# hillnls	<- try(nls())
			
			HILLEQ <<- hilleq
			#show(hilleq)
			if(class(hilleq)!="try-error") {
				if(hilleq$code<4){
					outdata$method	<- paste("Hill",npar, sep="")
					outdata$hilleq	<- hilleq
				}
			}
			return(outdata)
		}
		
		if (Options$analysis$hillchoice =="fixed") {
			outdata	<- hillfitfunc(3)
		}else{
			outdata	<- hillfitfunc(4)
			if(outdata$method=="fail"){
				outdata	<- hillfitfunc(3)
			}
		}
	}
	#show("Finished Hill Fitting")
	return(outdata)
}

