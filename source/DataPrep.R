DataPrep <- function(Vars,Track,Options, Bins, titr.current=1){
	show("DataPrep")
	#show(Bins)
	First 				<- Track$FirstDataprep
	Track$FirstDataprep <- F
	num.doses 			<- Vars$num.doses
	t3.dataprep 		<- tktoplevel()
	tkwm.title(t3.dataprep, "Select Parameters")
	
	#### TCL / TK VARIABLES #########################################################
	
	channel.labels	<- list()
	X.radio		<- list()
	Y.radio		<- list()
	Z.radio		<- list()
	name.entries.tk <- list()
	name.entry	<- list()
	rbval.tc	<- list()
	loglin.tc	<- list()
	
	if (First==T){
		# Set up 
		rbval.tc[["x"]]	<- tclVar(1)
		rbval.tc[["y"]]	<- tclVar(2)
		rbval.tc[["z"]]	<- tclVar(3)
		loglin.tc[["x"]]	<- tclVar(1)
		loglin.tc[["y"]]	<- tclVar(1)
		loglin.tc[["z"]]	<- tclVar(1)
		titr.current	<- 1
		Track[["titr.check"]]		<- c(1,rep(0,times=Vars$num.titr-1))
			#Use this to determine if values have been entered for a particular titration
			
		for(titr.idx in 1:Vars$num.titr){
			channels	<- names(Vars$Datafiles[[titr.idx]][[1]])
			
			Vars$Channels.orig[[titr.idx]] 	<- channels
				#These are the names of the channels from the header of the exported txt file
			
			Vars$Channels[[titr.idx]]		<- channels
				# Names entered by user
			
			Vars$xyzLabels[[titr.idx]]		<- list(x="",y="",z="")
			
			# Set up a tracking list to use when coming back
			#Vars$Chan.track[[titr.idx]]		<- list() # Delete after next run of Main.R
			#length(Vars$Chan.track[[titr.idx]])	<- length(channels)
			#names(Vars$Chan.track[[titr.idx]])	<- channels
			#for(ch.idx in 1:length(channels)){
			#	Vars$Chan.track[[titr.idx]][[ch.idx]]	<- list(xyz="none", islog=T, label=channels[ch.idx], original=channels[ch.idx])
			#}
			#Vars$Chan.track[[titr.idx]][[1]]$xyz	<- "x"
			#Vars$Chan.track[[titr.idx]][[2]]$xyz	<- "y"
			#Vars$Chan.track[[titr.idx]][[3]]$xyz	<- "z"
		}
		
	}else{
		for (xyz in c("x","y","z")){
			idx <- match(Vars$xyzLabels[[titr.current]][[xyz]],Vars$Channels[[titr.current]])
			rbval.tc[[xyz]] <- tclVar(idx)
			loglin.tc[[xyz]] <- tclVar(as.numeric((Vars$islog[[xyz]])))
		}
	}
	
	#show("Initial entry setup complete")
	
	
	#### DATA MANAGEMENT FUNCTIONS ###################################################
	
	StoreVars <- function(titr.current) {
		
		#show(sprintf("Storing Vars for titr %i", titr.current))
			
		#Take the names from the entries to store in Vars$Channels
		for (ch.idx in 1:length(Vars$Channels.orig[[titr.current]])){
			Vars$Channels[[titr.current]][[ch.idx]]			<<- tclvalue(name.entries.tk[[ch.idx]])
			#Vars$chan.track[[titr.current]][[ch.idx]]$label	<<- tclvalue(name.entries.tk[[ch.idx]])
		}

		
		# Take x, y and z selections and log choices, store in xyzLabels
		for (xyz in c("x","y","z")){
			Vars$islog[[xyz]]	<<- as.logical(as.numeric(tclvalue(loglin.tc[[xyz]])))
			if (!Vars$islog[[xyz]]) {
				Options$plots[[paste("max",xyz, sep="")]] <<- 1000
			}
			idx <- as.numeric(tclvalue(rbval.tc[[xyz]]))
			#show(sprintf("%s column index is %i", xyz, idx))
			Vars$xyzLabels[[titr.current]][[xyz]]	<<- Vars$Channels[[titr.current]][[idx]]
		}
		
		# Propagate the names and selections from titration 1 to the rest, unless they've already been named
		if(titr.current==1){
			for(titr.idx in which(Track$titr.check!=1)){	
				Vars$Channels[[titr.idx]] 	<<- Vars$Channels[[1]]
				Vars$xyzLabels[[titr.idx]]	<<- Vars$xyzLabels[[1]]
			}
		}
		
		
		ScsVars<<-Vars
		#show("Stored")
	}

	
	# TITRATION NAVIGATION  #############################################################################		
	titr.nav.fr	<- tkframe(t3.dataprep)
	OnNextTitr	<- function(){
		StoreVars(titr.current)
		Track$titr.check[titr.current]<-1
		titr.current<-titr.current+1
		if(titr.current>Vars$num.titr){titr.current<-1}
		
		tkdestroy(t3.dataprep)
		DataPrep(Vars,Track,Options,Bins,titr.current)
	}
	
	OnPrevTitr	<- function(){
		StoreVars(titr.current)
		Track$titr.check[titr.current]<-1
		titr.current<-titr.current-1
		if(titr.current==0){titr.current<-Vars$num.titr}
		
		tkdestroy(t3.dataprep)
		DataPrep(Vars,Track,Options,Bins,titr.current)
	}
	#show(titr.current)
	current.label.tk<- tclVar(Vars$prefix[[titr.current]])
	next.titr.but	<- tkbutton(titr.nav.fr, text=" >> ", command=OnNextTitr)
	prev.titr.but	<- tkbutton(titr.nav.fr, text=" << ", command=OnPrevTitr)
	titr.curr.label	<- tklabel(titr.nav.fr, text=Vars$prefix[[titr.current]])
	tkgrid(prev.titr.but, tklabel(titr.nav.fr, text="  "), titr.curr.label, tklabel(titr.nav.fr, text="  "), next.titr.but)
		
	channels.fr	<- tkframe(t3.dataprep)

	## Set up header
	channel.head<- tklabel(channels.fr, text="   Channel   ")
	name.head	<- tklabel(channels.fr, text="   Label   ")
	X.head		<- tklabel(channels.fr, text="  X  ")
	Y.head		<- tklabel(channels.fr, text="  Y  ")
	Z.head		<- tklabel(channels.fr, text="  Z  ")
	spacer		<- tklabel(channels.fr, text="   ")
	tkgrid(spacer)
	tkgrid(channel.head, name.head, X.head, Y.head, Z.head, spacer)
	
	
	# Set up rows with channels gleaned from header of input files,
	# a box to enter an alternate label and radio buttons to specify
	# whether each channel will be X, Y or Z
	
	for (channel.idx in 1:length(Vars$Channels.orig[[titr.current]])) {
		#show("filling channels")
		X.radio[[channel.idx]] 		<- tkradiobutton(channels.fr)
		Y.radio[[channel.idx]] 		<- tkradiobutton(channels.fr)
		Z.radio[[channel.idx]] 		<- tkradiobutton(channels.fr)
		channel.labels[[channel.idx]]	<- tklabel(channels.fr, text=Vars$Channels.orig[[titr.current]][channel.idx])
		if (First==T){
			name.entries.tk[[channel.idx]] 	<- tclVar(Vars$Channels.orig[[titr.current]][channel.idx])
		}else{
			if (channel.idx <= length(Vars$Channels[[titr.current]])){
			name.entries.tk[[channel.idx]] <- tclVar(Vars$Channels[[titr.current]][[channel.idx]])
			}
			else {
				name.entries.tk[[channel.idx]] <- tclVar(Vars$Channels.orig[[titr.current]][channel.idx])
			}
		}		
		name.entry[[channel.idx]] 	<- tkentry(channels.fr, textvariable=name.entries.tk[[channel.idx]], width = 15)
		tkconfigure(X.radio[[channel.idx]],variable=rbval.tc[["x"]],value=channel.idx)
		tkconfigure(Y.radio[[channel.idx]],variable=rbval.tc[["y"]],value=channel.idx)
		tkconfigure(Z.radio[[channel.idx]],variable=rbval.tc[["z"]],value=channel.idx)
		tkgrid(channel.labels[[channel.idx]], name.entry[[channel.idx]], X.radio[[channel.idx]], Y.radio[[channel.idx]], Z.radio[[channel.idx]], tklabel(channels.fr, text="   "))	
		#show(1)
		tkgrid.configure(name.entry[[channel.idx]], sticky="ew")
	}

	## Create log/lin checkboxes
	xlog.chkbut <- tkcheckbutton(channels.fr)
	ylog.chkbut <- tkcheckbutton(channels.fr)             
	zlog.chkbut <- tkcheckbutton(channels.fr)
	tkconfigure(xlog.chkbut, variable=loglin.tc$x)
	tkconfigure(ylog.chkbut, variable=loglin.tc$y)
	tkconfigure(zlog.chkbut, variable=loglin.tc$z)
	loglabel <- tklabel(channels.fr, text=" Log:")
	tkgrid(tklabel(channels.fr, text="   "), loglabel, xlog.chkbut, ylog.chkbut, zlog.chkbut)
	#show(2)
	tkgrid.configure(loglabel, sticky="e")
	
	
	
	
	
	
	## PROGRAM NAVIGATION  #########################################################
	
	# Specify functions of navigation buttons
	
	
	OnNextToInitPlot <- function(){
		StoreVars(titr.current)
		finished	<- F
		
		for(titr.idx in 1:Vars$num.titr){
			
			if (length(Vars$xyzLabels[[titr.idx]])==length(unique(Vars$xyzLabels[[titr.idx]]))){
				# Check that unique x, y and z value were chosen
				#show(Vars$xyzLabels[[titr.idx]])
				
				x.idx <- match(Vars$xyzLabels[[titr.idx]][["x"]],Vars$Channels[[titr.idx]])
				y.idx <- match(Vars$xyzLabels[[titr.idx]][["y"]],Vars$Channels[[titr.idx]])
				z.idx <- match(Vars$xyzLabels[[titr.idx]][["z"]],Vars$Channels[[titr.idx]])
				
				#show(sprintf("x.idx = %i; y.idx=%i; z.idx=%i",x.idx,y.idx,z.idx))
				
				parseFile	<- function(dfile){
					## File Parsing function - extracts x y and z columns
					x <- dfile[[x.idx]]
					y <- dfile[[y.idx]]
					z <- dfile[[z.idx]]
					
					d <- data.frame(x,y,z)
					if(Vars$islog$x){
						d <- subset(d, (x>=1))
						d$x <- log10(d$x)
					}
					if(Vars$islog$y){
						d <- subset(d, (y>=1))
						d$y <- log10(d$y)
					}
					return(d)
				}
				
				
				# Read files, parse out the x, y and z columns
				separator	<- ifelse(grepl("csv", Vars$Filenames[[1]][[1]]), ",", "\t") # Prep separator for csv or text files
				
				parsetime <- system.time(
				Vars$Datas[[titr.idx]] <- foreach (dose.idx=1:num.doses) %do% {
					infile	<- read.delim(paste(Vars$Filenames[[titr.idx]][[dose.idx]]),header=TRUE, colClasses="numeric", sep=separator)
					parsed	<- parseFile(infile)
					return(parsed)
				}
				)[3]
				show(sprintf("reading and parsing time was %.2f sec", parsetime))
				
				ScsVars 	<<- Vars
				ScsTrack 	<<- Track
				tkdestroy(t3.dataprep)
				finished	<- T
			}else{
				show("XYZ error")
				tkmessageBox(message=" Please choose distinct parameters for X, Y and Z. ", icon="warning", type="ok")
				finished<-F
				break()
			}
		}
		if (finished){
			ScsVars	<<- Vars
			#show(Bins)
			BatchPrep(Vars,Track,Options,Bins)
		}
		
		
	}
	

	
	OnBackToEntry <- function(){

		for(titr.idx in 1:Vars$num.titr){
			StoreVars(titr.idx)	
		}
		
		#show("Prior to return to file entry")
		tkdestroy(t3.dataprep)
		ScsTrack 	<<- Track
		EnterFiles(Vars,Track,Options,Bins)
	}
	
	
	# Specify navigation buttons
	prog.nav.fr	<- tkframe(t3.dataprep)
	Next.InitialPlot.but 	<- tkbutton(prog.nav.fr,text="  Next  ", command=OnNextToInitPlot)
	#Next.Batch.but			<- tkbutton(prog.nav.fr,text=" Batch  ", command=OnBatch)
	Back.ToEntry.but 		<- tkbutton(prog.nav.fr,text="  Back  ", command=OnBackToEntry)
	
	tkbind(t3.dataprep, "<Return>", OnNextToInitPlot)
	
	# Layout of buttons
	tkgrid(tklabel(prog.nav.fr, text="   "))
	tkgrid(Back.ToEntry.but,Next.InitialPlot.but) #, Next.Batch.but)
	#show(3)
	tkgrid.configure(Next.InitialPlot.but, sticky="e")
	tkgrid.columnconfigure(prog.nav.fr, 1, weight=1)
	
	
	##########  INSERT FRAMES INTO WINDOW #############################################
	tkgrid(tklabel(t3.dataprep, text="   "))
	tkgrid(titr.nav.fr)
	tkgrid(ttkseparator(t3.dataprep, orient="horizontal"), sticky="ew", columnspan=1)
	tkgrid(channels.fr)
	tkgrid(ttkseparator(t3.dataprep, orient="horizontal"), sticky="ew", columnspan=1)
	tkgrid(prog.nav.fr)
	tkgrid(tklabel(t3.dataprep, text="   "))
}
