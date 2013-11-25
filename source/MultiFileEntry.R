# Use this to get files:

startBox	<- function(Vars,Track,Options,Bins){
	#show("startBox")
	t1.start	<- tktoplevel()
	tkwm.title(t1.start, "Start")
	loadfiles.fr	<- tkframe(t1.start)
	j=0
	
	# Allows loading of a saved analysis
	OnLoad	<- function(){
		savedFile <- tclvalue(tkgetOpenFile(filetypes="{{ScatterSlice files} {.scs}}"))
		tkdestroy(t1.start)
		load(savedFile)
		Analysis3D(saveData$Answers,saveData$Vars,saveData$Options,saveData$Bins,saveData$Track)
	}
	load.label	<- tklabel(loadfiles.fr,text="  Load saved analysis ")
	load.but	<- tkbutton(loadfiles.fr,text=" Select...", command=OnLoad)
	tkgrid(load.label, load.but)
	
	newfiles.fr		<- tkframe(t1.start)
	new.label		<- tklabel(newfiles.fr, text="Enter new files")
	num.titr.label	<- tklabel(newfiles.fr, text=" Number of titrations: ")
	
	num.titr.tk		<- tclVar(Vars$num.titr)
	num.titr.entry	<- tkentry(newfiles.fr, width="10", textvariable=num.titr.tk)
	num.doses.label	<- tklabel(newfiles.fr, text=" Doses per titration: ")
	num.doses.tk	<- tclVar(Vars$num.doses)
	num.doses.entry	<- tkentry(newfiles.fr, width="10", textvariable=num.doses.tk)
	
	OnNewFiles	<- function(){
		Vars$num.titr	<- as.numeric(tclvalue(num.titr.tk))
		Vars$num.doses	<- as.numeric(tclvalue(num.doses.tk))
		
		if(is.numeric(Vars$num.titr)&is.numeric(Vars$num.doses)){
			if(Vars$num.doses <= 4){
				Options$analysis$method	<- "subtract"
			}
			
			tkdestroy(t1.start)
			for(section.idx in c("Datafiles","Filenames","Concs","Datas","units","Control.File","prefix")){
				length(Vars[[section.idx]])=Vars$num.titr
			}

			EnterFiles(Vars,Track,Options,Bins)
		}
	}
	new.but		<- tkbutton(newfiles.fr,text="  Go  ", command=OnNewFiles)
	
	tkbind(t1.start,"<Return>", OnNewFiles)
	tkgrid(new.label, sticky="w")
	tkgrid(num.titr.label, num.titr.entry)
	tkgrid(num.doses.label, num.doses.entry)
	tkgrid(new.but)
	
	tkgrid(loadfiles.fr)
	tkgrid(ttkseparator(t1.start, orient="horizontal"), sticky="ew", columnspan=1)
	tkgrid(newfiles.fr)
	
}





#### File selection and concentration specification dialog #####################
#show(sprintf("FRAME = %d", sys.nframe()))

EnterFiles <- function(Vars,Track,Options,Bins, titr.current=1){
	#show("EnterFiles")
	#show(Bins)
	First    		 	<- Track$FirstEntry
	Track$FirstEntry 	<- F
	num.doses 			<- Vars$num.doses
	num.titr			<- Vars$num.titr
	titr.current		<- 1
	
	#	Set up basic window

	#	Create variables and fill in lists to allow construction of
	#	file entry window.
	filename 		<- tclVar("")
	control.radio 	<- list()
	conc.entries 	<- list()
	file.buttons 	<- list()
	file.enter 		<- list()
	file.vals 		<- list()
	file.entries 	<- list()
	rbValue			<- tclVar(1)
	units.tk 	<- tclVar(Vars$units[[titr.current]])
	prefix.tk	<- tclVar(Vars$prefix[[titr.current]])
	
	
	
	## Trying to get the folder to stick
	#show(sprintf("FRAME = %d", sys.nframe()))
	
	
	env	<- environment()
	
	
	
	#### If this is the first time through, set up as appropriate.
	#### If not, use the previously-entered data. 
	if(First==T){
		#show(c("Running first-time entry setup"))
		rbValue 		<- tclVar(1)
		
		for(titr.idx in 1:num.titr){
			#show(sprintf("titration %i",titr.idx))
			for (dose.idx in 1:num.doses) {
				#show(sprintf("dose %i",dose.idx))
				file.vals[[dose.idx]] <- tclVar(0)
				file.entries[[dose.idx]] <- tclVar("...")
			
				
				Vars$Filenames[[titr.idx]][[dose.idx]] <- "..."
				Vars$Concs[[titr.idx]][[dose.idx]] <- 0
				Vars$units[[titr.idx]]	<- "pM"
				Vars$prefix[[titr.idx]]	<- as.character(titr.idx)
				Vars$Control.File[[titr.idx]]	<- 1
				Bins$params[titr.idx]	<- "0"
	}}}
	
	#	Resizes data variables to fit new number of files if this window
	#	has been accessed after clicking "Back" button and changing number
	# 	of files to input.
	if (First==FALSE){
		#show(paste("rerun titration=",titr.current))
		file.entries  	<- lapply(Vars$Filenames[[titr.current]],tclVar) #filename
		file.vals		<- lapply(Vars$Concs[[titr.current]], tclVar)
		rbValue <- tclVar(unique(min(c(Vars$Control.File[[titr.current]], Vars$num.doses))))
		
		if (num.doses>length(file.entries)){
			for (dose.idx in ((length(file.entries)+1):num.doses)){
				file.vals[[dose.idx]] <- tclVar(0)
				file.entries[[dose.idx]] <- tclVar("...")
				for(titr.idx in 1:num.titr){
					Vars$Filenames[[titr.idx]][[dose.idx]] <- "..."
					Vars$Concs[[titr.idx]][[dose.idx]] <- 0
		}}}
		
		if (num.doses<length(file.entries)){
			length(file.entries)=num.doses
			length(file.vals)=num.doses
			for(titr.idx in 1:num.titr){
				length(Vars$Filenames[[titr.idx]])=num.doses
				length(Vars$Concs[[titr.idx]])=num.doses
		}}
	}
	
	
	StoreVars <- function(titr.current) {
	# Function used to store in global variables
		#show(c("StoreVars"))
		#show(sprintf("FRAME = %d", sys.nframe()))
		
		Vars$Control.File[[titr.current]] <<- as.numeric(tclvalue(rbValue))
		
		
		errors <- list()
		errors$filename	<- 0
		errors$conc		<- 0

		for (i in 1:num.doses){
			filename 	<<- as.character(tclvalue(file.entries[[i]])) #filename
			#show("filename")
			#show(filename)
			if((grepl(".txt", filename) | grepl(".csv", filename)) < 1){
				show("Wrong formats detected")
				errors$filename <- 1
			}
			
			Vars$Filenames[[titr.current]][[i]] <<- filename
			
			fileconc	<- as.numeric(as.character(tclvalue(file.vals[[i]]))) #concentration
			#show("fileconc")
			#show(fileconc)
			if (is.na(fileconc)) {
				#show("Improper conc value detected")
				errors$conc	<- 1
			}
			else {
				Vars$Concs[[titr.current]][i] <<- fileconc
			}

			
		}
		
		Vars$errors		<<- errors
		Vars$units[[titr.current]]	<<- as.character(tclvalue(units.tk))
		Vars$prefix[[titr.current]]	<<- as.character(tclvalue(prefix.tk))
		ScsVars 		<<- Vars	
		ScsTrack 		<<- Track
		#show("Stored Vars")
	}


	
	MakeEntryWindow <- function(titr.current) {
		#show("MakeEntryWindow")
		#show("Vars$Concs in MakeEntryWindow")
		#show(Vars$Concs)
		#show("Vars$Control.File")
		#show(Vars$Control.File)
		#show(Bins)
		#show("top of MakeEntryWindow")
		#show(sprintf("Filefolder = %s", Filefolder))
		#show(sprintf("FRAME = %d", sys.nframe()))
		t2.entry <- tktoplevel()
		
		tkwm.title(t2.entry,"File Entry")
		
		#### Functions for buttons ####
		
		OnNextStep <- function(){
			#show("OnNextStep")
			#show(Bins)
			StoreVars(titr.current)
			headererror	<- F
			errmsg 		<- character(0)
			
			#If there are no obvious errors, try to read files, checking for errors along the way
			if((Vars$errors$filename + Vars$errors$conc) == 0){ 
				separator	<- ifelse(grepl("csv", Vars$Filenames[[1]][[1]]), ",", "\t") # Prep separator for csv or text files
				for (titr.idx in 1:num.titr){
					
					# Read beginning of files, do loads of error checking
					infiles <- foreach (dose.idx=1:num.doses) %do% {
						try(read.delim(paste(Vars$Filenames[[titr.idx]][[dose.idx]]),header=TRUE, colClasses="numeric", sep=separator, nrows=5))
					}
					
					# Error checking
					for(dose.idx in 1:length(infiles)){
						infile	<- infiles[[dose.idx]]
						if (attr(infile, "class") == "try-error"){
							errmsg 	<- sprintf("Could not read file %i in titration %i.  Error was as follows: \n %s \n Please check that csv or plain text (.csv or .txt) file of correct format was entered", dose.idx, titr.idx, infile[1])
							break
						}else{  # if no try-error, check that files have the same columns
							Vars$Datafiles[[titr.idx]][[dose.idx]] <- infile
							# store stats of first file for comparison 
							if (dose.idx == 1) {
								headernames 	<- names(infile)
								headerlength	<- length(headernames)
							}else{ # for files 2 to num.doses
								if (length(infile) != headerlength){
									errmsg <- "Error:  Files must have same number of columns"
								}else{ #if same length, check that columns are the same
									if (sum(headernames != names(infile)) > 0) {
										if (sum(is.na(as.numeric(a))) == 0) {
											headererror <- T
										}else{
											errmsg	<- "Error:  Column headings differ between files"
					}}}}}}}
				if (headererror) {
					tkmessageBox(message="Input files may not have column titles. xyzLabels were not compared for consistency between files.", icon="warning", type="ok")
				}
				if (length(errmsg) > 0) {
					tkmessageBox(message=errmsg,icon="warning",type="ok")
				}
				else { # if no errors, store variables and proceed to next step
					ScsVars 	<<- Vars
					ScsTrack 	<<- Track
					DataPrep(Vars,Track,Options,Bins, titr.current=1)
					tkdestroy(t2.entry)
				}
			}
						
			else { #if there are obvious errors, #show a message box
				if((Vars$errors$filename + Vars$errors$conc)==2){
					errmsg <- "2 errors: \n 1.Input only .txt or .csv files. \n 2.Enter only numbers in concentration fields."
				}
				else {
					if (Vars$errors$filename == 1) {
						errmsg <- "Input only plain .txt or .csv files."
						}
					if (Vars$errors$conc == 1) {
						errmsg <- "Enter only numbers in concentration fields"
						}
				}
				tkmessageBox(message=errmsg,icon="warning",type="ok")
			}
		}
				
	
		
		
		OnBack <- function(){
			for(titr.idx in 1:num.titr){
				StoreVars(titr.idx)
			}
			tkdestroy(t2.entry)
			startBox(Vars, Track,Options,Bins)
		}
		
		updateTclVars	<- function(titr.current){
			#show("updateTclVars")
			#show(Vars$Control.File)
			file.entries<<- lapply(Vars$Filenames[[titr.current]],tclVar)
			file.vals	<<- lapply(Vars$Concs[[titr.current]],tclVar)
			rbValue		<<- tclVar(Vars$Control.File[[titr.current]])
			units.tk 	<<- tclVar(Vars$units[[titr.current]])
			prefix.tk	<<- tclVar(Vars$prefix[[titr.current]])
		}
			
		OnNextTitr	<- function(){
			StoreVars(titr.current)
			titr.current<-titr.current+1
			if(titr.current>Vars$num.titr){titr.current<-1}
			updateTclVars(titr.current)
			tkdestroy(t2.entry)
			MakeEntryWindow(titr.current)
		}
		
		OnPrevTitr	<- function(){
			StoreVars(titr.current)
			titr.current<-titr.current-1
			if(titr.current==0){titr.current<-Vars$num.titr}
			updateTclVars(titr.current)
			tkdestroy(t2.entry)
			MakeEntryWindow(titr.current)
		}
		
		OnConc <- function(){
			#show("**OnConc *******************************************************************")
			StoreVars(titr.current)
			#show(Vars$Control.File)
			t.concs	<- tktoplevel()
			tkwm.title(t.concs, "Concentrations")
			tkgrid(tklabel(t.concs, text= " Fill concentration values "))
			
			# Entry box for top value
			topval.fr	<- tkframe(t.concs)
			topval.tk	<- tclVar(max(Vars$Concs[[titr.current]]))
			topval.entry	<- tkentry(topval.fr, textvariable=topval.tk, width=10)
			topval.label	<- tklabel(topval.fr, text=" Top Value:")
			tkgrid(topval.label,topval.entry)
			tkgrid(topval.fr)
			
			#Radio buttons to select dilution
			dil.sel.fr	<- tkframe(t.concs)
			rbval.tk	<- tclVar(1)
			by10.rb		<- tkradiobutton(dil.sel.fr)
			bysqrt10.rb	<- tkradiobutton(dil.sel.fr)
			byother.rb	<- tkradiobutton(dil.sel.fr)
			tkconfigure(by10.rb, variable=rbval.tk, value=1)
			tkconfigure(bysqrt10.rb, variable=rbval.tk, value=2)
			tkconfigure(byother.rb, variable=rbval.tk, value=3)
			
			# Labels for radio buttons
			by10.lb	<- tklabel(dil.sel.fr, text=" 10 ")
			bysqrt10.lb	<- tklabel(dil.sel.fr, text=" sqrt(10) ")
			byother.lb	<- tklabel(dil.sel.fr, text=" other ")
			
			#Entry box to input "other" dilution
			otherdil.tk	<- tclVar(2)
			otherdil.entry	<- tkentry(dil.sel.fr, textvariable=otherdil.tk, width=5)
			
			tkgrid(tklabel(dil.sel.fr, text=" Dilution "))
			tkgrid(by10.lb, by10.rb, tklabel(dil.sel.fr, text="   "))
			tkgrid(bysqrt10.lb, bysqrt10.rb, tklabel(dil.sel.fr, text="   "))
			tkgrid(byother.lb, byother.rb, otherdil.entry)
			tkgrid(tklabel(dil.sel.fr, text="   "))
			
			topbot.lb	<- tklabel(dil.sel.fr, text="  Highest value:")
			topbot.tk	<- tclVar(1)
			top.rb	<- tkradiobutton(dil.sel.fr)
			bot.rb	<- tkradiobutton(dil.sel.fr)
			tkconfigure(top.rb, variable=topbot.tk, value=1)
			tkconfigure(bot.rb, variable=topbot.tk, value=2)
			spacer <- tklabel(dil.sel.fr, text="  ")
			
			tkgrid(spacer, tklabel(dil.sel.fr, text="  Top  "), tklabel(dil.sel.fr, text="  Bottom  "))
			tkgrid(topbot.lb, top.rb, bot.rb)
			
			
			tkgrid(dil.sel.fr)
			
			
			
			fillseq <- function(dilution, topval,Vars, topbot){
				#show("Running fillseq() **********************************************************")
				#show(sprintf("Creating %i (Num.doses) dilution by 10^1/%i (dilution), starting at %.2f (topval)", Vars$num.doses-1, dilution, topval))
				if (dilution==3){
					dilution	<- as.numeric(tclvalue(otherdil.tk))
				}else{
					dilution	<- 10^(1/dilution)
				}
				
				dil.seq	<- topval/dilution^(0:(Vars$num.doses-2))
				dil.seq <- signif(dil.seq,4)
				#show("Dilution sequence calculated")
				#show(dil.seq)
				if(topbot==2){
					dil.seq	<- sort(dil.seq)
				}
				return(dil.seq)
			}
			
			
			
			# Function and button to specify sequence
			OnGenerate	<- function(){
				#show("OnGenerate() *************")
				#show(c("control file", Vars$Control.File))
				dilution	<- as.numeric(tclvalue(rbval.tk))
				topval		<- as.numeric(tclvalue(topval.tk))
				topbot		<- as.numeric(tclvalue(topbot.tk))
				doses		<- fillseq(dilution, topval, Vars, topbot)
				for (titr.idx in 1:Vars$num.titr){
					ctl	<- Vars$Control.File[[titr.idx]]
					Vars$Concs[[titr.idx]][-ctl]	<<- doses
					#Vars$Control.File	<<- Vars$Control.File
					updateTclVars(titr.idx)
				}
				
				tkdestroy(t.concs)
				tkdestroy(t2.entry)
				#show("Vars$Concs in OnGenerate")
				#show(Vars$Concs)
				EnterFiles(Vars,Track,Options,Bins, titr.current)
				#MakeEntryWindow(titr.current)
			}
			
			
			specify.seq.but	<- tkbutton(t.concs, text="Insert Sequence", command=OnGenerate)
			tkgrid(specify.seq.but)
		}
		
		
		
		
		#### Window objects
		
		# Titration info Frame
		titr.info.fr	<- tkframe(t2.entry)
		units.label <- tklabel(titr.info.fr, text=" units ")
		units.entry	<- tkentry(titr.info.fr, textvariable= units.tk, width=15)
		prefix.label	<- tklabel(titr.info.fr, text=" Titration title: ")
		prefix.entry	<- tkentry(titr.info.fr, textvariable=prefix.tk, width=40)
		nextTitr.but	<- tkbutton(titr.info.fr, text=" >> ", command=OnNextTitr)
		prevTitr.but	<- tkbutton(titr.info.fr, text=" << ", command=OnPrevTitr)
		tkgrid(prevTitr.but, prefix.label, prefix.entry, units.label, units.entry, nextTitr.but)
		
		#tkgrid.configure(units.label, sticky="e")
		#tkgrid.columnconfigure(t2.entry, 1, weight=1)
		
		
		# File entry frame
		filentry.fr	<- tkframe(t2.entry)
		Radio.head	<- tklabel(filentry.fr,text="  Control:  ")
		But.head 	<- tklabel(filentry.fr,text="         ")
		Filename.head	<- tklabel(filentry.fr,text="File")
		Conc.head.but	<- tkbutton(filentry.fr,text="Concentration", command=OnConc)
		tkgrid(But.head, Filename.head, Radio.head, Conc.head.but)
		
		#Radio.head	<- tklabel(t2.entry,text="  Control:  ")
		#But.head 	<- tklabel(t2.entry,text="         ")
		#Filename.head	<- tklabel( t2.entry,text="File")
		#Conc.head.but	<- tkbutton(t2.entry,text="Concentration", command=OnConc)
		#tkgrid(But.head, Filename.head, Radio.head, Conc.head.but)
		
		
		# Build rows
		MakeRow <- function(row.idx, filenamewidth){
						
			ChooseFile <- function(){
				tkdestroy(t2.entry)
				#show("ChooseFile")
				#show(sprintf("Filefolder = %s", Filefolder))
				#show(sprintf("FRAME = %d", sys.nframe()))
				
				fullnames	<- tkgetOpenFile(multiple=T, filetypes="{{plain text csv} {.txt .csv}}", initialdir=Filefolder)
				fullnames	<- tclvalue(fullnames)
				
				
				
				FILES	<<- fullnames #Save to global for troubleshooting
				filename.vec	<- c()
				#show("fullnames = all files selected")
				#show(fullnames)
				if(nchar(fullnames)!=0){	# Checks to be sure that the user didn't cancel
					
					##########################################################################
					## This section breaks up the returned filenames into individual files
					## since tkgetOpenFile() returns everything as one long string.
					## Unfortunately, the formatting is bizarre, so this is tricky.
					
					#show("Splitting file string")
					allwords 	<- strsplit(fullnames, split="/", fixed=TRUE)[[1]]
					drive	<- ""
					if(.Platform$OS.type=="windows"){
						drive	<- allwords[1]
						driveltrs	<- strsplit(drive, split="")[[1]]
						driveltrs	<- driveltrs[driveltrs != "{"]
						drive	<- paste(driveltrs, collapse="")
						#allwords	<- allwords[2:length(allwords)]
						}
					uqwords		<- unique(allwords)
					
					# determine whether csv or txt files were chosen
					if(grepl("csv", allwords[length(allwords)])){filetype	<- "csv"}
					if(grepl("txt", allwords[length(allwords)])){filetype 	<- "txt"}
					
					if(length(allwords)==length(uqwords)){
						## If one file is chosen, just have to get rid of possible brackets

						ltrs	<- strsplit(fullnames, split="")[[1]]
						if(ltrs[1]=="{") {ltrs<-ltrs[-c(1, length(ltrs))]}
						filename	<- paste(ltrs, collapse="")
						#show(filename)
						filename.vec[1] <- filename
						
					}else{
						# If multiple files, were chosen, dealing with the formatting of the filenames is a bastard...
						
						# Count the number of repeats by countint how many times words in the filenames are repeated
						# Assumes that the folder will be repeated the most.
						repct		<- c()
						for(word.idx in 1:length(uqwords)){
							repct[word.idx]	<- sum(uqwords[word.idx]==allwords)
						}
						filenr	<- max(repct)
						path	<- paste(drive, "/", paste(c(uqwords[repct==filenr]), collapse="/"), "/", sep="")
						
						fnams	<- strsplit(fullnames, split=path, fixed=TRUE)[[1]]
						#show(fnams)
						file.idx	<- 1
						for(fnam.idx in 1:length(fnams)){
							
							ltrs	<- strsplit(fnams[fnam.idx], split="")[[1]]
							#show(ltrs)
							ltrs	<- ltrs[ltrs != "{"]
							ltrs	<- ltrs[ltrs != "}"]
							#if(.Platform$OS.type=="windows"){
							#	# filenames will end with "C:" (or D: or whatever) so strip the last two characters
							#	length(ltrs) <-length(ltrs)-2
							#	}
							
							# If the name ends with a space, strip the space out
							spaces	<- which(ltrs == " ")
							if(length(spaces)>0){
								if(max(spaces)==length(ltrs)){
									ltrs	<- ltrs[-max(spaces)]
								}
							}
							
							filename	<- paste(ltrs, collapse="")
							if(nchar(filename)>0){
								filename.vec[file.idx] <- sprintf("%s%s",path,filename)
								file.idx	<- file.idx+1
							}
						}
					}
					#########################################################################
					
					
					
					
					## Extract the location of the files and store it for ease of repeat.
					fileloc <- gregexpr(paste("[^/]*.", filetype, sep=""), filename.vec[1])[[1]]	# Extract the location of the files
					filefolder.temp <- substr(filename.vec[1], 0, fileloc-2)	# Set this as the default for next selection
					assign("Filefolder", filefolder.temp, env=.GlobalEnv)
					
					
					
					
					## Take the files chosen, and put them in the appropriate places
					## checking along the way to be sure that the correct number were chosen
					FILENAME.VEC<<-filename.vec
					num.selected<- length(filename.vec)
					num.stated	<- num.titr*num.doses
					if(num.selected>1){
						## If more than one file selected
						#show("Counting and placing files")
						titr.idx	<- titr.current
						if((titr.idx-1) * num.doses + row.idx + num.selected - 1 > num.stated){
							## Check to be sure that we're not overflowing, trim if necessary
							length(filename.vec)= num.stated - row.idx + 1 - (titr.idx-1)*num.doses
							num.selected		= num.stated - row.idx + 1 - (titr.idx-1)*num.doses
							tkmessageBox(message="More files chosen than number of doses and titrations indicated. Some files will not be loaded", icon="warning", type="ok")
						}
						
						for(filevec.idx in 1:num.selected){
							## Place files in the right slots
							#show(sprintf("placing file %d", filevec.idx))
							dose.idx	<- (row.idx+filevec.idx-1)%%num.doses
							if(dose.idx==0){dose.idx<-num.doses}
							if((row.idx+filevec.idx-1)<=num.doses){
								#show(filename.vec[filevec.idx])
								
								tclvalue(file.entries[dose.idx])	<<- filename.vec[filevec.idx]
								Vars$Filenames[[titr.idx]][[dose.idx]] <<- filename.vec[filevec.idx]
							}else {
								#show(sprintf("dose.idx=%i",dose.idx))
								if(dose.idx==1){
									titr.idx	<- titr.idx+1
								}
								Vars$Filenames[[titr.idx]][[dose.idx]] <<- filename.vec[filevec.idx]
								#show(filename.vec[filevec.idx])
								
							}	
						}
					}else{
						## If only one file selected, store it in the right place
						Vars$Filenames[[titr.current]][[row.idx]] <<- filename.vec[1]
					}
				}
				
				EnterFiles(Vars,Track,Options,Bins, titr.current)
				ScsVars <<- Vars
				#show("End of ChooseFile")
				
			}
			
			
			file.buttons[[row.idx]] 	<- tkbutton(t2.entry,text="Choose...",command=ChooseFile)
			file.enter[[row.idx]] 		<- tkentry( t2.entry, textvariable=file.entries[[row.idx]], width = 110)
			control.radio[[row.idx]] 	<- tkradiobutton(t2.entry)
			conc.entries[[row.idx]] 	<- tkentry(t2.entry, textvariable=file.vals[[row.idx]], width=15)
			
			file.buttons[[row.idx]] 	<- tkbutton(filentry.fr,text="Choose...",command=ChooseFile)
			file.enter[[row.idx]] 		<- tkentry( filentry.fr, textvariable=file.entries[[row.idx]], width = 110)
			control.radio[[row.idx]] 	<- tkradiobutton(filentry.fr)
			conc.entries[[row.idx]] 	<- tkentry(filentry.fr, textvariable=file.vals[[row.idx]], width=15)
			
			tkconfigure(control.radio[[row.idx]],variable=rbValue,value=row.idx)
			tkgrid(file.buttons[[row.idx]], file.enter[[row.idx]], control.radio[[row.idx]], conc.entries[[row.idx]])
			tkgrid.configure(file.enter[[row.idx]], sticky = "ew")
			
		}
		for (i in 1:num.doses){
			MakeRow(i, filenamewidth)
		}
		
		# Navigation frame
		navigation.fr	<- tkframe(t2.entry)
		Nextstep.but 	<- tkbutton(navigation.fr,text=" Read files and continue ", command=OnNextStep)
		Back.but 		<- tkbutton(navigation.fr,text="  Back  ", command=OnBack)
		spacer			<- tklabel(navigation.fr,text="   ")
		tkgrid(Back.but, spacer,  Nextstep.but)
		
		
		# Window layout
		tkgrid(titr.info.fr)
		tkgrid(ttkseparator(t2.entry, orient="horizontal"), sticky="ew", columnspan=1)
		tkgrid(filentry.fr)
		tkgrid(ttkseparator(t2.entry, orient="horizontal"), sticky="ew", columnspan=1)
		tkgrid(navigation.fr)
		
		#tkgrid.configure(titr.info.fr, sticky = "ew")
		#tkgrid.configure(filentry.fr, sticky = "ew")
		#tkgrid.configure(navigation.fr, sticky = "ew")
		

	#show("end of MakeEntryWindow")
	
	}
	MakeEntryWindow(titr.current)
}

