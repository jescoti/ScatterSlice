## Log plot prep functions ############################################################

loglabels	<- function(label.range){
	#show("log labeling")
	
	bottom.rnd	<- floor(label.range[1])
	top.rnd		<- ceiling(label.range[2])
	
	tickslabels	<- list()
	label.final	<- expression()
	label.seq	<- bottom.rnd:top.rnd
	
	for (i in 1:length(label.seq)){
		label.final[i]	<- as.expression(bquote(.(10)^.(label.seq[i])))
	}
	
	# Find decimal sequence
	dec	<- log10(1:9)
	decticks.df	<- merge(label.seq,dec)
	decticks	<- decticks.df[,1]+decticks.df[,2]
	decticks	<- decticks[decticks >= label.range[1]]
	decticks	<- decticks[decticks <= label.range[2]]
	decticks	<- decticks[-which(decticks %in% label.seq)]
	
	tickslabels$ticks	<- label.seq
	tickslabels$labels	<- label.final
	tickslabels$decseq	<- decticks
	
	#show(tickslabels)
	
	return(tickslabels)
}

linlabels	<- function(label.range){
	#show("lin labeling")
	bottom.rnd	<- floor(label.range[1])
	top.rnd		<- ceiling(label.range[2])
	
	tickslabels	<- list()
	label.seq	<- pretty(label.range, n=5)
	
	# Add halfway points
	
	decticks	<- label.seq[2:length(label.seq)]-(label.seq[2]-label.seq[1])/2
	
	tickslabels$ticks	<- label.seq
	tickslabels$labels	<- label.seq
	tickslabels$decseq	<- decticks
	
	#show(tickslabels)
	
	return(tickslabels)
	
}

##############################################
# niceaxes - main function for plotting nice log or linear axes, with labels at a good distance
#
# parameters:
#	xy - "x" if for the x axis, "y" if for the y axis.
#	limits - vector - range of axes, expects c(min, max)
#	axtitle - string - title of axis
#	islog - boolean - is the axis log?


niceaxes	<- function(xy, limits, axtitle="", islog=T)  {
	
	
	
	# Sets up axes and titles for log graphs, for plotting log transformed data
	
	if (islog == T){
		tickslabels	<- loglabels(limits)
	}else{
		tickslabels	<- linlabels(limits)
	}
	
	
	#show(tickslabels)
	
	whichaxis	<- ifelse(xy == "x", 1, 2)
	
	axis(whichaxis, at=tickslabels$ticks, labels=F) # tick marks
	axis(whichaxis, at=tickslabels$decseq, labels=F, tcl=-0.3) # little tick marks
	axis(whichaxis, at=tickslabels$ticks, labels=tickslabels$labels, line=-.3, lwd=0) # labels
	mtext(axtitle, side=whichaxis, line=1.9, cex=1)
		
}
