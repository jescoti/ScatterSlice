# Function to determine color of boxes
ColrFind		<- function(zval, minval, maxval, flip, colorbits=32, colorpalette){
	#show(sprintf("Colorfind: minval=%.4g, maxval=%.4g, zval=%.4g", minval, maxval, zval))
	colr.rgb 	<- rgb(1,1,1)
	if(!(is.nan(minval)|is.nan(maxval)|is.nan(zval))){
		zrange	<- ifelse(minval == maxval, 1, maxval - minval)			
		colr.ratio	<- (zval-minval)/zrange
		colr.idx	<- ceiling(colr.ratio*colorbits)  
		if (zval >= maxval) { colr.idx <- colorbits }
		if (zval <= minval) { colr.idx <- 1 }			
		if (flip) {	colorpalette <- colorpalette[length(colorpalette):1] }
		colr.rgb	<- colorpalette[colr.idx]	
	}
	#show(sprintf("colr.idx=%.2g",colr.idx))
	return(colr.rgb)
}

neglog	<- function(x){
	if(x>1){ ans <- log10(x) }
	if((x >= -1) & (x <= 1)) { ans <- 0 }
	if(x < -1) { ans <- -1*log10(abs(x)) }
	return(ans)
}


Colorspec <- function (zval, minval, maxval, flipjet=F, colormap="jet", negvals.ok=T, logdisp=F) {
	#show(sprintf("Starting Colorspec: minval=%.2g, maxval=%.2g, zval=%.2g", minval, maxval, zval))
	colorbits	<- 32
	colr.rgb <- rgb(1,1,1)
	
	
	
	if(!(is.nan(minval)|is.nan(maxval)|is.nan(zval))){
		
		
		## Color palettes
		jet.colors 	<- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan",
    	                 "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
    	Jet <- jet.colors(colorbits)
    	neg.colors.f	<- colorRampPalette(c("black", "blue", "deepskyblue", "aquamarine")) 
		negcolors		<- neg.colors.f(33)[2:33]
		pos.colors.f	<- colorRampPalette(c("black", "yellow", "red", "indianred2"))
		poscolors		<- pos.colors.f(33)[2:33]
		heatcolors		<- heat.colors(colorbits)
		colorpalette	<- switch(colormap,
									jet = Jet, 
									blues = negcolors, 
									yellows = poscolors,
									heat = heatcolors
									)
		# These are here just in case we want to be able to change the color scale later
		
		if(negvals.ok==F | (sign(minval)>=0 & sign(maxval)==1)){
			## if all values are positive, or if neg values will be clipped.
			#	show("Positive only")
			if(minval == 0) { minval <- 1 } # minval will only possibly be negative for amplitude
			if(logdisp){
				minval	<- log10(minval)
				maxval	<- log10(maxval)
				zval	<- log10(zval)
			}
			colr.rgb	<- ColrFind(zval, minval, maxval, flipjet, colorbits, colorpalette)
		}else{
			if((sign(minval)==-1 & sign(maxval)==-1)){
				# If both are negative, still use the "jet" colormap.
				flip	<- T
				if(logdisp){
					## If the output should be displayed as log
					minval	<- neglog(minval)
					maxval	<- neglog(maxval)
					zval	<- neglog(zval)
				}
				colr.rgb	<- ColrFind(zval, minval, maxval, flip, colorbits, colorpalette)
				
				
			}else{
				## If bottom and top have different signs, use a +- color scheme
				disprange	<- c(minval, maxval)
				maxabv		<- max(abs(c(minval, maxval)))	 ## find which side (- or +) has greater magnitude.
				hasmaxabv	<- which(abs(disprange)==maxabv)
				
				if(logdisp){
					maxval	<- neglog(maxval)
					minval	<- neglog(minval)
					zval	<- neglog(zval)
				}
				
				fulllist	<- seq(minval, maxval, length.out=32)
				zero.idx	<- which.min(abs(fulllist)) # Closest to zero
				fulllist[zero.idx]	<- 0				# Just call it zero
				col.idx	<- which.min(abs(fulllist-zval))		# Find where the zval sits in the list
				
				#show("full list")
				#show(fulllist)
				show(sprintf("zval= %.3g, col.idx= %d", zval, col.idx))
				
				
				if(col.idx == zero.idx){
					colr.rgb	<- rgb(0,0,0)
					# If the number is close to zero, color the square black
				}
				if(col.idx > zero.idx){
					colr.rgb	<- poscolors[col.idx-zero.idx]
					# If the number is greater than zero, use the poscolors scale, start at 1
				}
				if(col.idx < zero.idx){
					colr.rgb	<- negcolors[zero.idx - col.idx]
					# If the number is less that zero, use the negcolors scale
				}
			}
		}
	}
	return(colr.rgb)
}
