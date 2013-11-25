Colorbar <- function(bottom=0, top=1, titlestring="", islog=T, resetpar=T, negvals.ok=F,
						colorscale="jet", flip=F) {
	
	show("Creating Colorbar")						
	if(resetpar){
		oldpar <- par(no.readonly = T)
	}
	par(bg="white")
	parplt3 <- par("plt")[3]
	parplt4 <- par("plt")[4]
	par("plt"=c(0,0.95,parplt3, parplt4))
	plot(c(0,1),c(0,1), pch="", xlab=NA, ylab=NA, xaxt="n", yaxt="n")
	
	show(sprintf("bottom, top values for colorbar: %.2g and %.2g", bottom, top))
	
	colorvec	<- vector("numeric", 32)
	for(c.idx in 1:32){
		colorvec[c.idx]	<- Colorspec(c.idx, 1, 32, logdisp=F)
	}
	
	if(sign(bottom)>=0 & sign(top)==1){
		#show("Both ends of colorbar are non-negative")
		nlabels	<- 5
		cb.labels.ex	<- vector("list", nlabels)
		## If the scale is log, make the labels evenly spaced in log
		## and display as 10^x
		if(islog){
			if(bottom==0){bottom=1}
			colorbar.labels <- round(seq(log10(bottom), log10(top), length.out=nlabels),1)
			for(idx in 1:nlabels){
				cb.labels.ex[[idx]]	<- as.expression(bquote(.(10)^.(colorbar.labels[idx])))
			}
		}else{
			colorbar.labels <- signif(seq(bottom, top, length.out=nlabels),2)
			for(idx in 1:nlabels){
				cb.labels.ex[[idx]]	<- as.expression(colorbar.labels[idx])
			}
		}
		
	}else{
		if(sign(bottom)==-1 & sign(top)==-1){
			#show("Both ends of colorbar are negative")
			nlabels	<- 5
			cb.labels.ex	<- vector("list", nlabels)
			flip	<- T
			## If the scale is log, make the labels evenly spaced in log
			## and display as 10^x
			
			if(islog){
				colorbar.labels <- round(seq(-1*neglog(bottom), -1*neglog(top), length.out=nlabels),1)
				for(idx in 1:nlabels){
					cb.labels.ex[[idx]]	<- as.expression(bquote(.(-10)^.(colorbar.labels[idx])))
				}
			}else{
				colorbar.labels <- signif(seq(bottom, top, length.out=nlabels),2)
				for(idx in 1:nlabels){
					cb.labels.ex[[idx]]	<- as.expression(colorbar.labels[idx])
				}
			}
			
			
			
			
		}else{
			## If bottom and top have different signs, use a +- color scheme
			#show("Colorbar from neg to pos")
			
			colorscale	<- "negpos"
			nlabels	<- 7
			cb.labels.ex	<- vector("list", nlabels)
			
			if(islog==T){
				bottom	<- neglog(bottom)
				top		<- neglog(top)
			}
			show(c(bottom, top))
			idx.vec	<- vector("numeric", 32)
			colorvalseq	<- seq(bottom, top, length.out=32)
			zero.idx	<- which.min(abs(colorvalseq))
			# If zero index is too close to the top or the bottom,
			re.start	<- F
			if(zero.idx <= 2){
				bottom	<- 1
				re.start<- T
			}
			if(zero.idx >= 30){
				top	<- -1
				re.start<-T
			}
			# Go back and use a regular color scheme
			if(re.start==T){
				Colorbar(bottom, top, titlestring, islog, resetpar, negvals.ok,	colorscale, flip)
			}
			
			# Otherwise, put in a neg to pos colorbar
			colorvalseq[zero.idx]	<- 0
			idx.vec[1:(zero.idx-1)] 	<- -1*(((zero.idx)-1):1)
			idx.vec[(zero.idx+1):32] 	<- 1:(32-zero.idx)
			show(idx.vec)
			for (idx in 1:32){
				colorvec[idx]	<- Colorspec(idx.vec[idx], idx.vec[1], idx.vec[32])
			}
			colorvec[zero.idx]	<- rgb(0,0,0) #black
			colorbar.labels <- signif(seq(bottom, top, length.out=nlabels),2)
			for(idx in 1:nlabels){
				cb.labels.ex[[idx]]	<- as.expression(colorbar.labels[idx])
			}
		}
			
		
	}
	#show(sprintf("bottom, top values for colorbar: %.2g and %.2g", bottom, top))	
	
	
	## Flipping keeps colors, but flips locations. All flipping is done here
	if(flip){
		label.Y.locations	<- seq(0.87, 0.03, length.out=nlabels)
	}else{
		label.Y.locations	<- seq(0.03, 0.87, length.out=nlabels)
	}
	cb.y.locations		<- 0.9*(0:32)/32
	
	## Colored rectangles
	for(idx in 1:32){
		#show(sprintf("point = %.4g", point))
		color	<- colorvec[idx]
		rect(.05, cb.y.locations[idx], .35, cb.y.locations[idx+1], col=color, lwd=0.5)
	}
	
	## Place numeric labels
	for(idx in 1:nlabels){
		text(x=0.4, y= label.Y.locations[idx], labels = cb.labels.ex[[idx]], adj=0)
	}
	
	## Title Colorbar
	text(x=0,y=.95, pos=4, labels=titlestring)
	
	
	
	## Set parameters back to orginal
	if(resetpar){
		par(oldpar) #(reset parameters)
	}
	show("Colorbar finished")
	
}
