#!/bin/bash
	
	R --no-save <<EOF
	install.packages("tkrplot", repos="http://cran.r-project.org")
	install.packages("Hmisc", repos="http://cran.r-project.org")
	install.packages("foreach", repos="http://cran.r-project.org")
	install.packages("doMC", repos="http://cran.r-project.org")
	install.packages("gplots", repos="http://cran.r-project.org")
	q()
	EOF
	
