StorePrefs	<- function(homedir="~", saveloc="~", filename="Preferences.Scs"){
	Preferences	<- list()
	Preferences$home.directory	<- homedir
	save(Preferences, file=paste(saveloc, "/", filename, sep=""))

}

SetPrefs	<- function(prefFile="~/Preferences.Scs"){
	load(prefFile)
}

GetPrefs	<- function(){
	prefFile	<- tkgetOpenFile(multiple=F, filetypes="{{Scs Pref Files} {.Scs}}", initialdir="~")
	SetPrefs(prefFile)
}


