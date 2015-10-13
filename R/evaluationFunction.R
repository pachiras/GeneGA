#RNAfold was called in linux environment
fold <-
function(x){
	fold_r=system(paste("echo ",x,"|RNAfold",sep=""),intern =T)
	nindex=grep("[\\(\\)]",s2c(fold_r[2]))
	nlength=length(nindex)
	index01=nindex[nlength]-1
	index02=nindex[nlength-1]+1
	fold_f=as.numeric(paste(s2c(fold_r[2])[index02:index01],collapse=""))
	return(fold_f)
}
evaluationFunction <-
function(xx,w = NULL,frontValue=NULL){
	if(is.null(frontValue)){
		free_energy=fold(paste(xx,collapse=""))
	} else {
		free_energy=fold(paste(c(s2c(frontValue),xx),collapse=""))
	}
	CAI=cai(s2c(paste(xx,collapse="")),w)
	return(c(free_energy,CAI))
}

evaluationFoldFunction <-
function(xx,frontValue=NULL){
	if(is.null(frontValue)){
		free_energy=fold(paste(xx,collapse=""))
	} else {
		free_energy=fold(paste(c(s2c(frontValue),xx),collapse=""))
	}
	return(free_energy)
}

evaluationFunction_internal <-
function(xx,w = NULL,frontValue=NULL,region_=NULL,ramp_value_=NULL){
	if(is.null(frontValue)){
		free_energy=fold(paste(xx,collapse=""))
	} else {
		free_energy=fold(paste(c(s2c(frontValue),xx),collapse=""))
	}
	CAI=cai(s2c(paste(xx,collapse=""))[1:(ramp_value_-region_[1]+1)],w)
	CAI_=cai(s2c(paste(xx,collapse=""))[(ramp_value_-region_[1]+2):(region_[2]-region_[1]+1)],w)
	return(c(free_energy,CAI,CAI_))
}


