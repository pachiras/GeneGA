evaluationFunction <-
function(xx,w = NULL,frontValue=NULL){
	if(is.null(frontValue)){
		free_energy=fold(paste(xx,collapse=""))[[2]]
	} else {
		free_energy=fold(paste(c(s2c(frontValue),xx),collapse=""))[[2]]
	}
	CAI=cai(s2c(paste(xx,collapse="")),w)
	return(c(free_energy,CAI))
}

evaluationFoldFunction <-
function(xx,frontValue=NULL){
	if(is.null(frontValue)){
		free_energy=fold(paste(xx,collapse=""))[[2]]
	} else {
		free_energy=fold(paste(c(s2c(frontValue),xx),collapse=""))[[2]]
	}
	return(free_energy)
}

evaluationFunction_internal <-
function(xx,w = NULL,frontValue=NULL,region_=NULL,ramp_value_=NULL){
	if(is.null(frontValue)){
		free_energy=fold(paste(xx,collapse=""))[[2]]
	} else {
		free_energy=fold(paste(c(s2c(frontValue),xx),collapse=""))[[2]]
	}
	CAI=cai(s2c(paste(xx,collapse=""))[1:(ramp_value_-region_[1]+1)],w)
	CAI_=cai(s2c(paste(xx,collapse=""))[(ramp_value_-region_[1]+2):(region_[2]-region_[1]+1)],w)
	return(c(free_energy,CAI,CAI_))
}


