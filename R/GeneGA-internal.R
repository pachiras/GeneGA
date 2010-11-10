setClass("GeneFoldGA",
	representation(
		seq="character",
		iters="numeric",
		popSize="numeric",
		crossoverRate="numeric",
		mutationChance="numeric",
		region="ANY",
		organism="character",
		eval_value="numeric",
		free_en="numeric",
		eval_value_set="numeric",
		eval_value_set02="numeric",
		population="matrix",
		ramp="ANY"
	)
)

setClass("GeneGA",
	representation(
		CAI_value="numeric",
		CAI_value_="numeric",
		free_en_set="numeric",
		CAI_value_set="numeric",
		CAI_value_set_="numeric",
		free_en_set02="numeric",
		CAI_value_set02="numeric",
		CAI_value_set02_="numeric"
	),
	contains="GeneFoldGA")

setGeneric("plotGeneGA", 
	function (x, ...)
	standardGeneric("plotGeneGA"))

setMethod("plotGeneGA",signature(x="GeneGA"),
	function(x,type="default"){
		if(type == "default"){
			range = range(x@eval_value_set,x@eval_value_set02)
			plot(x@eval_value_set,ylim=range,xlab=NA,ylab=NA,col=12,pch=20)
			par(new=TRUE)
			plot(x@eval_value_set02,ylim=range,xlab="Generation",ylab="Evaluation Value",col=2,pch=20)
		} else if(type == 1) {
			if(x@ramp==FALSE | x@ramp >= x@region[2] | x@region[1] > x@ramp ){
				range = range(x@CAI_value_set,x@CAI_value_set02)
				plot(x@CAI_value_set,ylim=range,xlab=NA,ylab=NA,col=12,pch=20)
				par(new=TRUE)
				plot(x@CAI_value_set02,ylim=range,xlab="Generation",ylab="CAI Value",col=2,pch=20)
			} else {
				range = range(x@CAI_value_set,x@CAI_value_set02,x@CAI_value_set_,x@CAI_value_set02_)
				plot(x@CAI_value_set,ylim=range,xlab=NA,ylab=NA,col=12,pch=20)
				par(new=TRUE)
				plot(x@CAI_value_set02,ylim=range,xlab=NA,ylab=NA,col=2,pch=20)
				par(new=TRUE)
				plot(x@CAI_value_set_,ylim=range,xlab=NA,ylab=NA,col=12,pch=20)
				par(new=TRUE)
				plot(x@CAI_value_set02_,ylim=range,xlab="Generation",ylab="CAI Value",col=2,pch=20)
				}
				
		} else if(type == 2) {
			range = range(x@free_en_set,x@free_en_set02)
			plot(x@free_en_set,ylim=range,xlab=NA,ylab=NA,col=12,pch=20)
			par(new=TRUE)
			plot(x@free_en_set02,ylim=range,xlab="Generation",ylab="Free Energy",col=2,pch=20)
		} else if(type == 3) {
			if(x@ramp==FALSE | x@ramp >= x@region[2] | x@region[1] > x@ramp ){
				xRange=range(x@CAI_value_set,x@CAI_value_set02)
				yRange=range(x@free_en_set,x@free_en_set02)
				plot(x@CAI_value_set,x@free_en_set,xlim=xRange,ylim=yRange,col=12,xlab=NA,ylab=NA,pch=20)
				par(new=TRUE)
				plot(x@CAI_value_set02,x@free_en_set02,xlim=xRange,ylim=yRange,
					col=2,xlab="CAI Value",ylab="Free Energy",pch=20)
			} else {
				stop("Type 3 is not available while region and ramp are intersect")
				}
		} else {
			stop(paste("Plot type",type,"is not a supported type",sep=" "))
	}
}
)

setMethod("plotGeneGA",signature(x="GeneFoldGA"),
	function(x){
		range = range(x@eval_value_set,x@eval_value_set02)
			plot(x@eval_value_set,ylim=range,xlab=NA,ylab=NA,col=12,pch=20)
			par(new=TRUE)
			plot(x@eval_value_set02,ylim=range,xlab="Generation",ylab="Evaluation Value",col=2,pch=20)
}
)

setMethod("show","GeneGA",
	function(object){
	output=paste(
	"  GA Settings:", "\n",
	"  Population size       = ", object@popSize, "\n",
	"  Number of Generations = ", object@iters, "\n",
	"  crossoverRate         = ", object@crossoverRate, "\n",
	"  Mutation Chance       = ", object@mutationChance, "\n","\n",
        sep="")
    	cat(output)
    	seq=object@seq
   	ramp=object@ramp
    	region=object@region
	if(ramp == FALSE){
		seq_first=GeneCodon(substr(seq,1,(region[1]-1)),organism=object@organism)
		seq_last=GeneCodon(substr(seq,(region[2]+1),nchar(seq)),organism=object@organism)
	} else if(ramp!=FALSE & region[2] <= ramp){
	  	seq_first=GeneCodon(substr(seq,1,(region[1]-1)),max=F,organism=object@organism)
	  	seq_last=paste(GeneCodon(substr(seq,(region[2]+1),ramp),max=F,organism=object@organism),
			GeneCodon(substr(seq,(ramp+1),nchar(seq)),organism=object@organism),sep="")
	} else if(ramp!=FALSE & region[1] > ramp){
	  	seq_first=paste(GeneCodon(substr(seq,1,ramp),max=F,organism=object@organism),
			GeneCodon(substr(seq,(ramp+1),(region[1]-1)),organism=object@organism),sep="")
	  	seq_last=GeneCodon(substr(seq,(region[2]+1),nchar(seq)),organism=object@organism)
	} else if(ramp!=FALSE & region[1] < ramp & region[2] > ramp){
	  	seq_first=GeneCodon(substr(seq,1,(region[1]-1)),max=F,organism=object@organism)
	  	seq_last=GeneCodon(substr(seq,(region[2]+1),nchar(seq)),organism=object@organism)
	}

	if (!is.null(object@population)){
		first_three=sort(unique(object@eval_value),decreasing=TRUE)[1:3]
		for(one in first_three){
     	 		obtained=object@eval_value%in%one
     	 	 	if(sum(obtained) == 1){
     	 	 		select_solution=paste(object@population[obtained,],collapse="")
     	 	 		select_solution=paste(seq_first,select_solution,seq_last,sep="")
				select_solution=paste(sapply(0:(nchar(select_solution)%/%70),
					function(x){substr(select_solution,(x*70+1),(x+1)*70)}),collapse="\n")
	     	 	if(ramp!=FALSE & region[1] < ramp & region[2] > ramp){
	     	 		cat(" evaluaton value = ",object@eval_value[obtained],"\n",
	     	 		"free energy     = ",object@free_en[obtained],"\n",
	     	 		"the two CAI value is ",object@CAI_value[obtained],object@CAI_value_[obtained],"\n",
					select_solution,"\n")
	     	 	}else{
		     	 	cat(" evaluaton value = ",object@eval_value[obtained],"\n",
		     	 	"free energy     = ",object@free_en[obtained],"\n",
		     	 	"CAI  value      = ",object@CAI_value[obtained],"\n",
					select_solution,"\n")
					}
				} else {
					select_solution = paste(object@population[obtained,][1,],collapse="")
     	 	 			select_solution=paste(seq_first,select_solution,seq_last,sep="")
					select_solution=paste(sapply(0:(nchar(select_solution)%/%70),
						function(x){substr(select_solution,(x*70+1),(x+1)*70)}),collapse="\n")
				if(ramp!=FALSE & region[1] < ramp & region[2] > ramp){
	     	 	 		cat(" evaluaton value = ",object@eval_value[obtained][1],"\n",
	     	 	 			"free energy     = ",object@free_en[obtained][1],"\n",
	     	 	 			"the two CAI value is ",object@CAI_value[obtained][1],"\t",
						object@CAI_value_[obtained][1],"\n",
						select_solution,"\n")
					}else{
						cat(" evaluaton value = ",object@eval_value[obtained][1],"\n",
	     	 	 				"free energy     = ",object@free_en[obtained][1],"\n",
	     	 	 				"CAI  value      = ",object@CAI_value[obtained][1],"\n",
							select_solution,"\n")
     	 	 	}
     	}
     	}
		}
}
)
setMethod("show","GeneFoldGA",
	function(object){
	output=paste(
	"  GA Settings:", "\n",
	"  Population size       = ", object@popSize, "\n",
	"  Number of Generations = ", object@iters, "\n",
	"  crossoverRate         = ", object@crossoverRate, "\n",
	"  Mutation Chance       = ", object@mutationChance, "\n","\n",
        sep="")
    	cat(output)
    	seq=object@seq
   	ramp=object@ramp
    	region=object@region
	if(ramp == FALSE){
	  	seq_first=GeneCodon(substr(seq,1,(region[1]-1)),organism=object@organism)
		seq_last=GeneCodon(substr(seq,(region[2]+1),nchar(seq)),organism=object@organism)
	} else if(ramp!=FALSE & region[2] <= ramp){
	  	seq_first=GeneCodon(substr(seq,1,(region[1]-1)),max=F,organism=object@organism)
	  	seq_last=paste(GeneCodon(substr(seq,(region[2]+1),ramp),max=F,organism=object@organism),
			GeneCodon(substr(seq,(ramp+1),nchar(seq)),organism=object@organism),sep="")
	} else if(ramp!=FALSE & region[1] > ramp){
	  	seq_first=paste(GeneCodon(substr(seq,1,ramp),max=F,organism=object@organism),
			GeneCodon(substr(seq,(ramp+1),(region[1]-1)),organism=object@organism),sep="")
	  	seq_last=GeneCodon(substr(seq,(region[2]+1),nchar(seq)),organism=object@organism)
	} else if(ramp!=FALSE & region[1] < ramp & region[2] > ramp){
	  	seq_first=GeneCodon(substr(seq,1,(region[1]-1)),max=F,organism=object@organism)
	  	seq_last=GeneCodon(substr(seq,(region[2]+1),nchar(seq)),organism=object@organism)
	}

	if (!is.null(object@population)){
		first_three=sort(unique(object@eval_value),decreasing=TRUE)[1:3]
     	 	for(one in first_three){
     	 		obtained=object@eval_value%in%one
     	 	 	if(sum(obtained) == 1){
     	 	 		select_solution=paste(object@population[obtained,],collapse="")
     	 	 		select_solution=paste(seq_first,select_solution,seq_last,sep="")
				select_solution=paste(sapply(0:(nchar(select_solution)%/%70),
					function(x){substr(select_solution,(x*70+1),(x+1)*70)}),collapse="\n")
				cat( "free energy     = ",object@free_en[obtained],"\n",select_solution,"\n")
			} else {
				select_solution = paste(object@population[obtained,][1,],collapse="")
     	 	 		select_solution=paste(seq_first,select_solution,seq_last,sep="")
				select_solution=paste(sapply(0:(nchar(select_solution)%/%70),
					function(x){substr(select_solution,(x*70+1),(x+1)*70)}),collapse="\n")
				cat( "free energy     = ",object@free_en[obtained][1],"\n",select_solution,"\n")
				}
			}
	}
}
)


