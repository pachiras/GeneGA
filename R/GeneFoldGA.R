GeneFoldGA <-
function(sequence=NULL,
		popSize=50,
		iters=100,
		crossoverRate=0.2,
		mutationChance=0.05,
		region=NULL,
		showGeneration=TRUE,
		frontSeq=NULL,
		organism="ec",
		ramp=FALSE,
		numcode=1
		){
	if(sum(s2c(toupper(sequence))%in%c("A","T","G","C")) != nchar(sequence)){
		stop("The input sequence must be nucleotide sequence only containg 'ATGC'")
	}
	if(is.null(region)) region=c(1,nchar(sequence))
	assign("translate",seqinr::translate,envir=.GlobalEnv)
	#construct a hash between amino acids and codons
	code_invert=hash()
	stardardSeq="GAGTCGTTCAGTAAACACTGTCAACGGCAGGCTCGAATATAACCCGAAACAACGTGGAATTTTCTCACCGTCCTTCCGTGCTCTCTGGGGGTTATGTAGCCAGCAACTGTAAGATACGGCAGCAGGAACTGAAAGTCAGTGATCGCGCGTGGATTACCTCATTATGACCTATTGCGCGATGCCGGTTCCATT"
	for(i in seq(1,192,3)){
		code_invert[[substr(stardardSeq,i,i+2)]]=translate(s2c(stardardSeq)[i:(i+2)],numcode = numcode)
	}
	code=invert(code_invert)
	#if not full sequence was specified, subtract the target sequence from given sequence
	seq=substr(sequence,region[1],region[2])
	if((region[2]-region[1]+1)%%3 != 0){
		stop("the given region must be the multiple of  3, please reassign the region")
		}
	amino_seq=translate(s2c(seq),numcode = numcode)
	#producing random population
	codonsize=nchar(seq)/3
	population=matrix(nrow=popSize,ncol=codonsize)
	for(i in 1:dim(population)[1]){
		for(j in 1:dim(population)[2]){
			population[i,j]=sample(code[[amino_seq[j]]],1)
		}
	}
	
	eval_value=rep(NA,popSize)
	free_en=rep(NA,popSize)

	#do iters
	eval_value_set=c()
	free_en_set=c()
	eval_value_set02=c()
	free_en_set02=c()
	for(iter_order in 1:iters){
		if(showGeneration){
			cat("The current generation is:",iter_order,"\n")
		}
		#re-evaluation
		for(i in 1:popSize){
			if(is.na(free_en[i])){
				free_en[i]=evaluationFoldFunction(population[i,],frontValue=frontSeq)
				}
		}
		#compute evaluation values
		eval_value=sapply(free_en,exp)
		#cat(eval_value,"\n",file="result.txt",sep=" ",append=TRUE)
		#store the mean eval value and maxium eval value into eval_value_set and eval_value_set02
		eval_value_set=append(eval_value_set,mean(free_en))
		eval_value_set02=append(eval_value_set02,max(free_en))
		#do iteration when iter is less than generation
		if(iter_order < iters){
			#selection process
			#firstly, select the population based on the integer part 
			#of the product of proportion of its evaluation value and popSize.
			new_pop=matrix(nrow=popSize,ncol=codonsize)
			eval_prop=rep(NA,popSize)
			eval_digit_prop=rep(NA,popSize)
			eval_accum=rep(NA,popSize)
			eval_value_new=rep(NA,popSize)
			free_en_new=rep(NA,popSize)
			eval_sum=sum(eval_value)
			for(i in 1:popSize){
				eval_prop[i]=eval_value[i]/eval_sum
			}
			num_each=floor(eval_prop*popSize)
			n=1
			for(num in 1:popSize){
				if(num_each[num] > 1){
					new_pop[n:(n+num_each[num]-1),]=matrix(rep(population[num,],
						times=num_each[num]),nrow=num_each[num],byrow=TRUE)
					free_en_new[n:(n+num_each[num]-1)]=rep(free_en[num],times=num_each[num])
					n=n+num_each[num]
				}
				else if(num_each[num] == 1){
					new_pop[n,]=population[num,]
					free_en_new[n]=free_en[num]
					n=n+1
			}
			}
			#secondly, select the population by the digit part of eval_prop*popSize 
			#using roulette algorithm
			eval_digit=sapply(eval_prop*popSize,function(x)(x-floor(x)))
			eval_sum_digit=sum(eval_digit)
			for(i in 1:popSize){
				eval_digit_prop[i]=eval_digit[i]/eval_sum_digit
			}
			for(i in 1:popSize){
				eval_accum[i]=sum(eval_digit_prop[1:i])
			}
			
			for(i in (sum(num_each)+1):popSize){
				random_prop=runif(1)
				for(j in 1:popSize){
					if(eval_accum[j] > random_prop){
						new_pop[i,]=population[j,]
						free_en_new[i]=free_en[j]
					 	break
					}
				}
			}
			
			# crossover
			#preserve the first ten biggest evaluation value to prevent their changing
			new_eval_value=sapply(free_en_new,exp)
			eval_value_index=sort(new_eval_value,decreasing=TRUE,index=TRUE)$ix
			selected=eval_value_index[10:popSize]
			sample_pop=sample(selected,round(crossoverRate*popSize))
			free_en_new[sample_pop]=NA
			i=1
			while(i < length(sample_pop)){
				crossOverPoint=sample(1:(codonsize-1),1)
				tt=new_pop[sample_pop[i],]
				new_pop[sample_pop[i],]=c(new_pop[sample_pop[i],1:crossOverPoint],
					new_pop[sample_pop[i+1],(crossOverPoint+1):codonsize])
				new_pop[sample_pop[i+1],]=c(new_pop[sample_pop[i+1],
					1:crossOverPoint],tt[(crossOverPoint+1):codonsize])
				i=i+2
			}

			#mutation
			#reassign the evaluation value
			for(i in selected){
				for(j in 1:codonsize){
					if(runif(1) <= mutationChance){
						new_pop[i,j]=sample(code[[code_invert[[new_pop[i,j]]]]],1);
						free_en_new[i]=NA;
					}
				}
			}
			population=new_pop
			free_en=free_en_new
		}
	}
	#report GA results
		results <- new("GeneFoldGA", seq=sequence,iters=iters,popSize=popSize,
				crossoverRate=crossoverRate,mutationChance=mutationChance,region=region,
				eval_value=eval_value,free_en=free_en,eval_value_set=eval_value_set,
				eval_value_set02=eval_value_set02,population=population,ramp=FALSE,organism=organism);
	return(results)
}

