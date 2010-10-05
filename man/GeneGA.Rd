\name{GeneGA}
\alias{GeneGA}
\title{
the function optimizes gene expression based on both mRNA secondary structure and codon usage bias
}
\description{
Optimizing gene expression considering both mRNA secondary structure and codon usage bias using genetic algorithm. The function has the default evaluation function with variables CAI and minimum free energy, the sum of squares of ranks of which is used as an evaluation value to pilot the evoluation. The optimum is the gene for which the evaluation value is maximum.
Results can be visualized with plotGeneGA and displayed with show.
}
\usage{
GeneGA(sequence = NULL, popSize = 50, iters = 150, crossoverRate = 0.2, 
	mutationChance = 0.05, region = NULL, organism = "ec", showGeneration = TRUE,
	frontSeq = NULL, ramp=FALSE, numcode=1)
}
\arguments{
  \item{sequence}{
the mRNA sequence to optimize
}
  \item{popSize}{
the population size
}
  \item{iters}{
the number of iteration
}
  \item{crossoverRate}{
the crossover rate of each generation. By default 0.2
}
  \item{mutationChance}{
mutation chance of the gene in the unit of codon
}
  \item{region}{
the region of sequence to optimize, the other part of sequence outside the region is  optimized by only considering codon usage bias.
}
  \item{organism}{
the organism  that gene due to express, the package contains codon usage bias information of almost 200 genomes.
}
  \item{showGeneration}{
show the generation the genetic algorithm progressing, the default value is TRUE.
}
  \item{frontSeq}{
frontSeq denotes the regulatory segment before the start codon.If frontSeq is specified, frontSeq will be considered when compute the minimum free energy.The default value is NULL.
}
  \item{ramp}{
ramp specifies the region with low translation efficiency.Generally, the first 90 to 150 bases are the ramp region of gene, which are deemed to be an evolutionarily conserved mechanism for controlling the efficiency of protein translation.Referring the reference for more detailed description.
}
  \item{numcode}{
The ncbi genetic code number for translation. By default the standard genetic code is used. Referring to the help page of "translate" function in seqinr package for details.
}
}
\value{
A GeneGA instance is returned. 
}
\references{
Liu L.,Kang L.S.,Chen Y.P.. (1993). Non-numerical parallel algorithms(The second volume)--genetic algorithms,Science Press(In Chinese).

Tuller T, Carmi A, Vestsigian K, et al. An Evolutionarily Conserved Mechanism for Controlling the Efficiency of Protein Translation. Cell 2010, 141:344-354.

}
\author{
Zhenpeng Li, Fei Li, Xiaochen Bo and Shengqi Wang
}
\seealso{\code{\link{GeneGA-class}}, \code{\link{GeneFoldGA}}, \code{\link{GeneCodon}},
\code{\link{show-methods}}, \code{\link{plotGeneGA-methods}}
}
\examples{
seqfile=system.file("sequence","EGFP.fasta",package="GeneGA")
seq=unlist(getSequence(read.fasta(seqfile),as.string=TRUE))
GeneGA.result=GeneGA(sequence=seq,popSize=40,iters=100,crossoverRate=0.4,
	mutationChance=0.1,region=c(1,60))
}
