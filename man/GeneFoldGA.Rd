\name{GeneFoldGA}
\alias{GeneFoldGA}
\title{
the function optimizes gene merely considering mRNA secondary structure.
}
\description{
Optimizing gene by only considering mRNA secondary structure using genetic algorithm. The function has the default evaluation function with minimum free energy as variable. The optimum is the gene for which the minimum free energy is maximum. Results can be visualized with plotGeneGA and displayed with show.
}
\usage{
GeneFoldGA(sequence = NULL, popSize = 50, iters = 100, crossoverRate = 0.2, 
         mutationChance = 0.05, region = NULL, showGeneration = TRUE, 
         frontSeq = NULL, organism="ec", ramp = FALSE,numcode=1)
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
mutation chance of the gene in the unit of codons
}
  \item{region}{
the region of sequence to optimize, the other part of sequence outside the region is  optimized only considering codon usage bias as well.
}
  \item{showGeneration}{
show the generation the genetic algorithm progressing, the default value is TRUE
}
  \item{frontSeq}{
frontSeq denotes the regulatory segment before the start codon. If frontSeq is specified, frontSeq will be considered when compute the minimum free energy. The default value is NULL.
}
  \item{organism}{
the organism  that gene due to express, the package contains codon usage bias information of almost 200 genomes.
}
  \item{ramp}{
the organism  that gene due to express, the package contains codon usage bias information of almost 200 genomes.
}
  \item{numcode}{
The ncbi genetic code number for translation. By default the standard genetic code is used. Referring to the help page of "translate" function in seqinr package for details.
}
}
\value{
  A GeneFoldGA instance is returned.
}
\references{
Liu L.,Kang L.S., Chen Y.P. (1993)Non-numerical parallel algorithms(The second volume)--genetic algorithms,Science Press(In Chinese)
}
\author{
Zhenpeng Li, Fei Li, Xiaochen Bo and Shengqi Wang
}

\seealso{\code{\link{GeneFoldGA-class}}, \code{\link{GeneGA}}, 
\code{\link{show-methods}}, \code{\link{plotGeneGA-methods}}
}

\examples{
seqfile=system.file("sequence","EGFP.fasta",package="GeneGA")
seq=unlist(getSequence(read.fasta(seqfile),as.string=TRUE))
GeneGA.result=GeneFoldGA(sequence=seq,popSize=40,iters=100,crossoverRate=0.3,
         mutationChance=0.05,region=c(1,42))
}
