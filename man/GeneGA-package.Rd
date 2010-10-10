\name{GeneGA-package}
\alias{GeneGA-package}
\docType{package}
\title{
Designing gene based on mRNA secondary structure and codon usage bias using Genetic algorithm, GeneGA includes the information of highly expressed genes of almost 200 genomes.
}
\description{
R based Genetic algorithm for gene expression optimization considering mRNA secondary structure and codon usage bias
}
\details{
\tabular{ll}{
Package: \tab GeneGA\cr
Type: \tab Package\cr
Version: \tab 0.99.1\cr
Date: \tab 2010-10-10\cr
License: \tab GPL version 2\cr
LazyLoad: \tab yes\cr
}
}
\author{
Zhenpeng Li, Fei Li, Xiaochen Bo and Shengqi Wang

Maintainer: Zhenpeng Li <zpli21@gmail.com>
}
\references{
Liu L.,Kang L.S.,Chen Y.P.. (1993). Non-numerical parallel algorithms(The second volume)--genetic algorithms,Science Press(In Chinese)

Tuller T, Carmi A, Vestsigian K, et al. An Evolutionarily Conserved Mechanism for Controlling the Efficiency of Protein Translation. Cell 2010, 141:344-354
}
\keyword{ package }
\seealso{\code{\link{GeneGA}}, \code{\link{GeneFoldGA}}, \code{\link{GeneCodon}}, \code{\link{GeneGA-class}}, \code{\link{GeneFoldGA-class}},
\code{\link{show-methods}}, \code{\link{plotGeneGA-methods}}, \code{\link{wSet}}
}
\examples{
seqfile=system.file("sequence","EGFP.fasta",package="GeneGA")
seq=unlist(getSequence(read.fasta(seqfile),as.string=TRUE))
GeneGA.result=GeneFoldGA(sequence=seq,popSize=40,iters=100,crossoverRate=0.2,
	mutationChance=0.05,region=c(1,42))
plotGeneGA(GeneGA.result)
show(GeneGA.result)
}
