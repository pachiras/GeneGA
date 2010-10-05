\name{show-methods}
\docType{methods}
\alias{show-methods}
\alias{show.GeneGA}
\alias{show.GeneFoldGA}
\alias{show,GeneFoldGA-method}
\alias{show,GeneGA-method}
\title{ show methods for displaying the results of GeneGA and GeneFoldGA }
\description{
Displaying the results of GeneGA and GeneFoldGA, the first three distinctive and optimum sequences can be returned, as well as their overall evaluation values, CAI values and minimum free energys.
}
\section{Methods}{
\describe{
\item{"GeneFoldGA"}{ Show method of GeneFoldGA }
\item{"GeneGA"}{ Show method of GeneGA }
}
}
\author{
Zhenpeng Li, Fei Li, Xiaochen Bo and Shengqi Wang
}
\examples{
seqfile=system.file("sequence","EGFP.fasta",package="GeneGA")
seq=unlist(getSequence(read.fasta(seqfile),as.string=TRUE))
GeneFoldGA.result=GeneFoldGA(sequence=seq,popSize=40,iters=150,crossoverRate=0.2,
	mutationChance=0.05,region=c(1,42))
show(GeneFoldGA.result)
}
\seealso{\code{\link{GeneGA-class}}, \code{\link{GeneFoldGA-class}}, \code{\link{GeneGA}},
\code{\link{GeneFoldGA}}, \code{\link{plotGeneGA-methods}}
}
\keyword{methods}
