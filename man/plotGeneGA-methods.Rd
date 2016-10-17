\name{plotGeneGA-methods}
\docType{methods}
\alias{plotGeneGA}
\alias{plotGeneGA-methods}
\alias{plotGeneGA.GeneGA}
\alias{plotGeneGA.GeneFoldGA}
\alias{plotGeneGA,GeneFoldGA-method}
\alias{plotGeneGA,GeneGA-method}
\title{ plotGeneGA methods of GeneGA and GeneFoldGA objects }
\description{
 plotGeneGA implement plotGeneGA methods for GeneGA and GeneFoldGA objects respectively.The functions visualize the variation of optimized and mean overall evaluation values and variable values during the progress that genetic algorithm performed. Furthermore, the plot also can be used to check the results whether well converged.
}
\section{Methods}{
\describe{

\item{x = "GeneFoldGA"}{ plotGeneGA method for GeneFoldGA}

\item{x = "GeneGA"}{ plotGeneGA method for GeneGA, it also has a parameter "type", which can be assigned with one of the three values, the default value is "default", which show the variation of  mean  and optimized overall evaluation value along with the generation. while 1, 2 display the variation of mean  and optimized overall CAI value and minimum free energy respectively, and 3 display the scatter plot of two variables--CAI and minimum free energy }
}}

\author{
Zhenpeng Li, Fei Li, Xiaochen Bo and Shengqi Wang
}

\examples{
seqfile=system.file("sequence", "EGFP.fasta", package="GeneGA")
seq=unlist(getSequence(read.fasta(seqfile), as.string=TRUE))
GeneGA.result=GeneFoldGA(sequence=seq, popSize=40, iters=100, crossoverRate=0.3, 
         mutationChance=0.05, region=c(1,42))
plotGeneGA(GeneGA.result)
}
\seealso{\code{\link{GeneGA-class}}, \code{\link{GeneFoldGA-class}}, \code{\link{GeneGA}},
\code{\link{GeneFoldGA}}, \code{\link{show-methods}}
}
\keyword{methods}
