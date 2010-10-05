\name{wSet}
\alias{wSet}
\docType{data}
\title{
Codon Adaptation Index (CAI) w tables
}
\description{
Codon Adaptation Index (CAI) w tables of almost 200 genomes,the computation is based on the highly expressed genes of these genomes.
}
\usage{data(wSet)}
\format{
  A data frame with 200 observations on the following 64 variables, which denote the 64 codes.
}
\details{
Codons can be accessed by row.names(wSet) and the genome names can be accessed by names(wSet).
}
\source{
The data were combined from two sources: the first part is the first three columns,which contains w tables masked from the caitab dataset of seqinr package,while the other 197 columns were computed based the predicted highly expressed genes of 197 genomes from HEG-DB database.
}
\references{
Puigbo P, Romeu A, Garcia-Vallve S. HEG-DB: a database of predicted highly expressed genes in prokaryotic complete genomes under translational selection. Nucl. Acids Res. 2008, 36:D524-527.

Charif D, Lobry J. 2007. SeqinR 1.0-2. A Contributed Package to the R Project for Statistical Computing Devoted to Biological Sequences Retrieval and Analysis. In: Structural Approaches to Sequence Evolution (Bastolla U, Porto M, Roman E, Vendruscolo M, eds.), Berlin Heidelberg: Springer, p207-232.
}
\examples{
data(wSet)
}
\seealso{ \code{\link{GeneCodon}}, \code{\link{GeneGA}}, \code{\link{GeneFoldGA}}
}
\keyword{datasets}
