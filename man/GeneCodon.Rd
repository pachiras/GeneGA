\name{GeneCodon}
\alias{GeneCodon}
\title{
Optimizing the codon usage of gene by replacing less usage codons with the most preference used codons or replacing the most preference used codons with less usage codons.
}
\description{
The function adjusts the codon usage of gene by replacing less usage codons(the most preference codons) with the most preference codons(less usage codons).
}
\usage{
GeneCodon(seq, organism = "ec", max=TRUE, scale=0.5, numcode= 1)
}
\arguments{
  \item{seq}{
the sequence to optimize
}
  \item{organism}{
the organism  that gene due to express, GeneGA contains codon usage bias information of almost 200 genomes
}
  \item{max}{
if max is TRUE, less usage codons will be replaced with the most preference codons, and vice versa.
}
  \item{scale}{
When max is FALSE, scale is used to assign the range of less used synonymous codons of each amino acid to sample. The default value is 0.5, it means each codon will be sampled from its 50\% less used synonymous codons.
}
  \item{numcode}{
The ncbi genetic code number for translation. By default the standard genetic code is used. Referring to the help page of "translate" function in seqinr package for details
}
}
\value{
This function returns the optimized sequence as string. If max is True, sequence consisting of most preference used codons will be returned, 
}
\author{
Zhenpeng Li, Fei Li, Xiaochen Bo and Shengqi Wang
}
\seealso{\code{\link{wSet}}, \code{\link{GeneGA}}
}
\examples{
seqfile=system.file("sequence","EGFP.fasta",package="GeneGA")
seq=unlist(getSequence(read.fasta(seqfile),as.string=TRUE))
GeneCodon(seq)
}
