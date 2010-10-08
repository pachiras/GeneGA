\name{GeneFoldGA-class}
\docType{class}
\alias{GeneFoldGA-class}
\alias{show,GeneFoldGA-method}
\title{Class "GeneFoldGA"}
\description{GeneFoldGA class for representing the GeneFoldGA results}
\section{Slots}{
  \describe{
    \item{\code{seq}:}{Object of class \code{"character"}. The mRNA sequence to optimize.}
    \item{\code{iters}:}{Object of class \code{"integer"}. The number of iteration.}
    \item{\code{popSize}:}{Object of class \code{"numeric"}. The population size.}
    \item{\code{crossoverRate}:}{Object of class \code{"numeric"}. The crossover rate of each generation. By default 0.2. }
    \item{\code{mutationChance}:}{Object of class \code{"numeric"}. Mutation chance of the gene in the unit of codon. }
    \item{\code{region}:}{Object of class \code{"ANY"}. The region of sequence to optimize, the other part of sequence outside the region is  optimized by only considering codon usage bias.}
    \item{\code{organism}:}{Object of class \code{"character"}. The organism  that gene due to express,the package contains codon usage bias information of almost 200 genomes.}
    \item{\code{eval_value}:}{Object of class \code{"numeric"}. The evaluation function values of all final population.}
    \item{\code{free_en}:}{Object of class \code{"numeric"}. The minimux free energy values of all final population.}
    \item{\code{eval_value_set}:}{Object of class \code{"numeric"}. The mean of evaluation function values of the population. }
    \item{\code{eval_value_set02}:}{Object of class \code{"numeric"}. The maximum of evaluation function values of the population. }
    \item{\code{population}:}{Object of class \code{"character"}. The final population undergone the genetic algorithm. }
    \item{\code{ramp}:}{Object of class \code{"ANY"}. Ramp specifies the region with low translation efficiency. Generally,the first 90 to 150 bases are the ramp region of gene,which are deemed to be an evolutionarily conserved mechanism for controlling the efficiency of protein translation. Referring to the reference for more detailed description. }
  }
}

\section{Methods}{
  \describe{
    \item{show}{\code{signature(object = "GeneFoldGA")}: Displaying the results of GeneFoldGA, the first three distinctive and optimum sequences can be returned, as well as the corresponding minimum free energys. }
    \item{plotGeneGA}{\code{signature(object = "GeneFoldGA")}: Visualizing the variation of maximum and mean variable values during the progress that genetic algorithm performed. Furthermore, the plot also can be used to check the results whether well converged. }
  }
}
\seealso{\code{\link{GeneFoldGA}}, \code{\link{GeneGA}}, 
\code{\link{show-methods}}, \code{\link{plotGeneGA-methods}}
}
\keyword{classes}

