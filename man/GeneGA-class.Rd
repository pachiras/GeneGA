\name{GeneGA-class}
\docType{class}
\alias{GeneGA-class}
\title{Class "GeneGA"}
\description{GeneGA class for representing the GeneGA results}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("GeneGA", ...)}.
}
\section{Slots}{
All the slots of the GeneFoldGA class, plus the following slots:  
    \describe{
    \item{\code{CAI_value}:}{Object of class \code{"numeric"}. If ramp and specified region are not intersecting, it denotes the CAI values of all final population, otherwise, it denotes the CAI values of the intersected region of all final population.}
    \item{\code{CAI_value_}:}{Object of class \code{"numeric"}. When ramp and specified region are intersecting, it denotes the CAI values of the specified region that not intersecting with the specified region of all final population.}
    \item{\code{free_en_set}:}{Object of class \code{"numeric"}. The mean of minimum free energy of the population.}
    \item{\code{CAI_value_set}:}{Object of class \code{"numeric"}. If ramp and specified region are not intersecting, it denotes the mean of CAI values of the population, otherwise, it denotes the mean CAI values of the the intersected region of the population.}
    \item{\code{CAI_value_set_}:}{Object of class \code{"numeric"}. When ramp and specified region are intersecting, it denotes the mean CAI values of the specified region that not intersecting with the specified region of the population.}
    \item{\code{free_en_set02}:}{Object of class \code{"numeric"}. The maximum of minimum free energy of the population.}
    \item{\code{CAI_value_set02}:}{Object of class \code{"numeric"}. If ramp and specified region are not intersecting, it denotes the maximum of CAI values of the population, otherwise, it denotes the maximum CAI values of the the intersected region of the population.}
    \item{\code{CAI_value_set02_}:}{Object of class \code{"numeric"}. When ramp and specified region are intersecting, it denotes the maximum CAI values of the specified region that not intersecting with the specified region of the population.}
  }
}

\section{Methods}{
  \describe{
    \item{show}{\code{signature(object = "GeneGA")}: ... }
    \item{plotGeneGA}{\code{signature(object = "GeneGA")}: ... }
  }
}
\seealso{\code{\link{GeneGA}}, \code{\link{GeneFoldGA}}, \code{\link{GeneCodon}},
\code{\link{show-methods}}, \code{\link{plotGeneGA-methods}}
}
\keyword{classes}

