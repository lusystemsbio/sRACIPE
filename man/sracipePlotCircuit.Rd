% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AllGenerics.R, R/methods.R
\name{sracipePlotCircuit}
\alias{sracipePlotCircuit}
\alias{sracipePlotCircuit,RacipeSE-method}
\title{Plot Gene Regulatory Circuit}
\usage{
sracipePlotCircuit(
  .object,
  plotToFile = FALSE,
  physics = TRUE,
  namedNodes = TRUE
)

\S4method{sracipePlotCircuit}{RacipeSE}(
  .object,
  plotToFile = TRUE,
  physics = TRUE,
  namedNodes = TRUE
)
}
\arguments{
\item{.object}{RacipeSE object
A list returned by \code{\link{sracipeSimulate}} function}

\item{plotToFile}{(optional) logical. Default \code{FALSE}. Whether to save
plots to a file.}

\item{physics}{(optional) logical. Default \code{TRUE}. Whether or not to
enable physics in the nodes of the visNetwork graph.}

\item{namedNodes}{(optional) logical. Default \code{TRUE}. Whether or not to
display gene names in the circuit visualization}
}
\value{
circuit plot
}
\description{
Plot Gene Regulatory Circuit to a file or output device using
visNetwork. Edge color coding: Transcription-"black", Protein
Degradation-"red", Signaling Interaction-"green".
}
\section{Related Functions}{


\code{\link{sracipeSimulate}},  \code{\link{sracipeKnockDown}},
\code{\link{sracipeOverExp}},  \code{\link{sracipePlotData}}
}

\examples{
data("demoCircuit")
\dontrun{
rSet <- sRACIPE::sracipeSimulate(circuit = demoCircuit, numModels = 20,
integrateStepSize = 0.1, simulationTime = 30)
sracipePlotCircuit(rSet, plotToFile = FALSE, physics = TRUE)
rm(rSet)
}
}
