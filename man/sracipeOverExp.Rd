% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AllGenerics.R, R/methods.R
\name{sracipeOverExp}
\alias{sracipeOverExp}
\alias{sracipeOverExp,RacipeSE-method}
\title{Perform in-silico over expression analysis}
\usage{
sracipeOverExp(
  .object,
  overProduction = 10,
  nClusters = 2,
  clusterOfInterest = 2,
  plotFilename = NULL,
  plotHeatmap = TRUE,
  plotBarPlot = TRUE,
  clusterCut = NULL,
  plotToFile = FALSE
)

\S4method{sracipeOverExp}{RacipeSE}(
  .object,
  overProduction = 10,
  nClusters = 2,
  clusterOfInterest = 2,
  plotFilename = NULL,
  plotHeatmap = TRUE,
  plotBarPlot = TRUE,
  clusterCut = NULL,
  plotToFile = FALSE
)
}
\arguments{
\item{.object}{RacipeSE object generated by
\code{\link{sracipeSimulate}} function.}

\item{overProduction}{(optional) Percentage to which production rate
decreases on knockdown. Uses a default value of 10 percent.}

\item{nClusters}{(optional) Number of clusters in the data. Uses a default
value of 2.}

\item{clusterOfInterest}{(optional) cluster number (integer)
to be used for arranging
the transcription factors}

\item{plotFilename}{(optional) Name of the output file.}

\item{plotHeatmap}{logical. Default TRUE. Whether to plot the heatmap or not.}

\item{plotBarPlot}{logical. Default TRUE. Whether to plot the barplot.}

\item{clusterCut}{integer or character. The cluster assignments.}

\item{plotToFile}{logical. Default FALSE.}
}
\value{
List containing fraction of models in different clusters
 in the original simulations and after knowcking down different genes.
Additionaly, it generates two pdf files in the results folder.
First is barplot
showing the percentage of different clusters in the original simulations
and after knocking down each gene. The second pdf contains the heatmap of
clusters after marking the models with cluster assignments.
}
\description{
Calculates the fraction of models in different clusters
with full parameter
range and on a subset of models with high production rate of a specific gene
 representing the over expression of the specific gene.
}
\section{Related Functions}{


\code{\link{sracipeSimulate}},  \code{\link{sracipeKnockDown}},
\code{\link{sracipeOverExp}},  \code{\link{sracipePlotData}},
}

\examples{
data("demoCircuit")
\dontrun{
rSet <- sRACIPE::sracipeSimulate(circuit = demoCircuit, numModels = 100,
plots=FALSE, plotToFile = FALSE)
rSet <- sRACIPE::sracipeNormalize(rSet)
}
}
