% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotMGPLoadings.R
\name{plotMGPLoadings}
\alias{plotMGPLoadings}
\title{Human MGP}
\usage{
plotMGPLoadings(obj, comp, type, colormap)
}
\arguments{
\item{obj}{runHumanMGP object}

\item{comp}{Which PLS component to visualize (default = 1)}

\item{type}{Type of plot (options: "single", "sum1", "sum2", "weighted"; default = "single")

single = plot individual loadings per gene (or per genetic PC); sum1 = abs. sum of genetic PCs per gene (stacked); sum2 = abs. sum of genetic PCs per gene (individual); weighted = abs. sum of genetic PCs, weighted by number of PCs, per gene}
}
\description{
Plot genetic loadings from "runHumanMGP" output
}
\author{
Hanne Hoskens
}
