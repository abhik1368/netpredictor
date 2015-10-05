\name{netpredictor-package}
\alias{netpredictor-package}
\alias{netpredictor}
\docType{package}
\title{Link prediction and analysis of Bipartite Networks}
\description{
The package analyses the Unipartite and Bipartite networks and also algorithms have 
been implemented for link prediction like Random walks with Restart, Network Based Inference and Heterogenous Graph based inference algorithms . More Graph based Machine Learning algorithms algorithms will be included in the upcoming versions.
}
\details{
\tabular{ll}{
Package: \tab netpredicter\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2015-10-02\cr
License: \tab MIT\cr
}
To install the netpredictor

libary(devtools)

install_github("netpredictor", "abhik1368")

}
\author{
Abhik Seal

Maintainer: abhik1368@gmail.com
}
\references{
\itemize{
\item {Chen X, et al. Drug–target interaction prediction by random walk on the heterogeneous network. Mol. BioSyst 2012;8:1970-1978.}
\item {Vanunu O, Sharan R. Proceedings of the German Conference on Bioinformatics. Germany: GI; 2008. A propagation-based algorithm for inferring gene-disease assocations; pp. 54–63.}
\item {Wang W, et al. Drug target predictions based on heterogeneous graph inference. Pac. Symp. Biocomput. 2013:53-64}
\item {Cheng F, et al. Prediction of drug-target interactions and drug repositioning via network-based inference. PLoS Comput. Biol. 2012;8:e1002503.}
\item {Borgatti, S. P. and Everett, M. G. (1997) Network analysis of 2--mode data. Social Networks \bold{19}, 243--269.}
 }
}

\seealso{
\code{\link{igraph},\link{biNetwalk},\link{nbiNet},\link{uNetwalk},\link{sig.net}}
}
