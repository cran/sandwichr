% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sandwich.model.r
\name{sandwich.model}
\alias{sandwich.model}
\title{Perform Sandwich model-based mapping}
\usage{
sandwich.model(object,
       sampling.attr,
       type="shp",
       ssh.id.col=NULL,
       ssh.weights=NULL)
}
\arguments{
\item{object}{When \code{type="shp"}, \code{object} is a list of three \code{sf} objects generated by \code{\link{load.data.shp}}, including a point \code{sf} object used as the sampling layer, a polygon \code{sf} object used as the SSH layer, and a polygon \code{sf} object used as the SSH layer. When \code{type="txt"}, \code{object} is a list of two data frames generated by \code{\link{load.data.txt}}, including a file linking sampling and SSH layers and a file linking reporting and SSH layers.}

\item{sampling.attr}{Text for the name of the attribute to be predicted in the sampling layer.}

\item{type}{Text for the type of input data. \code{type="shp"} denotes shapefiles, and \code{type="txt"} denotes text files. By default, \code{type="shp"}.}

\item{ssh.id.col}{Text for the column that specifies which stratum each sampling unit falls into (see \code{\link{load.data.txt}}). Set to \code{NULL} when \code{type="shp"}.}

\item{ssh.weights}{A \code{list} that specifies the strata in the SSH layer and their corresponding columns of weights in \code{reporting_ssh.file} (see \code{\link{load.data.txt}}).}
}
\value{
A \code{sandwich.ci} object that contains the estimated mean and standard deviation for each reporting unit.
}
\description{
Estimate the mean and standard error for each reporting unit using SSH-based spatial prediction.
}
\examples{
data(sim.data)
sim.sw <- sandwich.model(object=sim.data, sampling.attr="Value", type="shp")

}
\references{
Wang, J. F., Haining, R., Liu, T. J., Li, L. F., & Jiang, C. S. (2013). Sandwich estimation for multi-unit reporting on a stratified heterogeneous surface. \emph{Environment and Planning A}, 45(10), 2515-2534.
}
\seealso{
\code{\link{load.data.shp}}, \code{\link{load.data.txt}}
}
