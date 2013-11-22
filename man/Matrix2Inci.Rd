\name{Matrix2Inci}
\alias{Matrix2Inci}
\title{
Presence/Absence data transform into incidence frequency data
}
\description{
Transfer prensence/absence data into incidence frequency data.
}
\usage{
Matrix2Inci(mat.data)
}
\arguments{
  \item{mat.data}{
  a numerical matrix or a data frame, the presence/absence of each species is recorded in the matrix. If the species is presence, it is record 1. On the contrary, it is record 0. See \code{\link[ChaoEntropy]{Seed_data}}
  }
}
\value{
  a numeric vector which first entry is total number of sampling units, and followed by species incidence frequency.
}
\author{
Y.T. Wang ; Y.H. Lee ; K.S. Tseng ; Anne Chao\cr
Maintainer : Anne Chao \email{chao@stat.nthu.edu.tw}
}
\seealso{
  \code{\link[ChaoEntropy]{ChaoEntropy}}
}
\examples{
# load the prensence/absence data
data(Seed_data)

# Transform the data format
Matrix2Inci(Seed_data)
ChaoEntropy(data=Matrix2Inci(Seed_data), datatype="incidence", method="all",
            se=TRUE, nboot=200, conf=0.95)
}
