\name{Count2Abun}
\alias{Count2Abun}
\title{
Frequency counts data transform into abundance data
}
\description{
Transfer frequency counts data into abundance data.
}
\usage{
Count2Abun(count.data)
}
\arguments{
  \item{count.data}{
  input species frequency counts data: a numerical matrix or a data frame of two
columns. The first column is the frequency j=1, 2...; and the second column is frequency counts (fj). \cr\cr
  Example: there are 59 singletons, 9 doubletons and so on in Insects_data:
    \tabular{rlllll}{
      [1,]  \tab 1  \tab 59 \cr
      [2,]  \tab 2  \tab 9 \cr
      [3,]  \tab 3  \tab 3 \cr
      [4,]  \tab 4  \tab 2 \cr
      [5,]  \tab 5  \tab 2 \cr
      [6,]  \tab 6  \tab 2 \cr
      [7,]  \tab 11  \tab 1 \cr
    }
  }
}
\value{
a numeric vector of species abundance data.
}
\author{
Y.T. Wang ; Y.H. Lee ; K.S. Tseng ; Anne Chao\cr
Maintainer : Anne Chao \email{chao@stat.nthu.edu.tw}
}
\seealso{
  \code{\link[ChaoEntropy]{ChaoEntropy}}
}
\examples{
# load the frequency counts data
data(Insects_data)

# Transform the data format
Count2Abun(Insects_data)
ChaoEntropy(data=Count2Abun(Insects_data), datatype="abundance", method="all",
            se=TRUE, nboot=200, conf=0.95)
}
