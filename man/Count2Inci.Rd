\name{Count2Inci}
\alias{Count2Inci}
\title{
Incidence counts data transform into incidence frequency data
}
\description{
Transfer incidence counts data into incidence frequency data.
}
\usage{
Count2Inci(count.data, t)
}
\arguments{
  \item{count.data}{
  input species incidence counts data: a numerical matrix or a data frame of two
columns. The first column is the frequency j=1, 2...; and the second column is incidence counts (Qj). \cr\cr
  Example: there are 39 singletons, 18 doubletons and so on in Ant_data.count:
    \tabular{rlllll}{
      [1,]  \tab 1  \tab 39 \cr
      [2,]  \tab 2  \tab 18 \cr
      [3,]  \tab 3  \tab 13 \cr
      [..,]  \tab ...  \tab ... \cr
      [16,]  \tab 19  \tab 2 \cr
      [17,]  \tab 20  \tab 1 \cr
      [18,]  \tab 29  \tab 1 \cr
    }
  }
  \item{t}{
  the total number of sampling units. \cr
  e.g. See \code{\link[ChaoEntropy]{Ant_data}}: \code{t = 62} 
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
# load the incidence counts data
data(Ant_data.count)
t <- 62


# Transform the data format
Count2Inci(Ant_data.count, t)
ChaoEntropy(data=Count2Inci(Ant_data.count, t), datatype="incidence", method="all",
            se=TRUE, nboot=200, conf=0.95)
}
