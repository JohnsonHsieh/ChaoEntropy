\name{ChaoEntropy}
\alias{ChaoEntropy}
\title{
Estimation of Shannon entropy
}
\description{
\kbd{ChaoEntropy} is a function to provide a new Shannon entropy estimator proposed by Chao et al. (2013). \cr\cr
For individual-based (abundance) data, there are five other estimators which are  Chao and Shen (2003) estimator; Grassberger (2003) estimator; Zhal (1977) jackknife estimator; Zhang (2012) Hz* estimator; and the observed provided in comparison. \cr\cr
For sample-based (incidence) data, new Shannon entropy estimator is in comparison with the observed entropy.
}
\usage{
ChaoEntropy(data, datatype = c("abundance", "incidence"), method = c("all", "Chao",
            "ChaoShen", "Grassberger", "Jackknife", "Zhang", "Observed"),
            se = TRUE, nboot = 200, conf = 0.95)
}

\arguments{
  \item{data}{
  a vector of species abundance or incidence frequency. If \code{datatype = "incidence"}, then the input format of first entry should be total number of sampling units, and followed by species incidence frequency.
}
  \item{datatype}{
  the data type of input data. That is individual-based abundance data (\code{datatype = "abundance"}) or sample-based incidence data (\code{datatype = "incidence"}).
}
  \item{method}{
  the \code{method} constructed to estimate entropy (see Details)
}
  \item{se}{
  calculate bootstrap standard error and show confidence interval; default is \code{TRUE}.
}
  \item{nboot}{
  the number of bootstrap resampling times, default is \code{200}.
}
  \item{conf}{
  a positive number \eqn{\le} 1. "\code{conf}" specifies the confidence level for confidence interval. The default is \code{0.95}.
}
}
\details{
  \enumerate{
    \item If \code{datatype = "abundance"}:
    \itemize{
      \item \code{method = "all"} : all estimators below involved. The default is "\code{all}"
      \item \code{method = "Chao"} estimator, see Chao et al. (2013)
      \item \code{method = "ChaoShen"} estimator, see Chao and Shen (2003)
      \item \code{method = "Grassberger"} estimator, see Grassberger (2003)
      \item \code{method = "Jackknife"} estimator, see  Zhal (1977)
      \item \code{method = "Zhang"} estimator, see Zhang (2012)
      \item \code{method = "Observed"} estimator, the observed entropy estimator
    }
    \item If \code{datatype = "incidence"}: 
    \itemize{
      \item \code{method = "all"} : all estimators below involved. The default is "\code{all}"
      \item \code{method = "Chao"} estimator, see Chao et al. (2013) in Appendix S6
      \item \code{method = "Observed"} estimator, the observed entropy estimator
    }
  }
}

\value{
\kbd{ChaoEntropy} returns a table of various entropy estimators, their standard error and 95 \% confidence interval which the \code{method} you choose.
}
\references{
Chao, A., Wang, Y.T. & Jost, L. (2013) Entropy and the species accumulation curve: a novel entropy estimator via discovery rates of new species. To appear in Methods in Ecology and Evolution. \cr\cr
Chao, A. & Jost, L. (2012) Coverage-based rarefaction and extrapolation: standardizing samples by completeness rather than size. \emph{Ecology}, \bold{93}, 2533-2547. \cr\cr
Chao, A. & Shen, T.J. (2003) Nonparametric estimation of Shannon's index of diversity when there are unseen species. \emph{Environmental and Ecological Statistics}, \bold{10}, 429-443. \cr\cr
Grassberger, P. (2003) Entropy estimates from insufficient samplings. URL \url{www.arxiv.org.} arXiv:physics/0307138v2. Updated 2008. \cr\cr
Zahl, S. (1977) Jackknifing an index of diversity. \emph{Ecology}, \bold{58}, 907-913. \cr\cr
Zhang, Z. (2012) Entropy estimation in Turing's perspective. \emph{Neural Computation}, \bold{24}, 1368-1389. \cr\cr
Chao, A. (1984) Nonparametric estimation of the number of classes in a population. \emph{Scandinavian Journal of Statistics}, \bold{11}, 265-270. \cr\cr
Chao, A. (1987) Estimating the population size for capture-recapture data with unequal catchability. \emph{Biometrics}, \bold{43}, 783-791.
}
\author{
Y.T. Wang ; Y.H. Lee ; K.S. Tseng ; Anne Chao\cr
Maintainer : Anne Chao \email{chao@stat.nthu.edu.tw}
}
\seealso{
  \code{\link[ChaoEntropy]{Count2Abun}}, 
  \code{\link[ChaoEntropy]{Count2Inci}}, 
  \code{\link[ChaoEntropy]{Matrix2Inci}}
}

\examples{
# load the individual-base (abundacne) data
data(Spider_data)

# Estimation of Shannon entropy
ChaoEntropy(Spider_data, datatype="abundance", method="all", se=TRUE, nboot=200, conf=0.95)

# load the sample-base (incidence) data
data(Ant_data)

# Estimation of Shannon entropy
ChaoEntropy(Ant_data, datatype="incidence", method="all", se=TRUE, nboot=200, conf=0.95)
}
\keyword{entropy}
\keyword{Diversity}
