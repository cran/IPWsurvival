\name{IPWsurvival-package}
\alias{IPWsurvival-package}
\alias{IPWsurvival}
\docType{package}
\title{
Adjusted Kaplan-Meier Estimator and Log-Rank Statistic.
}
\description{
This package computes adjusted Kaplan-Meier estimator and log-rank statistic by using inverse probability weighting (IPW).
}
\details{
In observational studies, the presence of confounding factors is common and the comparison of different groups of subjects requires adjustment. In the presence of survival data, this adjustment can be achieved with a multivariate model. A recent alternative solution is the use of adjusted survival curves and log-rank test based on inverse probability weighting (IPW). By using the approach proposed by Xie and Liu (2005), we illustrate the usefulness of such methodology by studying the patient and graft survival of kidney transplant recipients according to the expanded donor criteria (ECD).

\tabular{ll}{
Package: \tab IPWsurvival\cr
Type: \tab Package\cr
Version: \tab 0.5\cr
Date: \tab 2017-03-20\cr
License: \tab GPL (>=2)\cr
LazyLoad: \tab yes\cr
}

\tabular{rl}{
adjusted.KM \tab Compute adjusted survival curves by weighting the indvidual contributions  \cr
 \eqn{  } \tab by the inverse of the probability to be in the group. \cr
adjusted.LR \tab Propose the log-rank test adapted to the corresponding adjusted \cr
 \eqn{  }  \tab  survival curves. \cr
DIVAT \tab These data extracted from the DIVAT data bank of kidney \cr
 \eqn{  }  \tab  transplant recipients. \cr
}
}
\author{
Y. Foucher <Yohann.Foucher@univ-nantes.fr>

and F. Le Borgne <fleborgne@idbc.fr>
}
\references{
Le Borgne F, Giraudeau B, Querard AH, Giral M and Foucher Y. Comparisons of the performances of different statistical tests for time-to-event analysis with confounding factors: practical illustrations in kidney transplantation. Statistics in medicine. 30;35(7):1103-16, 2016. <doi:10.1002/ sim.6777>

Cole S and Hernan M. Adjusted survival curves with inverse probability weights. Computer methods and programs in biomedicine, 75(1):45-49, 2004. <doi:10.1016/j.cmpb.2003.10.004>

Jun Xie and Chaofeng Liu. Adjusted Kaplan-Meier estimator and log-rank test with inverse probability of treatment weighting for survival data. Statistics in medicine, 24(20):3089-3110, 2005. <doi:10.1002/sim.2174>
}

\seealso{ URL: http://www.labcom-risca.fr }
