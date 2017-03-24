\name{adjusted.KM}
\alias{adjusted.KM}
\title{
Adjusted Survival Curves by Using IPW and Kaplan-Meier Estimators.
}
\description{
This function allows for computing adjusted survival curves by weighting the individual contributions by the inverse of the probability to be in the group. The user enters individual survival data and the weights previously calculated (by using logistic regression for instance). The usual Kaplan-Meier estimator is adapted in order to obtain the adjusted survival.
}
\usage{
adjusted.KM(times, failures, variable, weights)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{times}{A numeric vector with the follow up times.}
  \item{failures}{A numeric vector with the event indicators (0=right censored, 1=event).}
  \item{variable}{A numeric vector with the binary variable under interest (only two groups).}
  \item{weights}{The weights for correcting the contribution of each individual. By default, the weights are all equalled to 1 and the survival curves correspond to the usual Kaplan-Meier estimator.}
  }
\details{
For instance, \code{w} may be equal to 1/\code{p}, where \code{p} is the estimated probability of the individual to be in its group. The probabilities \code{p} are often estimated by a logistic regression in which the dependent binary variable is the group. The possible confounding factors are the explanatory variables of this model. 
}
\value{
\item{Table}{This data frame presents the survival probabilities (\code{survival}) in each group (\code{variable}) according to the \code{times}. The number of individuals at risk (\code{n.risk}) and the number of observed events are also provided (\code{n.event}).}
}
\references{
Le Borgne F, Giraudeau B, Querard AH, Giral M and Foucher Y. Comparisons of the performances of different statistical tests for time-to-event analysis with confounding factors: practical illustrations in kidney transplantation. Statistics in medicine. 30;35(7):1103-16, 2016. <doi:10.1002/ sim.6777>

Cole S and Hernan M. Adjusted survival curves with inverse probability weights. Computer methods and programs in biomedicine, 75(1):45-49, 2004. <doi:10.1016/j.cmpb.2003.10.004>

}
\author{
Y. Foucher <Yohann.Foucher@univ-nantes.fr> and

F. Le Borgne <florent.le-borgne@etu.univ-nantes.fr>
}

\examples{

data(DIVAT)

# Kaplan-Meier estimator
res.km <- adjusted.KM(times=DIVAT$times, failures=DIVAT$failures,
  variable=DIVAT$ecd, weights=NULL)
plot(NULL,xlim=c(0,13),ylim=c(0,1),ylab="Graft and patient survival",
 xlab="Time post-transplantation (years)")
lines(res.km$times[res.km$variable==1], res.km$survival[res.km$variable==1],
type="s",col=2,lty=2,lwd=2)
lines(res.km$times[res.km$variable==0], res.km$survival[res.km$variable==0],
type="s",col=1,lty=2,lwd=2)

# adjusted Kaplan-Meier estimator by IPW
Pr0 <- glm(ecd ~ 1, family = binomial(link="logit"), data=DIVAT)$fitted.values[1]
Pr1 <- glm(ecd ~ age + hla + retransplant, data=DIVAT,
 family=binomial(link = "logit"))$fitted.values
W <- (DIVAT$ecd==1) * (1/Pr1) + (DIVAT$ecd==0) * (1)/(1-Pr1)
res.akm <- adjusted.KM(times=DIVAT$times, failures=DIVAT$failures,
 variable=DIVAT$ecd, weights=W)
lines(res.akm$times[res.akm$variable==1], res.akm$survival[res.akm$variable==1],
 type="s",col=2,lwd=2)
lines(res.akm$times[res.akm$variable==0], res.akm$survival[res.akm$variable==0],
 type="s",col=1,lwd=2)

nb.risk1<-function(x) {sum(DIVAT$times[DIVAT$ecd==0]>x)}
nb.risk2<-function(x) {sum(DIVAT$times[DIVAT$ecd==1]>x)}
segments(x0=0, y0=0.1, x1=13, y1=0.1) 
text(x=6, y=0.12, "number of at-risk patients", cex=0.8)
tps <- seq(1,12,by=1)
text(x=tps, y=rep(0.07,length(tps)), as.character(sapply(tps, FUN="nb.risk1")),
 cex=0.8, col=1)
text(x=tps, y=rep(0.02,length(tps)), as.character(sapply(tps, FUN="nb.risk2")),
 cex=0.8, col=2)
legend("topright", legend=c("Unadjusted estimator for SCD",
 "Adjusted estimator for SCD", "Unadjusted estimator for ECD",
 "Adjusted estimator for ECD"), col=c(1,1,2,2),
 lty=c(2,1,2,1), lwd=2, cex=0.8)
}
