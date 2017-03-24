\name{DIVAT}
\alias{DIVAT}
\docType{data}
\title{
A Sample From the DIVAT Data Bank of Kidney Transplant Recipients.
}
\description{
The sample were extracted from the DIVAT cohort. The nephrology department of the Nantes University Hospital has developed DIVAT: a database which includes the monitoring of medical records for kidney transplant recipients. For more details or requesting the data for your research, please visit the following website: www.divat.fr/en.
}
\usage{data(DIVAT)}
\format{
A data frame with the 4 following variables:
  \describe{
    \item{\code{age}}{This numeric vector provides the age of the recipient at the transplantation (in years).}
	\item{\code{hla}}{This numeric vector provides the indicator of transplantations with at least 4 HLA incompatibilities between the donor and the recipient (1 for high level and 0 otherwise).}
    \item{\code{retransplant}}{This numeric vector provides the indicator of re-transplantation (1 for more than one transplantation and 0 for first kidney transplantation).}
    \item{\code{ecd}}{The Expended Criteria Donor (1 for transplantations from ECD and 0 otherwise).}
    \item{\code{times}}{This numeric vector is the follow up times of each patient.}
    \item{\code{failures}}{This numeric vector is the event indicator (0=right censored, 1=event). An event is considered when return in dialysis or patient death with functioning graft is observed.}
  }
}
\details{
This database is constituted by 1912 patients transplanted between January 1996 and December 2013 and followed in the prospective DIVAT cohort from Nantes University hospital. Adults receiving a first or second kidney transplant alone from a deceased heart beating donor were included. In kidney transplantation, donors are classified into two categories, the expanded criteria donor and the standard criteria donor. ECD are defined by widely accepted criteria, which includes donors older than 60 years of age or 50-59 years of age with two of the following characteristics: history of hypertension, cerebrovascular accident as the cause of death or terminal serum creatinine higher than 1.5 mg/dL.
}
\references{
Le Borgne F, Giraudeau B, Querard AH, Giral M and Foucher Y. Comparisons of the performances of different statistical tests for time-to-event analysis with confounding factors: practical illustrations in kidney transplantation. Statistics in medicine. 30;35(7):1103-16, 2016. <doi:10.1002/ sim.6777>
}

\examples{

data(DIVAT)

# Compute the non-adjusted Cox PH model
cox.ECD0<-coxph(Surv(times, failures) ~ ecd, data=DIVAT)
summary(cox.ECD0) # Hazard Ratio = 1.97
}
