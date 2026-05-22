% Generated manually to match current exports
\name{set_reference_point_defaults}
\alias{set_reference_point_defaults}
\title{Set Package-Level Reference Point Defaults}
\usage{
set_reference_point_defaults(biomass_target = NULL, baseline = c("auto", "B0", "K"))
}
\arguments{
\item{biomass_target}{Optional numeric vector of default biomass fractions
to use when a function call omits \code{biomass_target}. Pass
\code{NULL} explicitly to clear the default.}

\item{baseline}{Optional default biomass baseline, one of
\code{"auto"}, \code{"B0"}, or \code{"K"}.}
}
\value{
Invisibly returns the updated defaults as a named list.
}
\description{
Sets package-wide defaults for depletion target reporting. These defaults
are used by functions such as \code{calculate_reference_points()},
\code{profile_likelihood()}, \code{bayesian_fit()},
\code{summary.ProductionModel()}, and \code{print.ProductionModel()}
when \code{biomass_target} and/or \code{baseline} are not supplied.
}
