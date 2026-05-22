% Generated manually to match current exports
\name{project_forward}
\alias{project_forward}
\alias{print.pt_projection}
\title{Forward Projection Under Assumed Harvest Control}
\usage{
project_forward(
  model_fit,
  horizon = 10L,
  control,
  control_type = c("F", "catch"),
  n_sim = 1000L,
  process_error = c("historical", "normal", "none"),
  historical_window = 10L,
  seed = NULL,
  probs = c(0.05, 0.5, 0.95)
)

\method{print}{pt_projection}(x, ...)
}
\arguments{
\item{model_fit}{A fitted \code{ProductionModel} object.}

\item{horizon}{Integer projection length in years.}

\item{control}{Numeric scalar or vector of length \code{horizon}. Interpreted
as fishing mortality for \code{control_type = "F"} or total catch for
\code{control_type = "catch"}.}

\item{control_type}{Character, one of \code{"F"} or \code{"catch"}.}

\item{n_sim}{Number of Monte Carlo simulation trajectories.}

\item{process_error}{Character, one of \code{"historical"}, \code{"normal"},
or \code{"none"}.}

\item{historical_window}{Number of most recent historical years used when
\code{process_error = "historical"}. Default is 10.}

\item{seed}{Optional integer random seed.}

\item{probs}{Quantile probabilities for summarized outputs.}

\item{x}{A \code{pt_projection} object.}

\item{...}{Additional arguments (unused).}
}
\value{
\code{project_forward()} returns an object of class \code{pt_projection}
containing summary statistics, simulation arrays, and settings.
\code{print.pt_projection()} returns the object invisibly.
}
\description{
Projects a fitted model forward in time under assumed F or catch controls,
with uncertainty propagated through process error. By default, process
deviations are sampled from reconstructed historical process errors using the
last 10 years.
}
