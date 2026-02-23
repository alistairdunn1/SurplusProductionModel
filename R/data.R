#' Simulated Ross Sea Toothfish Catch Data
#'
#' Simulated catch data for Antarctic toothfish (\emph{Dissostichus mawsoni})
#' in the Ross Sea region (CCAMLR Area 88.1/88.2), spanning 1997--2024.
#' The catch trajectory reflects a fishery that developed in the late 1990s,
#' expanded through the 2000s, and was progressively constrained by catch limits.
#'
#' @format A data frame with 28 rows and 2 columns:
#' \describe{
#'   \item{year}{Integer, calendar year (1997--2024)}
#'   \item{catch}{Numeric, total catch in tonnes}
#' }
#'
#' @source Simulated based on approximate characteristics of the Ross Sea
#'   toothfish fishery. Not real data.
#'
#' @examples
#' data(ross_sea_catch)
#' head(ross_sea_catch)
#' plot(ross_sea_catch$year, ross_sea_catch$catch,
#'   type = "b",
#'   xlab = "Year", ylab = "Catch (t)"
#' )
"ross_sea_catch"

#' Simulated Ross Sea Toothfish CPUE Data (Multi-Index)
#'
#' Simulated CPUE (catch per unit effort) indices for Antarctic toothfish
#' in the Ross Sea region, containing two abundance indices: a longline
#' fishery CPUE series and a biennial trawl survey index.
#'
#' @format A data frame with approximately 40 rows and 3 columns:
#' \describe{
#'   \item{year}{Integer, calendar year}
#'   \item{cpue}{Numeric, CPUE index value (units vary by index)}
#'   \item{label}{Character, index identifier (\code{"Longline"} or \code{"Trawl"})}
#' }
#'
#' @source Simulated based on approximate characteristics of Ross Sea
#'   toothfish abundance indices. Not real data.
#'
#' @seealso \code{\link{ross_sea_cpue_simple}} for a single-index version
#'
#' @examples
#' data(ross_sea_cpue)
#' head(ross_sea_cpue)
"ross_sea_cpue"

#' Simulated Ross Sea Toothfish CPUE Data (Single Index)
#'
#' Simplified single-index CPUE time series for Antarctic toothfish
#' in the Ross Sea, containing only the longline fishery CPUE.
#' Suitable for basic single-area model fitting examples.
#'
#' @format A data frame with 27 rows and 2 columns:
#' \describe{
#'   \item{year}{Integer, calendar year (1998--2024)}
#'   \item{cpue}{Numeric, longline CPUE (kg/hook)}
#' }
#'
#' @source Simulated. See \code{\link{ross_sea_cpue}} for multi-index version.
#'
#' @examples
#' data(ross_sea_cpue_simple)
#' plot(ross_sea_cpue_simple$year, ross_sea_cpue_simple$cpue,
#'   type = "b",
#'   xlab = "Year", ylab = "CPUE (kg/hook)"
#' )
"ross_sea_cpue_simple"

#' True Simulation Parameters for Ross Sea Toothfish Dataset
#'
#' The parameter values used to generate the simulated Ross Sea toothfish
#' datasets. Useful for validating model estimation by comparing fitted
#' parameters against known true values.
#'
#' @format A named list with elements:
#' \describe{
#'   \item{r}{Intrinsic growth rate (0.08)}
#'   \item{K}{Carrying capacity in tonnes (80000)}
#'   \item{m}{Pella-Tomlinson shape parameter (2.0, i.e. Schaefer model)}
#'   \item{q_LL}{Longline catchability (0.00035)}
#'   \item{q_Trawl}{Trawl survey catchability (0.0008)}
#'   \item{sigma_proc}{Process error SD (0.05)}
#'   \item{sigma_obs_LL}{Longline observation error SD (0.20)}
#'   \item{sigma_obs_Trawl}{Trawl observation error SD (0.25)}
#'   \item{B0}{Initial biomass in tonnes (72000)}
#' }
#'
#' @source Simulation parameters chosen to approximate Antarctic toothfish
#'   life-history characteristics.
#'
#' @examples
#' data(ross_sea_true_params)
#' str(ross_sea_true_params)
"ross_sea_true_params"

#' True Biomass Trajectory for Ross Sea Toothfish Simulation
#'
#' The true (simulated) biomass trajectory used to generate the Ross Sea
#' toothfish datasets. Useful for comparing model-estimated biomass
#' against known truth.
#'
#' @format A data frame with 28 rows and 2 columns:
#' \describe{
#'   \item{year}{Integer, calendar year (1997--2024)}
#'   \item{biomass}{Numeric, true biomass in tonnes}
#' }
#'
#' @source Simulated from a Pella-Tomlinson model with process error.
#'
#' @examples
#' data(ross_sea_true_biomass)
#' plot(ross_sea_true_biomass$year, ross_sea_true_biomass$biomass,
#'   type = "l",
#'   xlab = "Year", ylab = "Biomass (t)", main = "True Biomass Trajectory"
#' )
"ross_sea_true_biomass"
