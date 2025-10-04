#' Generic Functions for Production Model Package
#'
#' This file contains generic function definitions for methods used with
#' ProductionModel objects.
#'
#' @name generics
NULL

#' Extract Model Parameters
#'
#' Generic function to extract parameters from model objects.
#'
#' @param object A model object
#' @param ... Additional arguments passed to methods
#' @return Parameter values
#' @export
setGeneric("parameters", function(object, ...) standardGeneric("parameters"))

#' Extract Model Results
#'
#' Generic function to extract results from fitted model objects.
#'
#' @param object A model object
#' @param ... Additional arguments passed to methods
#' @return Model results
#' @export
setGeneric("results", function(object, ...) standardGeneric("results"))

#' Extract Model Data
#'
#' Generic function to extract data from model objects.
#'
#' @param object A model object
#' @param ... Additional arguments passed to methods
#' @return Model data
#' @export
setGeneric("model_data", function(object, ...) standardGeneric("model_data"))

#' Check Model Fitted Status
#'
#' Generic function to check if a model has been fitted.
#' Note: This creates a generic for the base R 'fitted' function.
#'
#' @param object A model object
#' @param ... Additional arguments passed to methods
#' @return Logical indicating fitted status
#' @export
if (!isGeneric("fitted")) {
  setGeneric("fitted", function(object, ...) standardGeneric("fitted"))
}
