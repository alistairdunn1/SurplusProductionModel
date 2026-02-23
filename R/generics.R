#' S3 Generic Functions for Production Model Package
#'
#' This file contains S3 generic function definitions for methods used with
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
parameters <- function(object, ...) UseMethod("parameters")

#' Extract Model Results
#'
#' Generic function to extract results from fitted model objects.
#'
#' @param object A model object
#' @param ... Additional arguments passed to methods
#' @return Model results
#' @export
results <- function(object, ...) UseMethod("results")

#' Extract Model Data
#'
#' Generic function to extract data from model objects.
#'
#' @param object A model object
#' @param ... Additional arguments passed to methods
#' @return Model data
#' @export
model_data <- function(object, ...) UseMethod("model_data")
