# SurplusProductionModel 0.1.0

## Initial Release

### New Features
* Implemented Pella-Tomlinson surplus production model with RTMB integration
* Added state-space framework with process and observation error
* Included special cases for Schaefer (m=2) and Fox (m=1) models
* Created comprehensive data validation functions
* Implemented reference point calculations (MSY, BMSY, FMSY)
* Added model diagnostics and visualization tools
* Included example Ross Sea toothfish dataset
* Set up testing framework with >80% coverage target

### Package Structure
* Core S4 classes: `ProductionModel`
* Main functions: `fit_pella_tomlinson_model()`, `calculate_reference_points()`
* Data validation: `validate_cpue_data()`, `validate_catch_data()`
* Diagnostics: `plot_model_fit()`, `residual_analysis()`

### Documentation
* Complete roxygen2 documentation for all exported functions
* Getting started vignette
* Pella-Tomlinson theory vignette
* Comprehensive README with examples

### Testing
* Unit tests for all major functions
* Integration tests for full workflows
* Performance benchmarks
* Cross-validation tools

This release provides the MVP foundation for the ATO rTMB project, focusing on robust surplus production modeling for Antarctic toothfish stock assessment.