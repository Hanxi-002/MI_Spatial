## MI\_Spatial / ALLFibro

###Notes
  * HF v Conrol in all fibroblast.

###Input
  * Geomx V3 (analyzed Geomx data)
  * Data/x.csv: x matrix for SLIDE
  * Data/y.csv: y matrix for SLIDE

###Scripts
  * get_data.R: extract x an y from GeoMx\_V3.RDS.
  * run_SLIDE.yaml: yaml file for the first step of SLIDE (optimize SLIDE).
  * run_optimSLIDE.R: R code to run optimizeSLIDE.
  * optimSLIDE_shell.sh: submit run_optimSLIDE.R as a job.
  * run_SLIDEcv.yaml: yamil file for SLIDEcv.
  * run_SLIDEcv.R: R code to run SLIDEcv.
  * SLIDEcv_shell.sh: submit run_SLIDEcv.R as a job.

### Output
  * reults/SLIDE_results: all results from running SLIDE (optimizeSLIDE and SLIDEcv).