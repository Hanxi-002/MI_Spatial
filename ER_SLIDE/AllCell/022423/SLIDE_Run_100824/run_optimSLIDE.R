library(SLIDE)

yaml_path = "/ix/djishnu/Hanxi/MI_Spatial/ER_SLIDE/AllCell/022423/SLIDE_Run_100824/run_SLIDE.yaml"
input_params <- yaml::yaml.load_file(yaml_path)

SLIDE::checkDataParams(input_params)
SLIDE::optimizeSLIDE(input_params, sink_file = FALSE)
