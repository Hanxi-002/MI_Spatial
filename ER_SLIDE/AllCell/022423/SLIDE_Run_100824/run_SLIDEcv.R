library(SLIDE)

#yaml_path = "/ix/djishnu/Hanxi/MI_Spatial/ER_SLIDE/AllCell/022423/run_SLIDE.yaml"
#input_params <- yaml::yaml.load_file(yaml_path)

SLIDEcv("/ix/djishnu/Hanxi/MI_Spatial/ER_SLIDE/AllCell/022423/SLIDE_Run_100824/run_SLIDEcv.yaml",nrep = 25, k = 10)
