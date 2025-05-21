library(SLIDE)

#yaml_path = "/ix/djishnu/Hanxi/MI_Spatial/ER_SLIDE/Within_Region/091124_hf_CCR2/run_SLIDEcv.yaml"
#input_params <- yaml::yaml.load_file(yaml_path)

SLIDEcv("/ix/djishnu/Hanxi/MI_Spatial/ER_SLIDE/Within_Region/091124_hf_CCR2/run_SLIDEcv.yaml",nrep = 50, k = 10)