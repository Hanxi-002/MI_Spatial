# AllCell_HFvControl_011123

## Background
This analysis is repeated to include more genes than the previous experiment [here](https://github.com/Hanxi-002/Multiomics_Integration/tree/main/Dutta_Spatial/AllCell/121922). In the previous experiment, I used a threshold of 5% to exclude genes which resulted in ~ 3200 genes. In this experiment, we are using the threshold of 3% which resulted in around 4300 genes. 


## Links
* [ER Input Data](https://pitt-my.sharepoint.com/personal/xiaoh_pitt_edu/_layouts/15/onedrive.aspx?id=%2Fpersonal%2Fxiaoh%5Fpitt%5Fedu%2FDocuments%2FMultiOmic%2FDutta%5FSpatial%2FER%5FSLIDE%2FAllCell%2F011123%2FData) 

* [VanillaER Output]()

* [SLIDE Output]()

## Notes
1. We didn't run the entire pipeline. From the previous analysis, we learned that the number of genes varies a lot by different delta increments. And the large cluster numbers give rise to the long computational times for cross validation comparison for different delta and lambdas. For this analysis, we are using point delta values and running pipeline step 3 to simply get the cross validation performance between the model and the permutation.

