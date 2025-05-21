# Within\_Region\_091124\_CCR2

## Notes
* Label HF AFs as CCR2+ and CCR2- in macrophages in that same ROI. Then run SLIDE on the AFs.

### Input Files:
* Raw Data: started from begning becasue CCR2 is filtered out using the previous pre-processing object. 

### Scripts:
* getting_hi_low_af.R: the script that visualize ccr2 expression in HF mac regions, then subset AF regions into hi-CCR2 and low-CCR2.

### Output Files:
* low(hi)_ccr2_af.rds: rds files contain the Nanostring objects of low and hi ccr2 active fibroblast rds. 
* Data/x(y).csv: the raw expr matrix and y. 0 as low and 1 as hi CCR2.




