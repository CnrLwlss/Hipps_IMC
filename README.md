# Detecting respiratory chain deficiency in osteoblasts of older patients 
## Hipps et al. (2021)

Data and analysis from Hipps et al. (2021)(in preparation)

This analysis starts from tabular protein expression data, which is the output from image analysis described in the manuscript.  It calculates the proportion of osteoblasts which are RC deficient and generates some plots to visualise the results.

To execute analysis: open R, set working directory to "analysis" directory in this repository, make sure that you have all the required dependencies (mclust, umap, MASS and sp) and execute the code in the "big_cytof.R" script.  

Executing the script should produce two .pdf files, one which contains the plots shown in Figure 3 & 4 (along with several pages of other plots: "02_Osteocalcin_VDAC_mask_n_stripped.pdf") and another which contains exactly Figure 4 from the manuscript ("Fig4.pdf").