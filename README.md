# Detecting respiratory chain deficiency in osteoblasts of older patients 

Data and analysis from Hipps et al. (2021, in preparation).

## Background
The analysis in this repository starts from tabular protein expression data, which is the output from image analysis described in the manuscript.  It calculates the proportion of osteoblasts which are respiratory chain deficient and generates some plots to visualise the results.

## Image analysis
Using Volocity various masks of pixels are created based upon the intensity values of different channels. Cells are identified using a nuclear marker and cytoplasm defined by the Osteocalcin channel marker (these channels are pre-labelled when combining the raw images into multi image tiffs using Nikon elements). Within this cell mask VDAC is used to define mitochondrial mass as an image mask.

Combining these two masks gives only VDAC positive pixels within osteoblast cells defined by the Osteocalcin cell marker. Using these mask volocity will create a csv file with intensity values for each individual cell’s mitochondrial mass positive pixels. It is possible to get minimum and maximum values, mean values, pixel count and standard deviation within the area. For the purpose of this research we used the mean values of intensity for analysis.

## Statistical analysis
To execute analysis: open R, set working directory to "analysis" directory in this repository, make sure that you have all the required dependencies (mclust, umap, MASS and sp) and execute the code in the "big_cytof.R" script.  

Executing the script should produce two .pdf files, one which contains the plots shown in Figure 3 & 4 (along with several pages of other plots: "02_Osteocalcin_VDAC_mask_n_stripped.pdf") and another which contains exactly Figure 4 from the manuscript ("Fig4.pdf").
