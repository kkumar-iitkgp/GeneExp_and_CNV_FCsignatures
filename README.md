# GeneExp_and_CNV_FCsignatures
 
Scripts for reproducing results for association between AHBA Gene expression spatial patterns and FC signatures for 16p11.2
deletion and 22q11.2 deletion. 
Please cite: [Moreau, Clara*, Urchs, Sebastian*, et al. "Neuropsychiatric mutations delineate functional brain connectivity dimensions contributing to autism and schizophrenia." BioRxiv (2019): 862615](https://www.biorxiv.org/content/10.1101/862615v1).


## Dependencies
The code is written and tested in [MATLAB](https://www.mathworks.com/products/matlab.html) R2019b, and also includes the data required to reproduce the reported analysis and stats. The scripts can be run directly to reproduce the Partial Least Square Regression (PLSR) analysis and Correlation per Gene analysis. In addition a R script is included to make Histogram plots using [ggplot2](https://ggplot2.tidyverse.org/) and [ggrepel](https://cran.r-project.org/web/packages/ggrepel/vignettes/ggrepel.html).

## Analysis
To call Partial Least Square Regression (PLSR) and Correlation Per Gene analysis.
>**_Run:_** `script_call_PLSR_and_CorrPerGene.m`

This will call the following scripts: 

>**_1:_** `code/script1_call_PLSR_nodal_and_regional.m`

and then 

>**_2:_** `code/script2_call_CorrPerGene.m`

The PLSR results (Percentage Variance (PCTVAR) explained and p-value), as well as CorrPerGene results (Pearson r, p-value, FDR p-value) are saved as .xlsx files in the `data` folder. 
Finally, the histogram plots for CorrPerGene with Gene names for 16p11.2 (or 22q11.2) region genes can be made using 

>**_3:_** `code/scriptR_plot_HistogramCorrPerGene_16p22q.R`

and are saved in `plots` folder. Most of the data as well as output of the analysis are already included in `data` and `plots` folders for reference. 

## Additional citation for references
GeneExpression data: We used abagen toolbox ([Markello et al. 2020](https://abagen.readthedocs.io/en/stable/index.html)) and the guidelines in (Arnatkevic̆iūtė, Aurina et al. 2019) to align gene expression data in the adult human cortex from the Allen Human Brain Atlas ([AHBA](https://human.brain-map.org/)) dataset (Hawrylycz et al. 2012) to [MIST64 brain parcellation](https://doi.org/10.12688/mniopenres.12767.2). Please cite the following references if you use the gene expression data:
* Markello, Ross, Shafiei, Golia, Zheng, Ying-Qiu, Mišić, Bratislav. abagen: A toolbox for the Allen Brain Atlas genetics data. Zenodo; 2020. Available from: https://doi.org/10.5281/zenodo.3726257. 
* Arnatkevic̆iūtė, Aurina, Fulcher, Ben D., Fornito, Alex (2019). A practical guide to linking brain-wide gene expression and neuroimaging data. NeuroImage, 189, 353-367. https://doi.org/10.1016/j.neuroimage.2019.01.011
* Hawrylycz, Michael J., et. al (2012). An anatomically comprehensive atlas of the adult human brain transcriptome. Nature, 489(7416), 391–399. https://doi:10.1038/nature11405

In addition, for PLSR analysis we reffered to [Morgan et al 2019](https://github.com/SarahMorgan/Morphometric_Similarity_SZ). 
* Morgan, Sarah E., et al. "Cortical patterning of abnormal morphometric similarity in psychosis is associated with brain expression of schizophrenia-related genes." Proceedings of the National Academy of Sciences 116.19 (2019): 9604-9609. https://doi.org/10.1073/pnas.1820754116

