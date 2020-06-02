# GeneExp_and_CNV_FCsignatures
 
Scripts for reproducing results for association between AHBA Gene expression spatial patterns and FC signatures for 16p11.2
deletion and 22q11.2 deletion. 
Please cite: [Moreau, Clara, et al. "Neuropsychiatric mutations delineate functional brain connectivity dimensions contributing to autism and schizophrenia." BioRxiv (2019): 862615](https://www.biorxiv.org/content/10.1101/862615v1).

## Dependencies
The code is written in [MATLAB](https://www.mathworks.com/products/matlab.html), and also includes the data required to reproducing the analysis. The scripts can be run directly to reproduce the Partial Least Square Regression (PLSR) analysis and Correlatio per Gene analysis. In addition a R script is included to make Histogram plots using [ggplot2](https://ggplot2.tidyverse.org/) and [ggrepel](https://cran.r-project.org/web/packages/ggrepel/vignettes/ggrepel.html).
## Additional citation for references
GeneExpression data: We used abagen toolbox ([Markello et al. 2020](https://abagen.readthedocs.io/en/stable/index.html)) and the guidelines in (Arnatkevic̆iūtė, Aurina et al. 2019) to align gene expression data in the adult human cortex from the Allen Human Brain Atlas ([AHBA](https://human.brain-map.org/)) dataset (Hawrylycz et al. 2012) to MIST64 brain parcellation. Please cite the following references if you use the gene expression data:
* Markello, Ross, Shafiei, Golia, Zheng, Ying-Qiu, Mišić, Bratislav. abagen: A toolbox for the Allen Brain Atlas genetics data. Zenodo; 2020. Available from: https://doi.org/10.5281/zenodo.3726257. 
* Arnatkevic̆iūtė, Aurina, Fulcher, Ben D., Fornito, Alex (2019). A practical guide to linking brain-wide gene expression and neuroimaging data. NeuroImage, 189, 353-367. https://doi.org/10.1016/j.neuroimage.2019.01.011
* Hawrylycz, Michael J., et. al (2012). An anatomically comprehensive atlas of the adult human brain transcriptome. Nature, 489(7416), 391–399. https://doi:10.1038/nature11405

In addition, for PLSR analysis we reffered to [Morgan et al 2019](https://github.com/SarahMorgan/Morphometric_Similarity_SZ). 
* Morgan, Sarah E., et al. "Cortical patterning of abnormal morphometric similarity in psychosis is associated with brain expression of schizophrenia-related genes." Proceedings of the National Academy of Sciences 116.19 (2019): 9604-9609. https://doi.org/10.1073/pnas.1820754116

