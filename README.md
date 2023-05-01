# TME-mutational-signature-interaction

Introduction

Research on the associations between mutational signatures and tumor microenvironment (TME) in a pan-cancer context is limited. To address this, we analyzed over 8,000 TCGA tumor samples, using machine learning to explore the relationship between mutational signatures and TME, and develop a risk score for predicting patient survival outcomes. We also constructed an interaction model to examine their influence on cancer prognosis. Our analysis revealed varying associations between mutational signatures and TME. Risk scores based on specific mutational signatures showed strong pan-cancer survival stratification ability. We proposed a novel approach to predict TME cell types using genome-derived mutational signatures when transcriptome data is unavailable. Certain mutational signatures and their interaction with immune cells significantly impact clinical outcomes in specific cancer types. In conclusion, our study highlights the importance of considering both mutational signatures and immune phenotypes in cancer research and their implications for personalized cancer treatments and immunotherapy.


## Analysis Contents

- `Linear_regression_model.R`: R script for discovering correlations between mutational signatures and the tumor microenvironment.
- `get_riskscore.R`: R script for calculating overall survival risk score based on immune infiltration-associated mutational signatures.
-  `predicted_TME.R`: R script for predicting the level of immune infiltration using cancer mutational signatures.
- `identified_interaction_effects.R`: R script for performing statistical interaction tests to identify survival effects of immunophenotypes altered by mutational signatures.





## Citation
If you use data, results or conclusion from this work, please cite:
Luo et al. Unveiling the Interplay between Mutational Signatures and Tumor Microenvironment: A Pan-Cancer Analysis.

## LICENSE

GPL-3
