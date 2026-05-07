# Overview
As part of my Master’s Thesis work, I derived relevant features from UAV LiDAR and multispectral data, combined them, and subsequently applied a Random Forest algorithm to map standing and downed deadwood in the Harz region of Germany.
Multiple LiDAR-derived structural features were generated using the lidR package, with all processing conducted in R. In parallel, several vegetation indices were calculated from the multispectral imagery in JupyterLab. These datasets were then combined for further analysis, including classification performance assessment, feature importance analysis, and the generation of prediction maps in JupyterLab.
Classification performance was evaluated across three datasets: LiDAR-only, multispectral-only, and fusion (LiDAR + multispectral), and three scenarios: standing deadwood, downed deadwood, and a multi-class scenario including both deadwood types.
Overall, the study highlights the potential of UAV-based LiDAR-multispectral fusion for improving deadwood detection and spatial mapping in disturbed forests. 

# Data
The UAV data was pre-processed by the Geoinformation in Environmental Planning Lab team (Technische Universität Berlin) and then handed over to me.

# Code
Order of code
1. *R script for LiDAR metrics.R*: Code for deriving LiDAR-based structural features
2. *Additional VIs.ipynb*: Code for calculating different VIs relevant to deadwood detection
3. *Fusion model.ipynb*: Code for LiDAR-multispectral fusion model alongisde further analysis for deadwood mapping 

I have added extensive comments and structured the code with clear Markdown sectioning to improve readability and organization.

# Contact
Feel free to connect with me on LinkedIn: www.linkedin.com/in/kirankaphleep


