# Translucency_prediction

This repository contains the code for the following paper: "Kiyokawa, H., Nagai, T., Yamauchi, Y., and Kim, J., The perception of translucency from surface gloss, Vision Research, 108140, (2022)". https://doi.org/10.1016/j.visres.2022.108140

The code runs on MATLAB (2021a) and uses some functions from Psychotoolbox (http://psychtoolbox.org/) and matlabPyrTools (https://github.com/LabForComputationalVision/matlabPyrTools)

# How to use
1. Download Psychotoolbox and matlabPyrTools and add their folders to matlab path.
2. Calculate image features. Run OF_mapping.m and Feature_Calculator_highlights_shading.m.
3. Predict human response. Run Model_prediction.m.
