## Analysis code MMSplice manuscript
MMSplice predict variant effect on PSI, PSI5, PSI3, splicing efficiency and pathogenicity. This repository contains analysis code for MMSplice manuscript. MMSplice version `mmsplice==0.2.0`.

## Software dependencies
The analysis was done with python and ploting was done with R. Required packages are listed in env.yaml. 
1. Install miniconda or anaconda.
2. Run: `conda env create -f env.yaml`. This will install a new conda environment `mmaplice-manuscript`.
3. Activate the environment: `source activate mmaplice-manuscript`. Run the analysis code (in jupyter-notebooks) in this environment. 
4. Install COSSMO from https://github.com/PSI-Lab/COSSMO

## Folder structure
In general, the jupyter-notebooks were used to process, analysis the data and perform predictions. The result table are saved and read by R scripts for plotting. 
- `script/Figure2`: Code to generate figure for the Vex-seq and MFASS data. 
- `script/Figure3`: Code to generate figure for GTEx A5SS amd A3SS variants. 
- `script/Figure4`: Code to generate figure for MaPSy data. 
- `script/Figure5`: Code to generate figure for ClinVar data. Result table of other models were obtained from [Kipoi manuscript] (https://github.com/kipoi/manuscript)
- `data/`: contains raw and processed data
