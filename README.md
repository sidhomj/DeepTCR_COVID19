# Deep learning identifies antigenic determinants of severe SARS-CoV-2 infection within T-cell repertoires

## About
SARS-CoV-2 infection is characterized by a highly variable clinical course with patients experiencing asymptomatic infection all the way to requiring critical care support. This variation in clinical course has led physicians and scientists to study factors that may predispose certain individuals to more severe clinical presentations in hopes of either identifying these individuals early in their illness or improving their medical management. We sought to understand immunogenomic differences that may result in varied clinical outcomes through analysis of T-cell receptor sequencing (TCR-Seq) data in the open access ImmuneCODE database. We identified two cohorts within the database that had clinical outcomes data reflecting severity of illness and utilized DeepTCR, a multiple-instance deep learning repertoire classifier, to predict patients with severe SARS-CoV-2 infection from their repertoire sequencing. We demonstrate that patients with severe infection have repertoires with higher T-cell responses associated with SARS-CoV-2 epitopes and identify the epitopes that result in these responses. Our results provide evidence that the highly variable clinical course seen in SARS-CoV-2 infection is associated to certain antigen-specific responses.

## Publication
For full description of methods, please refer to the following manuscript:

[Sidhom, JW., Baras, A.S. Deep learning identifies antigenic determinants of severe SARS-CoV-2 infection within T-cell repertoires. Sci Rep 11, 14275 (2021).](https://www.nature.com/articles/s41598-021-93608-8)

## Prepare Data
Since the repertoire data is large, it is not included in this repository. To replicate the anlayses found here, follow the following instructions.
 
 - Download the full repertoire data [here](https://clients.adaptivebiotech.com/pub/covid-2020). There is a file (INCOV067-AC-3_TCRB.tsv) that is corrupted from the download link (as of the date this repository was created) that needs to be downloaded separately from the ImmuneAccess portal and used to replace the corrupted one. 
 - All downloaded repertoire data should be put into a directory labeled "repertoires_raw" underneath Data/ImmuneCODE directory.
 - Run reformat_files.py under scripts/preprocessing to reformat the repertoire data to inputs that DeepTCR can use. Doing this should create another directory underneath Data/ImmuneCODE called repertoires.
 - Run organize_files.py under scripts/preprocessing to sort the files into their respective cohorts. Doing this should create another directory underneath Data/ImmuneCODE called repertoires_org.
 
 This steps should prepare and format the data in a way that the other analysis scripts in the repository can use.
 
## DeepTCR & other dependencies

DeepTCR v2.0.10 with python 3.8 was used to run the analyses provided in this repository. All other dependencies to run analyses can be found under requirements.txt.

